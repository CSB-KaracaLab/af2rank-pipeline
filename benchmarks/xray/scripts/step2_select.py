"""Step 2: filter, novelty check, redundancy reduction and length-matched sampling.

Reads the GraphQL metadata cached by step1_fetch.py ($XRAY_SCRATCH, default ./cache),
runs the extra RCSB queries it needs (cached in the same folder), and writes
xray_dimers.csv, funnel.csv, pool_nonredundant.csv and length_match.txt next to
this script.
"""
import csv, glob, json, math, os, random, re, statistics, sys
import concurrent.futures as cf

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rcsb_util as u

SCR = os.environ.get("XRAY_SCRATCH", os.path.join(HERE, "cache"))
# target set: curated NMR benchmark, rows with status == "kept" (39 homo + 45 hetero)
# The curated NMR table used for size matching (a working file, not in this repository); set NMR_CURATED_CSV to its path.
NMR_CSV = os.environ.get("NMR_CURATED_CSV", "nmr_benchmark_curated.csv")
CUTOFF = "2021-09-30T23:59:59Z"
SEED = 20261121
MIN_CHAIN, MIN_TOT, MAX_TOT = 30, 69, 900
MIN_SHORT_HETERO = 30  # v3 admitted peptides (8); v4 is back to >= 30 for every chain

# pdbx_modification_feature categories/types that are tolerated on the dimer chains.
# Everything else (phosphorylation, glycosylation, acetylation, lipidation, covalent
# ligands/cofactors, chromophores, non-standard residues/linkages, ...) excludes the entry.
OK_MOD_CATEGORIES = {"Disulfide bridge", "Metal coordination"}  # interchain disulfides are caught separately
OK_MOD_TYPES = {"Selenomethionine"}
# source organisms that mark a chain as a synthetic construct / camelid nanobody (flag column only)
FLAG_ORGANISMS = ["synthetic construct", "Lama glama", "Vicugna pacos", "Camelus dromedarius", "Camelus bactrianus"]
COVALENT_CONN = {"covalent bond", "disulfide bridge"}


def t(attr, op, val):
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": attr, "operator": op, "value": val}}


ACC_ATTR = "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
XRAY_PRE = [t("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
            t("rcsb_accession_info.initial_release_date", "less_or_equal", CUTOFF)]
ANY_PRE = [t("rcsb_accession_info.initial_release_date", "less_or_equal", CUTOFF)]


def q_entries(nodes):
    return {"query": {"type": "group", "logical_operator": "and", "nodes": nodes},
            "return_type": "entry"}


# ----------------------------------------------------------------- caching helpers
def cached(name, fn):
    f = os.path.join(SCR, name)
    if os.path.exists(f):
        return json.load(open(f))
    v = fn()
    json.dump(v, open(f, "w"))
    return v


def cached_dict(name, keys, fn):
    """Dict cache that only queries keys not seen before."""
    f = os.path.join(SCR, name)
    have = json.load(open(f)) if os.path.exists(f) else {}
    missing = [k for k in keys if k not in have]
    if missing:
        have.update(fn(missing))
        json.dump(have, open(f, "w"))
    return have


def pmap(fn, items, workers=6):
    with cf.ThreadPoolExecutor(workers) as ex:
        return list(ex.map(fn, items))


def n_ops(expr):
    """Number of operators in a pdbx_struct_assembly_gen.oper_expression."""
    groups = re.findall(r"\(([^)]*)\)", expr) or [expr]
    total = 1
    for g in groups:
        n = 0
        for item in g.split(","):
            item = item.strip()
            if "-" in item:
                a, b = item.split("-")
                n += int(b) - int(a) + 1
            elif item:
                n += 1
        total *= n
    return total


# ----------------------------------------------------------------- parse one entry
def parse_entry(e):
    """Return (dimer dict, reject reason or None)."""
    asm = next(a for a in e["assemblies"] if a["rcsb_assembly_container_identifiers"]["assembly_id"] == "1")
    inst = {}   # label asym -> (entity record, instance record)
    for pe in e["polymer_entities"]:
        for i in pe["polymer_entity_instances"]:
            inst[i["rcsb_polymer_entity_instance_container_identifiers"]["asym_id"]] = (pe, i)
    chains = []  # (asym, oper_expression)
    for g in asm["pdbx_struct_assembly_gen"]:
        k = n_ops(g["oper_expression"])
        for a in g["asym_id_list"]:
            if a in inst:
                chains += [(a, g["oper_expression"])] * k
    d = {"pdb_id": e["rcsb_id"], "assembly_id": "1",
         "resolution": min(e["rcsb_entry_info"]["resolution_combined"]),
         "release_date": e["rcsb_accession_info"]["initial_release_date"][:10],
         "assembly_details": asm["pdbx_struct_assembly"]["details"]}
    if len(chains) != 2:
        return d, "assembly1_not_2_chains"
    (a1, op1), (a2, op2) = chains
    pe1, i1 = inst[a1]
    pe2, i2 = inst[a2]
    if any(p["entity_poly"]["rcsb_entity_polymer_type"] != "Protein" for p in (pe1, pe2)):
        return d, "assembly1_not_2_chains"
    eid1 = pe1["rcsb_polymer_entity_container_identifiers"]["entity_id"]
    eid2 = pe2["rcsb_polymer_entity_container_identifiers"]["entity_id"]
    seq1 = pe1["entity_poly"]["pdbx_seq_one_letter_code_can"]
    seq2 = pe2["entity_poly"]["pdbx_seq_one_letter_code_can"]
    homo = eid1 == eid2 or seq1 == seq2
    sym = a1 == a2
    def acc(pe):
        ids = pe["rcsb_polymer_entity_container_identifiers"].get("uniprot_ids") or []
        if not ids:
            ids = [r["database_accession"] for r in
                   pe["rcsb_polymer_entity_container_identifiers"].get("reference_sequence_identifiers") or []
                   if r["database_name"] == "UniProt"]
        return sorted(set(ids))
    auth = lambda i: i["rcsb_polymer_entity_instance_container_identifiers"]["auth_asym_id"]
    d.update({
        "type": "homo" if homo else "hetero",
        "symmetry_mate": sym,
        "label_asym": [a1, a2], "opers": [op1, op2],
        "chain_auth": [auth(i1), auth(i2)],
        "entity_ids": [eid1, eid2],
        "entity_rcsb": [pe1["rcsb_id"], pe2["rcsb_id"]],
        "uniprot": [acc(pe1), acc(pe2)],
        "seq": [seq1, seq2],
        "description": [pe1["rcsb_polymer_entity"]["pdbx_description"], pe2["rcsb_polymer_entity"]["pdbx_description"]],
        "res_chain": [i1["rcsb_polymer_instance_info"]["modeled_residue_count"] or 0,
                      i2["rcsb_polymer_instance_info"]["modeled_residue_count"] or 0],
    })
    d["total"] = sum(d["res_chain"])
    iface = [x["rcsb_interface_info"]["interface_area"] for x in asm.get("interfaces") or []
             if x["rcsb_interface_info"] and x["rcsb_interface_info"].get("interface_area")]
    d["interface_area"] = round(max(iface), 1) if iface else ""

    # -- mutations
    for pe in (pe1, pe2):
        if (pe["entity_poly"]["rcsb_mutation_count"] or 0) > 0 or (pe["rcsb_polymer_entity"]["pdbx_mutation"] or "").strip() not in ("", "?", "."):
            return d, "mutation"
    # -- modified residues / PTMs
    for pe in (pe1, pe2):
        if set(pe["entity_poly"]["rcsb_non_std_monomers"] or []) - {"MSE"}:
            return d, "modified_residue_or_ptm"
    for i in (i1, i2):
        for m in i["pdbx_modification_feature"] or []:
            if m["type"] in OK_MOD_TYPES or m["category"] in OK_MOD_CATEGORIES:
                continue
            return d, "modified_residue_or_ptm"
    # -- covalent interchain links (direct, or through a ligand bonded to both chains)
    poly_asyms = set(inst)
    lig = [set(), set()]
    for k, (a, i) in enumerate(((a1, i1), (a2, i2))):
        other = a2 if k == 0 else a1
        for c in i["rcsb_polymer_struct_conn"] or []:
            if c["connect_type"] not in COVALENT_CONN:
                continue
            ends = [c["connect_partner"], c["connect_target"]]
            pa = c["connect_partner"]["label_asym_id"]
            ta = c["connect_target"]["label_asym_id"]
            if not sym and {pa, ta} == {a1, a2}:
                return d, "covalent_interchain_link"
            if sym and pa == ta == a1 and any((x.get("symmetry") or "1_555") != "1_555" for x in ends):
                return d, "covalent_interchain_link"
            for x in (pa, ta):
                if x not in poly_asyms:
                    lig[k].add(x)
    if not sym and lig[0] & lig[1]:
        return d, "covalent_interchain_link"
    # -- size
    lo, hi = sorted(d["res_chain"])
    if lo < (MIN_CHAIN if homo else MIN_SHORT_HETERO) or hi < MIN_CHAIN or not (MIN_TOT <= d["total"] <= MAX_TOT):
        return d, "size"
    # -- the two chains must touch: RCSB computed a protein-protein interface for assembly 1
    if not iface:
        return d, "no_interface"
    return d, None


# ----------------------------------------------------------------- RCSB lookups
GQ_GROUPS = '''query($ids:[String!]!){ polymer_entities(entity_ids:$ids){ rcsb_id
  rcsb_polymer_entity_group_membership{group_id aggregation_method similarity_cutoff} } }'''


def fetch_groups(entity_ids):
    out = {}
    ids = sorted(entity_ids)
    batches = [ids[k:k + 200] for k in range(0, len(ids), 200)]
    for res in pmap(lambda b: u.graphql(GQ_GROUPS, {"ids": b})["polymer_entities"], batches):
        for pe in res:
            g = {}
            for m in pe["rcsb_polymer_entity_group_membership"] or []:
                if m["aggregation_method"] == "sequence_identity":
                    g[str(int(m["similarity_cutoff"]))] = m["group_id"]
            out[pe["rcsb_id"]] = g
    return out


GQ_ORGS = '''query($ids:[String!]!){ polymer_entities(entity_ids:$ids){ rcsb_id
  rcsb_entity_source_organism{ncbi_scientific_name} } }'''


def fetch_orgs(entity_ids):
    out = {}
    for pe in u.graphql(GQ_ORGS, {"ids": list(entity_ids)})["polymer_entities"]:
        names = sorted({o["ncbi_scientific_name"] or "" for o in pe["rcsb_entity_source_organism"] or []})
        out[pe["rcsb_id"]] = ";".join(n for n in names if n)
    return out


GQ_CIT = '''query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id
  rcsb_primary_citation{journal_abbrev year pdbx_database_id_DOI} } }'''


def fetch_citations(entry_ids):
    ids = sorted(entry_ids)
    out = {}
    batches = [ids[k:k + 100] for k in range(0, len(ids), 100)]
    for res in pmap(lambda b: u.graphql(GQ_CIT, {"ids": b})["entries"], batches):
        for e in res:
            c = e["rcsb_primary_citation"] or {}
            out[e["rcsb_id"]] = {"journal": c.get("journal_abbrev") or "", "year": c.get("year") or "",
                                 "doi": c.get("pdbx_database_id_DOI") or ""}
    return out


def no_paper(c):
    """'To be published', or no journal and no DOI."""
    j = (c.get("journal") or "").lower()
    return "to be published" in j or (not j and not c.get("doi"))


def count(nodes):
    return u.search_count(q_entries(nodes))


def ident_node(ident):
    kind, val = ident.split(":", 1)
    if kind == "acc":
        return t(ACC_ATTR, "exact_match", val)
    if kind == "g100":
        return t("rcsb_polymer_entity_group_membership.group_id", "exact_match", val)
    if kind == "motif":
        return {"type": "terminal", "service": "seqmotif",
                "parameters": {"value": val, "pattern_type": "simple", "sequence_type": "protein"}}
    raise ValueError(ident)


def ident_count(ident, pre):
    """Number of pre-cutoff entries (filter list `pre`) containing this chain identifier."""
    return count([ident_node(ident)] + pre)


def together(key, pre):
    """key = 'acc:P12345|g100:678_100' -> number of pre-cutoff entries containing both."""
    return count([ident_node(x) for x in key.split("|")] + pre)


def pre_any_acc(acc):
    return count([t(ACC_ATTR, "exact_match", acc)] + ANY_PRE)


# ----------------------------------------------------------------- matching
def nmr_lengths():
    r = [x for x in csv.DictReader(open(NMR_CSV)) if x.get("status", "kept") == "kept"]
    L = {"homo": [], "hetero": []}
    for x in r:
        L["homo" if x["dimer"] == "homodimer" else "hetero"].append(sum(map(int, x["residues_scored"].split("+"))))
    return L


def nmr_short_chains():
    r = [x for x in csv.DictReader(open(NMR_CSV)) if x.get("status", "kept") == "kept"]
    return [min(map(int, x["residues_scored"].split("+"))) for x in r if x["dimer"] == "heterodimer"]


def nmr_hetero_targets():
    """(total, shorter chain) per NMR heterodimer."""
    r = list(csv.DictReader(open(NMR_CSV)))
    out = []
    for x in r:
        if x["dimer"] == "heterodimer":
            v = list(map(int, x["residues_scored"].split("+")))
            out.append((sum(v), min(v)))
    return out


def cost_total(p, tgt):
    return abs(math.log(p["total"]) - math.log(tgt))


def cost_total_short(p, tgt):
    """|dlog total| + |dlog shorter chain|; tgt = (total, shorter chain)."""
    return abs(math.log(p["total"]) - math.log(tgt[0])) + abs(math.log(min(p["res_chain"])) - math.log(tgt[1]))


def match(pool, targets, rng, cost=cost_total):
    """Optimal 1:1 assignment of pool dimers to NMR targets minimising `cost` (tie-break random)."""
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    pool = pool[:]
    rng.shuffle(pool)
    if len(pool) <= len(targets):
        return pool
    jitter = np.array([rng.random() * 1e-6 for _ in pool])
    C = np.array([[cost(p, tg) for p in pool] for tg in targets]) + jitter[None, :]
    r, c = linear_sum_assignment(C)
    return [pool[j] for j in c]


def ks(a, b):
    xs = sorted(set(a) | set(b))
    F = lambda s, x: sum(v <= x for v in s) / len(s)
    return max(abs(F(a, x) - F(b, x)) for x in xs)


def summary(L):
    L = sorted(L)
    q = statistics.quantiles(L, n=4) if len(L) > 1 else [L[0]] * 3
    return f"n={len(L)} min={L[0]} Q1={q[0]:.0f} median={statistics.median(L):.0f} Q3={q[2]:.0f} max={L[-1]}"


# ----------------------------------------------------------------- main
def main():
    entries = []
    for f in sorted(glob.glob(os.path.join(SCR, "meta_*.json"))):
        entries += json.load(open(f))
    cand = json.load(open(os.path.join(SCR, "candidates.json")))["ids"]
    funnel = [("search: X-ray, res<=2.5 A, released after 2021-09-30, assembly 1 = 2 protein chains, no NA polymer in entry", len(cand))]

    reasons = {}
    passed = []
    for e in entries:
        d, why = parse_entry(e)
        reasons[why] = reasons.get(why, 0) + 1
        if why is None:
            passed.append(d)
    n = len(entries)
    for step in ["assembly1_not_2_chains", "mutation", "modified_residue_or_ptm", "covalent_interchain_link", "size", "no_interface"]:
        n -= reasons.get(step, 0)
        label = {"assembly1_not_2_chains": "assembly 1 = one entity x2 or two entities x1 (protein)",
                 "mutation": "no engineered mutation (rcsb_mutation_count=0, no pdbx_mutation)",
                 "modified_residue_or_ptm": "no non-standard residue except MSE, no PTM/covalent modification feature",
                 "covalent_interchain_link": "no covalent/disulfide link between the two chains",
                 "size": f"homo: each chain >= {MIN_CHAIN}; hetero: longer chain >= {MIN_CHAIN}, shorter >= {MIN_SHORT_HETERO} modeled residues; total {MIN_TOT}-{MAX_TOT}",
                 "no_interface": "RCSB reports an interface between the two chains in assembly 1"}[step]
        funnel.append((label, n))
    print("entries parsed:", len(entries), reasons)

    # ---- group memberships (30% / 100% sequence clusters)
    ents = {x for d in passed for x in d["entity_rcsb"]}
    groups = cached_dict("groups.json", sorted(ents), fetch_groups)

    # ---- "hetero" entries whose two entities are the same protein (same UniProt accession):
    #      same 100% cluster and lengths within 20%  -> construct variants, treat as homodimer;
    #      otherwise (e.g. serpin + cleaved RCL peptide, protease + propeptide) -> fragments of one
    #      polypeptide, not a two-protein dimer -> excluded.
    kept = []
    for d in passed:
        if d["type"] == "hetero" and set(d["uniprot"][0]) & set(d["uniprot"][1]):
            g0 = groups.get(d["entity_rcsb"][0], {}).get("100")
            g1 = groups.get(d["entity_rcsb"][1], {}).get("100")
            a, b = sorted(d["res_chain"])
            if g0 and g0 == g1 and a >= 0.8 * b:
                d["type"] = "homo"
                d["homo_note"] = "two entities, same UniProt and 100% cluster"
            else:
                continue
        kept.append(d)
    passed = kept
    funnel.append(("two distinct proteins for heterodimers (same-UniProt fragment pairs dropped; same-construct pairs -> homo)", len(passed)))

    # ---- novelty
    # A chain is identified by its UniProt accession(s) and its RCSB 100% sequence cluster. Peptides
    # too short to be clustered (< 10 residues) and without UniProt fall back to an exact sequence-motif
    # search (substring of any entity), which is stricter than 100% identity.
    def idents(d, k):
        ids = [("acc", a) for a in d["uniprot"][k]]
        g = groups.get(d["entity_rcsb"][k], {}).get("100")
        if g:
            ids.append(("g100", g))
        if not g and not d["uniprot"][k]:
            ids.append(("motif", d["seq"][k]))
        return ids

    all_ids = sorted({i for d in passed for k in (0, 1) for i in idents(d, k)})
    xr = cached_dict("ident_pre_xray.json", [f"{k}:{v}" for k, v in all_ids],
                     lambda ks: dict(zip(ks, pmap(lambda x: ident_count(x, XRAY_PRE), ks))))

    def seen(d, k, table):
        hits = [f"{a}:{b}" for a, b in idents(d, k) if table[f"{a}:{b}"] > 0]
        return bool(hits), hits

    # any-method counts only needed for chains that are new by X-ray
    new_ids = sorted({f"{a}:{b}" for d in passed for k in (0, 1) if not seen(d, k, xr)[0] for a, b in idents(d, k)})
    anym = cached_dict("ident_pre_any.json", new_ids,
                       lambda ks: dict(zip(ks, pmap(lambda x: ident_count(x, ANY_PRE), ks))))
    anym_full = {**xr, **anym}   # X-ray hits are also any-method hits

    def pair_keys(d):
        ia, ib = idents(d, 0), idents(d, 1)
        return ia, ib, [f"{ka}:{va}|{kb}:{vb}" for ka, va in ia for kb, vb in ib]

    both_x = [d for d in passed if d["type"] == "hetero" and seen(d, 0, xr)[0] and seen(d, 1, xr)[0]]
    pk = sorted({k for d in both_x for k in pair_keys(d)[2] if k.split("|")[0] != k.split("|")[1]})
    tog_x = cached_dict("pair_together_xray.json", pk, lambda ks: dict(zip(ks, pmap(lambda x: together(x, XRAY_PRE), ks))))
    tog_a = cached_dict("pair_together_any.json", pk, lambda ks: dict(zip(ks, pmap(lambda x: together(x, ANY_PRE), ks))))

    def basis(d, k):
        return ", ".join(f"{a}:{b}" for a, b in idents(d, k))

    novel_x, novel = [], []
    for d in passed:
        x0, x1 = seen(d, 0, xr)[0], seen(d, 1, xr)[0]
        a0, a1 = seen(d, 0, anym_full)[0], seen(d, 1, anym_full)[0]
        d["novelty_basis"] = basis(d, 0) if d["type"] == "homo" else f"A: {basis(d, 0)}; B: {basis(d, 1)}"
        ia, ib, keys = pair_keys(d)
        if d["type"] == "homo":
            if x0 or x1:
                continue
            novel_x.append(d)
            if a0 or a1:
                continue
            d["novelty"] = "no pre-cutoff structure (any method)"
            novel.append(d)
            continue
        # -- step 1: X-ray rule (as in v1)
        if x0 and x1:
            if set(ia) & set(ib) or any(tog_x[k] > 0 for k in keys):
                continue
        novel_x.append(d)
        # -- step 2: a chain that is new by X-ray must have no earlier entry of any method,
        #            and "never together" is checked across all methods
        if (not x0 and a0) or (not x1 and a1):
            continue
        if x0 and x1:
            if any(tog_a[k] > 0 for k in keys):
                continue
            d["novelty"] = "both pre-cutoff, never together (any method)"
        elif not x0 and not x1:
            d["novelty"] = "both partners: no pre-cutoff structure (any method)"
        elif not x1:
            d["novelty"] = "partner B: no pre-cutoff structure (any method)"
        else:
            d["novelty"] = "partner A: no pre-cutoff structure (any method)"
        novel.append(d)
    funnel.append(("novelty by X-ray: UniProt AND 100% cluster (motif search for unclustered peptides); hetero pairs per Savas rule", len(novel_x)))
    funnel.append(("novelty across methods: new chains have no pre-cutoff entry of any method; never-together checked across methods", len(novel)))
    print("novel:", len(novel), sum(d["type"] == "homo" for d in novel), "homo")

    # ---- published only: drop entries whose primary citation is "To be published" (or has no journal and no DOI)
    cits = cached_dict("citations.json", sorted({d["pdb_id"] for d in novel}), fetch_citations)
    for d in novel:
        d["citation"] = cits.get(d["pdb_id"], {})
    novel = [d for d in novel if not no_paper(d["citation"])]
    funnel.append(("published: primary citation is not 'To be published' (and has a journal or DOI)", len(novel)))

    # ---- redundancy: one per 30% cluster (homo) / unordered pair of 30% clusters (hetero)
    def g30(d, k):
        g = groups.get(d["entity_rcsb"][k], {}).get("30")
        if g:
            return g
        # unclustered short peptide: group by UniProt accession, else by exact sequence
        return ("uniprot:" + "/".join(d["uniprot"][k])) if d["uniprot"][k] else ("seq:" + d["seq"][k])
    best = {}
    for d in novel:
        d["cluster30"] = [g30(d, 0), g30(d, 1)]
        key = (d["type"], tuple(sorted(set(d["cluster30"]))) if d["type"] == "homo" else tuple(sorted(d["cluster30"])))
        cur = best.get(key)
        rank = (d["resolution"], d["release_date"], d["pdb_id"])
        if cur is None or rank < (cur["resolution"], cur["release_date"], cur["pdb_id"]):
            best[key] = d
    nr = sorted(best.values(), key=lambda d: d["pdb_id"])
    funnel.append(("non-redundant: best resolution per 30% cluster (homo) / 30% cluster pair (hetero)", len(nr)))
    homo = [d for d in nr if d["type"] == "homo"]
    het = [d for d in nr if d["type"] == "hetero"]
    funnel.append(("  of which homodimers", len(homo)))
    funnel.append(("  of which heterodimers", len(het)))

    # ---- heterodimers: at most one entry per 30% cluster of the LONGER chain (a domain family such as
    #      PDZ bound to different peptides is kept once); best resolution, then earlier release, then PDB ID
    def long_k(d):
        a, b = d["res_chain"]
        if a != b:
            return 0 if a > b else 1
        la, lb = len(d["seq"][0]), len(d["seq"][1])
        if la != lb:
            return 0 if la > lb else 1
        return 0 if d["cluster30"][0] <= d["cluster30"][1] else 1

    def one_per_long_cluster(ds):
        best = {}
        for d in ds:
            key = d["cluster30"][long_k(d)]
            rank = (d["resolution"], d["release_date"], d["pdb_id"])
            cur = best.get(key)
            if cur is None or rank < (cur["resolution"], cur["release_date"], cur["pdb_id"]):
                best[key] = d
        return sorted(best.values(), key=lambda d: d["pdb_id"])

    nr_het_all = het
    het = one_per_long_cluster(het)
    funnel.append(("heterodimers: one per 30% cluster of the longer chain", len(het)))
    nr = sorted(homo + het, key=lambda d: d["pdb_id"])

    # ---- length-matched selection on log total length (homo and hetero separately)
    L = nmr_lengths()
    Ls = nmr_short_chains()
    rng = random.Random(SEED)
    sel_homo = match(homo, L["homo"], rng)
    sel_het = match(het, L["hetero"], rng)
    sel = sorted(sel_homo + sel_het, key=lambda d: (d["type"] != "homo", d["total"], d["pdb_id"]))
    # variant: both chains must have a UniProt mapping (drops synthetic binders / designed proteins);
    # the long-chain rule is re-applied after the UniProt restriction
    uni = lambda d: all(d["uniprot"])
    homo_u = [d for d in homo if uni(d)]
    het_u = one_per_long_cluster([d for d in nr_het_all if uni(d)])
    rng_u = random.Random(SEED)
    sel_u = sorted(match(homo_u, L["homo"], rng_u) + match(het_u, L["hetero"], rng_u),
                   key=lambda d: (d["type"] != "homo", d["total"], d["pdb_id"]))
    funnel.append(("[variant] pool with UniProt on both chains: homo", len(homo_u)))
    funnel.append(("[variant] pool with UniProt on both chains: hetero (one per long-chain 30% cluster)", len(het_u)))
    funnel.append(("length-matched selection (homo)", len(sel_homo)))
    funnel.append(("length-matched selection (hetero)", len(sel_het)))

    # ---- informational: pre-cutoff entries of ANY method per UniProt accession of the selected chains
    sel_accs = sorted({a for d in sel + sel_u for ul in d["uniprot"] for a in ul})
    any_pre = cached_dict("acc_pre_any.json", sel_accs, lambda ks: dict(zip(ks, pmap(pre_any_acc, ks))))
    # source organism of the selected entities -> synthetic construct / nanobody flag
    sel_ents = sorted({x for d in sel + sel_u for x in d["entity_rcsb"]})
    orgs = cached_dict("orgs.json", sel_ents, fetch_orgs)

    def flag(d):
        why = []
        for k, lab in ((0, "A"), (1, "B")):
            o = orgs.get(d["entity_rcsb"][k], "")
            hit = [n for n in FLAG_ORGANISMS if n in o.split(";")]
            if hit:
                why.append(f"{lab}: {hit[0]}")
            elif not d["uniprot"][k]:
                why.append(f"{lab}: no UniProt")
        return "; ".join(why)

    # ---- outputs
    cols = ["pdb_id", "assembly_id", "type", "chain_A_auth", "chain_B_auth", "symmetry_mate", "entity_A", "entity_B",
            "uniprot_A", "uniprot_B", "residues_A", "residues_B", "residues_total", "resolution_A", "release_date",
            "novelty", "novelty_basis", "cluster30_A", "cluster30_B", "cluster100_A", "cluster100_B",
            "pre_cutoff_entries_any_method_A", "pre_cutoff_entries_any_method_B",
            "interface_area_A2", "assembly_details", "note", "description_A", "description_B",
            "organism_A", "organism_B", "synthetic_or_nanobody",
            "no_paper", "citation_journal", "citation_year", "citation_doi"]

    def row(d, full=True):
        anyc = [sum(any_pre.get(a, 0) for a in ul) if ul else "" for ul in d["uniprot"]]
        return {"pdb_id": d["pdb_id"], "assembly_id": d["assembly_id"], "type": d["type"],
                "chain_A_auth": d["chain_auth"][0], "chain_B_auth": d["chain_auth"][1],
                "symmetry_mate": "yes" if d["symmetry_mate"] else "no",
                "entity_A": d["entity_ids"][0], "entity_B": d["entity_ids"][1],
                "uniprot_A": ";".join(d["uniprot"][0]), "uniprot_B": ";".join(d["uniprot"][1]),
                "residues_A": d["res_chain"][0], "residues_B": d["res_chain"][1], "residues_total": d["total"],
                "resolution_A": round(d["resolution"], 2), "release_date": d["release_date"],
                "novelty": d.get("novelty", ""), "novelty_basis": d.get("novelty_basis", ""),
                "cluster30_A": d["cluster30"][0], "cluster30_B": d["cluster30"][1],
                "cluster100_A": groups.get(d["entity_rcsb"][0], {}).get("100", ""),
                "cluster100_B": groups.get(d["entity_rcsb"][1], {}).get("100", ""),
                "pre_cutoff_entries_any_method_A": anyc[0], "pre_cutoff_entries_any_method_B": anyc[1],
                "interface_area_A2": d["interface_area"], "assembly_details": d["assembly_details"],
                "note": d.get("homo_note", ""),
                "description_A": d["description"][0], "description_B": d["description"][1],
                "organism_A": orgs.get(d["entity_rcsb"][0], "") if full else "",
                "organism_B": orgs.get(d["entity_rcsb"][1], "") if full else "",
                "synthetic_or_nanobody": flag(d) if full else "",
                "no_paper": no_paper(d.get("citation", {})),
                "citation_journal": d.get("citation", {}).get("journal", ""),
                "citation_year": d.get("citation", {}).get("year", ""),
                "citation_doi": d.get("citation", {}).get("doi", "")}

    with open(os.path.join(HERE, "xray_dimers.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, cols)
        w.writeheader()
        for d in sel:
            w.writerow(row(d))
    with open(os.path.join(HERE, "xray_dimers_uniprot_only.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, cols)
        w.writeheader()
        for d in sel_u:
            w.writerow(row(d))
    with open(os.path.join(HERE, "pool_nonredundant.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, [c for c in cols if not c.startswith(("pre_cutoff_entries_any", "organism", "synthetic"))], extrasaction="ignore")
        w.writeheader()
        for d in nr:   # organisms / any-method counts / flag are only fetched for the selected sets
            w.writerow(row(d, full=False))
    with open(os.path.join(HERE, "funnel.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["step", "count"])
        w.writerows(funnel)

    lines = [f"target: {os.path.basename(NMR_CSV)} (status == kept): {len(L['homo'])} homo, {len(L['hetero'])} hetero", ""]
    for typ, s in (("homo", sel_homo), ("hetero", sel_het)):
        a = [d["total"] for d in s]
        lines.append(f"{typ:6s} NMR   : {summary(L[typ])}")
        lines.append(f"{typ:6s} X-ray : {summary(a)}   KS D={ks(a, L[typ]):.3f}")
    allx = [d["total"] for d in sel]
    alln = L["homo"] + L["hetero"]
    lines.append(f"all    NMR   : {summary(alln)}")
    lines.append(f"all    X-ray : {summary(allx)}   KS D={ks(allx, alln):.3f}")
    lines.append("")
    lines.append("hetero shorter-chain length (modeled residues)")
    lines.append(f"  NMR   : {summary(Ls)}")
    a = [min(d["res_chain"]) for d in sel_het]
    lines.append(f"  X-ray : {summary(a)}   KS D={ks(a, Ls):.3f}")
    lines.append(f"  shorter chain < 30: NMR {sum(v < 30 for v in Ls)}/{len(Ls)}, X-ray {sum(v < 30 for v in a)}/{len(a)}")
    lines.append("")
    lines.append("[variant: UniProt on both chains, xray_dimers_uniprot_only.csv]")
    for typ in ("homo", "hetero"):
        a = [d["total"] for d in sel_u if d["type"] == typ]
        lines.append(f"{typ:6s} X-ray : {summary(a)}   KS D={ks(a, L[typ]):.3f}")
    a = [min(d["res_chain"]) for d in sel_u if d["type"] == "hetero"]
    lines.append(f"hetero shorter chain X-ray : {summary(a)}   KS D={ks(a, Ls):.3f}   shorter chain < 30: {sum(v < 30 for v in a)}/{len(a)}")
    lines.append("")
    lines.append(f"selected entries flagged synthetic_or_nanobody: main {sum(bool(flag(d)) for d in sel)}/{len(sel)} "
                 f"(homo {sum(bool(flag(d)) for d in sel_homo)}, hetero {sum(bool(flag(d)) for d in sel_het)}); "
                 f"uniprot_only {sum(bool(flag(d)) for d in sel_u)}/{len(sel_u)}")
    open(os.path.join(HERE, "length_match.txt"), "w").write("\n".join(lines) + "\n")
    for f in funnel:
        print(f"{f[1]:7d}  {f[0]}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
