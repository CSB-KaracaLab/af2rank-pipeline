#!/usr/bin/env python3
"""Build AF2Rank inputs for the X-ray dimer benchmark in the nmr99_inputs layout.

For each row of ../xray_dimers.csv, read ../raw/<ID>-assembly<N>.cif (untrusted data, parsed only) and write
  <out>/<id>/target.fasta            >ID_A / >ID_B, canonical SEQRES trimmed to the modelled span
  <out>/<id>/models/<id>_model1.pdb  protein ATOM records, chains A/B, heavy atoms, first altloc, MSE -> MET
and <out>/prep_report.json.

Usage: python -I prep_targets.py <xray_dimers.csv> <raw dir> <out dir>
Needs gemmi.
"""
from __future__ import annotations
import csv
import json
import sys
from pathlib import Path

import gemmi

STD = "ACDEFGHIKLMNPQRSTVWY"
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
         "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
         "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def parent_letter(mon_id: str) -> str | None:
    info = gemmi.find_tabulated_residue(mon_id)
    if info is None:
        return None
    code = info.one_letter_code.upper()
    return code if code in STD else None


def atom_name_field(name: str, element: str) -> str:
    if len(name) >= 4:
        return name[:4]
    return (" " + name).ljust(4) if len(element) == 1 else name.ljust(4)


def pdb_line(serial, name, resname, chain, seqnum, icode, pos, occ, b, element):
    return (f"ATOM  {serial:5d} {atom_name_field(name, element)} {resname:>3s} {chain}{seqnum:4d}{icode}   "
            f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}{occ:6.2f}{b:6.2f}          {element:>2s}  ")


def entity_tables(block):
    can, seqs = {}, {}
    for row in block.find("_entity_poly.", ["entity_id", "pdbx_seq_one_letter_code_can"]):
        can[row[0]] = "".join(gemmi.cif.as_string(row[1]).split())
    for row in block.find("_entity_poly_seq.", ["entity_id", "num", "mon_id"]):
        seqs.setdefault(row[0], {}).setdefault(int(row[1]), []).append(row[2])
    return can, seqs


def chain_residues(chain, subchain_entity, entity_id, notes):
    """Protein residues of entity_id with a CA, first altloc, MSE converted. Returns list of residue dicts."""
    out = []
    for res in chain:
        if subchain_entity.get(res.subchain) != entity_id:
            continue
        name = res.name
        if name not in THREE and name != "MSE":
            info = gemmi.find_tabulated_residue(name)
            if info is None or not info.is_amino_acid():
                continue
            notes.append(f"chain {chain.name}: non-standard residue {name} {res.seqid.num}{res.seqid.icode.strip()} kept as-is")
        if res.find_atom("CA", "*") is None:
            continue
        atoms = []
        for at in res:
            el = at.element.name.upper()
            aname = at.name
            if name == "MSE" and aname == "SE":
                aname, el = "SD", "S"
            atoms.append((aname, at.pos, at.occ, at.b_iso, el))
        out.append(dict(resname="MET" if name == "MSE" else name, src_name=name, num=res.seqid.num,
                        icode=res.seqid.icode if res.seqid.icode != "\x00" else " ", label_seq=res.label_seq,
                        atoms=atoms))
    return out


def main(csv_path, raw_dir, out_dir):
    raw_dir, out_dir = Path(raw_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report, problems = {}, {}
    rows = list(csv.DictReader(open(csv_path)))
    for r in rows:
        pid = r["pdb_id"].upper()
        asm = r["assembly_id"] or "1"
        notes = []
        path = raw_dir / f"{pid}-assembly{asm}.cif"
        doc = gemmi.cif.read(str(path))
        block = doc.sole_block()
        st = gemmi.make_structure_from_block(block)
        st.setup_entities()
        st.remove_hydrogens()
        st.remove_alternative_conformations()  # keeps the first altloc (also point microheterogeneity)
        model = st[0]
        subchain_entity = {sc: ent.name for ent in st.entities if ent.entity_type == gemmi.EntityType.Polymer
                           for sc in ent.subchains}
        can, seqs = entity_tables(block)
        ent_a, ent_b = r["entity_A"], r["entity_B"]
        poly_chains = [ch for ch in model if any(subchain_entity.get(res.subchain) in (ent_a, ent_b) for res in ch)]
        if r["symmetry_mate"] == "yes":
            assert r["chain_A_auth"] == r["chain_B_auth"]
            cands = [ch for ch in poly_chains if ch.name == r["chain_A_auth"] or ch.name.startswith(r["chain_A_auth"] + "-")]
            if len(cands) != 2:
                problems[pid] = f"symmetry mate: expected 2 copies of chain {r['chain_A_auth']}, found {[c.name for c in cands]}"
                continue
            cands.sort(key=lambda c: (c.name != r["chain_A_auth"], c.name))  # original chain first
            chA, chB = cands
        else:
            # A deposited chain may enter the assembly only as a symmetry copy, named e.g. "B-2" (7VIM, 8KDX).
            def pick(auth):
                hits = [ch for ch in poly_chains if ch.name == auth] or \
                       [ch for ch in poly_chains if ch.name.split("-")[0] == auth]
                return hits[0] if len(hits) == 1 else None
            chA, chB = pick(r["chain_A_auth"]), pick(r["chain_B_auth"])
            if chA is None or chB is None:
                problems[pid] = f"chains {r['chain_A_auth']},{r['chain_B_auth']} not both in assembly ({[c.name for c in poly_chains]})"
                continue
            if chB.name != r["chain_B_auth"] or chA.name != r["chain_A_auth"]:
                notes.append(f"assembly holds the chains as {chA.name},{chB.name} (symmetry copy of a deposited chain)")
        if len(poly_chains) != 2:
            notes.append(f"assembly has {len(poly_chains)} chains of the dimer entities: {[c.name for c in poly_chains]}")
        resA = chain_residues(chA, subchain_entity, ent_a, notes)
        resB = chain_residues(chB, subchain_entity, ent_b, notes)

        # Canonical sequence per entity, non-standard letters -> parent
        nonstd = []
        full = {}
        for ent in {ent_a, ent_b}:
            s = list(can[ent])
            if len(s) != len(seqs[ent]):
                notes.append(f"entity {ent}: can length {len(s)} != entity_poly_seq length {len(seqs[ent])}")
            for i, ch in enumerate(s):
                if ch not in STD:
                    mons = seqs[ent].get(i + 1, ["?"])
                    par = parent_letter(mons[0])
                    nonstd.append({"entity": ent, "position": i + 1, "mon_id": mons[0], "letter": ch,
                                   "replaced_by": par or "?"})
                    s[i] = par or "X"
            full[ent] = "".join(s)

        def span(res):
            ls = [x["label_seq"] for x in res if x["label_seq"] is not None]
            return min(ls), max(ls)

        sa, sb = span(resA), span(resB)
        if ent_a == ent_b:  # homodimer of one entity: one span covering both copies, identical sequences
            lo, hi = min(sa[0], sb[0]), max(sa[1], sb[1])
            sa = sb = (lo, hi)
        seqA = full[ent_a][sa[0] - 1:sa[1]]
        seqB = full[ent_b][sb[0] - 1:sb[1]]
        if "X" in seqA + seqB:
            notes.append("unresolved non-standard letter left as X")

        # Check coordinate residue names against the SEQRES letters
        for lab, res, ent in (("A", resA, ent_a), ("B", resB, ent_b)):
            mism = [f"{x['src_name']}{x['num']}{x['icode'].strip()}->{full[ent][x['label_seq'] - 1]}" for x in res
                    if THREE.get(x["resname"], "?") != full[ent][x["label_seq"] - 1]]
            if mism:
                notes.append(f"chain {lab}: {len(mism)} residue(s) differ from SEQRES: {', '.join(mism[:10])}")

        tdir = out_dir / pid.lower()
        (tdir / "models").mkdir(parents=True, exist_ok=True)
        (tdir / "target.fasta").write_text(f">{pid}_A\n{seqA}\n>{pid}_B\n{seqB}\n")
        lines, serial = [], 0
        for lab, res in (("A", resA), ("B", resB)):
            for x in res:
                for aname, pos, occ, b, el in x["atoms"]:
                    serial += 1
                    lines.append(pdb_line(serial, aname, x["resname"], lab, x["num"], x["icode"], pos, occ, b, el))
        lines.append("END")
        (tdir / "models" / f"{pid.lower()}_model1.pdb").write_text("\n".join(lines) + "\n")

        report[pid] = {
            "models": 1,
            "chains": [chA.name, chB.name],
            "lengths": [len(seqA), len(seqB)],
            "homodimer": seqA == seqB,
            "nonstandard": nonstd,
            "seqres_lengths": [len(full[ent_a]), len(full[ent_b])],
            "seqres_span": [list(sa), list(sb)],
            "modelled_residues": [len(resA), len(resB)],
            "symmetry_mate": r["symmetry_mate"] == "yes",
            "notes": notes,
        }
        print(pid, report[pid]["chains"], report[pid]["lengths"], report[pid]["modelled_residues"],
              "csv:", r["residues_A"], r["residues_B"], "; ".join(notes))
    (out_dir / "prep_report.json").write_text(json.dumps(report, indent=1) + "\n")
    print(f"{len(report)} targets written; problems: {problems or 'none'}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
