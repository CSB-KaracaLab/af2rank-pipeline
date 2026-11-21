"""Step 1: candidate search + metadata download (RCSB search + GraphQL).

Writes cache files into $XRAY_SCRATCH (default: ./cache):
  candidates.json  -- assembly identifiers returned by the search query
  meta_*.json      -- GraphQL metadata per batch of entries
"""
import json, os, sys, concurrent.futures as cf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rcsb_util as u

SCR = os.environ.get("XRAY_SCRATCH", os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache"))
os.makedirs(SCR, exist_ok=True)
CUTOFF = "2021-09-30"


def t(attr, op, val):
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": attr, "operator": op, "value": val}}


SEARCH_QUERY = {
    "query": {"type": "group", "logical_operator": "and", "nodes": [
        t("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
        t("rcsb_entry_info.resolution_combined", "less_or_equal", 2.5),
        t("rcsb_accession_info.initial_release_date", "greater", CUTOFF + "T23:59:59Z"),
        t("rcsb_assembly_info.polymer_entity_instance_count", "equals", 2),
        t("rcsb_assembly_info.polymer_entity_instance_count_protein", "equals", 2),
        t("rcsb_entry_info.polymer_entity_count_nucleic_acid", "equals", 0),
        t("rcsb_assembly_container_identifiers.assembly_id", "exact_match", "1"),
    ]},
    "return_type": "assembly",
}

GQ = '''query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id
 rcsb_accession_info{initial_release_date}
 rcsb_entry_info{resolution_combined polymer_entity_count_nucleic_acid branched_entity_count}
 exptl{method}
 assemblies{ rcsb_assembly_container_identifiers{assembly_id}
   pdbx_struct_assembly{details oligomeric_count}
   pdbx_struct_assembly_gen{asym_id_list oper_expression}
   rcsb_assembly_info{polymer_entity_instance_count polymer_entity_instance_count_protein polymer_entity_count_protein}
   interfaces{ rcsb_interface_info{interface_area num_interface_residues}
               rcsb_interface_partner{interface_partner_identifier{asym_id entity_id}} } }
 polymer_entities{ rcsb_id
   entity_poly{rcsb_mutation_count rcsb_non_std_monomers rcsb_entity_polymer_type rcsb_sample_sequence_length pdbx_seq_one_letter_code_can}
   rcsb_polymer_entity{pdbx_mutation pdbx_description}
   rcsb_polymer_entity_container_identifiers{entity_id uniprot_ids reference_sequence_identifiers{database_name database_accession}}
   rcsb_cluster_membership{cluster_id identity}
   polymer_entity_instances{ rcsb_polymer_entity_instance_container_identifiers{asym_id auth_asym_id}
     rcsb_polymer_instance_info{modeled_residue_count}
     pdbx_modification_feature{type category label_comp_id modified_residue_label_comp_id modified_residue_label_asym_id label_asym_id}
     rcsb_polymer_struct_conn{connect_type role connect_partner{label_asym_id label_comp_id symmetry} connect_target{label_asym_id label_comp_id symmetry}} } } } }'''


def main():
    cfile = os.path.join(SCR, "candidates.json")
    if not os.path.exists(cfile):
        ids = u.search_all(SEARCH_QUERY)
        json.dump({"query": SEARCH_QUERY, "ids": ids}, open(cfile, "w"))
    ids = json.load(open(cfile))["ids"]
    entries = sorted({i.split("-")[0] for i in ids})
    print("assemblies:", len(ids), "entries:", len(entries))
    B = 50
    batches = [(k, entries[k:k + B]) for k in range(0, len(entries), B)]

    def run(kb):
        k, b = kb
        f = os.path.join(SCR, f"meta_{k:06d}.json")
        if os.path.exists(f):
            return
        d = u.graphql(GQ, {"ids": b})
        json.dump(d["entries"], open(f, "w"))

    with cf.ThreadPoolExecutor(6) as ex:
        for n, _ in enumerate(ex.map(run, batches)):
            if n % 20 == 0:
                print("batch", n, "/", len(batches), flush=True)
    print("done")


if __name__ == "__main__":
    main()
