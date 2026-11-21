# Trypsin–BPTI example (PDB 2PTC)

`model.pdb` is [PDB 2PTC](https://doi.org/10.2210/pdb2PTC/pdb) as deposited (Marquart et al., 1983): bovine β-trypsin (chain E) bound to the bovine pancreatic trypsin inhibitor (chain I), with waters and a calcium ion that the pipeline's cleaning removes. `target.fasta` holds the SEQRES sequences (223 and 58 residues).

```bash
af2rank-pipeline run-target --target 2ptc --fasta examples/2ptc/target.fasta \
  --models examples/2ptc/model.pdb --params params --strict-clean --out runs/2ptc
```

`results/scores_all_weights.csv` lists the scores for all 17 template-capable AF2 weights (multimer weights with the 20-glycine linker), from a run of the Colab notebook's cells on a laptop CPU (71 min in total). Figure 1 of the [report](../../docs/af2rank_all_weights.md) plots them.

The PDB archive data are [CC0](https://www.rcsb.org/pages/usage-policy).
