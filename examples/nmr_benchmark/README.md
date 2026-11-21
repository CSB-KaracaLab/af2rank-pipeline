# Curated NMR dimer benchmark: one conformer per ensemble

AF2Rank scores for 67 protein dimers solved by solution NMR, each scored with all 17 template-capable AlphaFold2 weights. The set is a curated subset of the 99-dimer NMR benchmark of B. Savaş's PhD thesis (chapter 4), whose interfaces AlphaFold2 did not see in training. The curation notes are in [`CURATION.md`](CURATION.md).

**One conformer only.** Every score here comes from **model 1 (the first conformer) of each NMR ensemble**. A run over every deposited conformer is in progress, and its results will be added and curated the same way.

## How the dimers were scored

- **Inputs:** model 1 of each deposited ensemble, cleaned by the pipeline's strict cleaning. Sequences come from SEQRES, trimmed to the modelled span; the FASTA files are in [`fasta/`](fasta).
- **Settings:** full resolution (template sequence and side chains kept), 0 recycles, one iteration. This was the best-performing condition in the thesis benchmark.
- **Weights:** AF2-ptm models 1–2 and AF2.1, AF2.2 and AF2.3 multimer models 1–5. For multimer weights, the two chains are joined by a 20-glycine linker that is masked in the template and removed before scoring (see the [report](../../docs/af2rank_all_weights.md)).
- **Scores:** `results/scores_model1_all_weights.csv` holds pLDDT, pTM, ipTM, actifpTM, DockQ input→output and interface RMSD input→output for each dimer and weight (1,139 rows). DockQ input→output compares the input model with AlphaFold's output. It measures how far AlphaFold moves the interface, not agreement with another experimental structure.
- **Hardware:** MareNostrum 5, H100 GPUs, with the Slurm scripts in [`../marenostrum5`](../marenostrum5).

## Sub-benchmarks

| Sub-benchmark | Dimers | What it contains | Structures |
| --- | --- | --- | --- |
| NMR-homo | 29 | symmetric homodimers without long disordered segments | [figure](figures/structures_homo.png) |
| NMR-hetero | 14 | heterodimers without long disordered segments | [figure](figures/structures_hetero.png) |
| NMR-flexible | 24 | dimers with long disordered segments (10 homo, 14 hetero) | [figure](figures/structures_flexible.png) |

The PDB IDs of each sub-benchmark are listed in `kept_homo.txt`, `kept_hetero.txt` and `kept_flexible.txt`. `targets.csv` has one row for each of the original 99 dimers, giving its status, sub-benchmark, exclusion reason and properties (interface size, chain-B spread across the ensemble, backbone RMSD, flags).

![ipTM against DockQ input→output, one row per sub-benchmark and one column per weight set](figures/iptm_vs_dockq_by_subbenchmark.png)

**ipTM against DockQ input→output** for model 1 of each dimer. Each point is one dimer scored with one weight. ρ is the Spearman correlation over all points in a panel. In the NMR-flexible row, homodimers are hollow and heterodimers filled.

**Means for each sub-benchmark** (ipTM / DockQ input→output):

| Weight set | NMR-homo | NMR-hetero | NMR-flexible |
| --- | --- | --- | --- |
| AF2-ptm | 0.70 / 0.80 | 0.71 / 0.77 | 0.64 / 0.76 |
| AF2.1 multimer | 0.57 / 0.64 | 0.62 / 0.64 | 0.53 / 0.48 |
| AF2.2 multimer | 0.69 / 0.64 | 0.74 / 0.63 | 0.65 / 0.49 |
| AF2.3 multimer | 0.68 / 0.71 | 0.72 / 0.68 | 0.62 / 0.59 |

With the multimer weights, DockQ input→output is about 0.1 lower for the NMR-flexible dimers than for the other two sets, while ipTM differs less. The heterodimer and flexible sets are small (14 and 24 dimers, five points per dimer and multimer weight set), so treat their correlations as rough estimates.

![Structures of the NMR-hetero sub-benchmark](figures/structures_hetero.png)

**NMR-hetero structures** (model 1; first chain teal, second amber). The homodimer and flexible sets are in [`figures/structures_homo.png`](figures/structures_homo.png) and [`figures/structures_flexible.png`](figures/structures_flexible.png).

## Curation in brief

The thesis selected 99 NMR dimers released between March 2016 and December 2025. Each has an interface that AlphaFold2 could not have seen, no mutations, post-translational modifications or covalent interfaces, and secondary structure at the interface. From these, 32 were removed:

- **No published paper** (RCSB primary citation "To be published"): 7 dimers.
- **Not an experimentally determined complex:** HADDOCK docking or complexes assembled from earlier structures (8 dimers), and one complex built from lanthanide pseudocontact shifts only.
- **Lipid-modified peptide at the interface:** 7NQC.
- **Peptide partner** (a chain shorter than 30 residues): 16 dimers.

[`CURATION.md`](CURATION.md) lists every exclusion with its reason, the kept IDs and the decision log.

## Files

| Path | Contents |
| --- | --- |
| `CURATION.md` | selection criteria, kept and excluded dimers, decision log |
| `targets.csv` | the original 99 dimers with status, sub-benchmark, reason and properties |
| `kept_*.txt` | PDB IDs per sub-benchmark |
| `fasta/` | trimmed FASTA files of the 67 kept dimers (model coordinates: model 1 of each PDB entry) |
| `results/scores_model1_all_weights.csv` | scores for 67 dimers × 17 weights |
| `figures/` | the scatterplot and one page of structures per sub-benchmark |
| `scripts/` | `plot_subsets_grid.py` (scatterplot), `render_nmr99.py` (PyMOL thumbnails), `make_a4_page.py` (structure pages) |

Redraw the scatterplot from this folder:

```bash
python scripts/plot_subsets_grid.py results/scores_model1_all_weights.csv targets.csv figures/iptm_vs_dockq_by_subbenchmark.png
```

The PDB archive data are [CC0](https://www.rcsb.org/pages/usage-policy). Per-target scores of the original 99 dimers under all 510 thesis conditions are at [busrasavas/AF2Rank-PPI](https://github.com/busrasavas/AF2Rank-PPI).
