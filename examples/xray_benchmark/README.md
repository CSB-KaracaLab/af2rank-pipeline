# Unseen X-ray dimer benchmark

X-ray protein dimers whose interfaces no AlphaFold2 version could have seen in training, curated as a control for the [NMR dimer benchmark](../nmr_benchmark). The crystal structures pass the same novelty test and the same quality rules as the curated NMR set.

- **X-ray-hetero:** 53 heterodimers. 27 are matched in size to the NMR heterodimers; 26 are larger natural complexes, since X-ray can determine larger complexes than NMR.
- **X-ray-homo:** 39 homodimers, matched in size to the NMR homodimers.

The full criteria, filtering counts and decision log are in [`CURATION.md`](CURATION.md).

**Status:** AF2Rank scores (17 weights, 20-glycine linker for multimer weights) are complete for all 53 heterodimers in two settings:
- full resolution, 0 recycles ([below](#scores-x-ray-hetero));
- the AF2Rank paper's masked setting ([below](#scores-in-the-af2rank-papers-masked-setting)).

X-ray-homo has not been scored yet.

## Curation logic

1. **Unseen by every AlphaFold2 version.** Released after 30 September 2021, the PDB cutoff of AF2-multimer v3. Novelty is judged per protein, by UniProt accession plus 100% sequence cluster:
   - homodimers: the protein had no earlier structure;
   - heterodimers: at least one partner had no earlier structure, or both had earlier structures but never together in one entry.

   Earlier structures of any method (X-ray, NMR, cryo-EM) count.
2. **Same quality rules as the curated NMR set:**
   - published (the RCSB primary citation is a paper, not "To be published");
   - no engineered mutations, modified residues, glycans or covalent cofactors, and no covalent link between the chains;
   - no peptides: every chain has at least 30 modelled residues;
   - natural proteins only: no synthetic constructs, designed binders or nanobodies.
3. **Clean two-chain assemblies.** X-ray at 2.5 Å or better. Biological assembly 1 is exactly two protein chains with an interface reported by RCSB, and no nucleic acid.
4. **Non-redundant.** One entry per 30%-identity cluster (homodimers) or cluster pair (heterodimers). Heterodimers are also limited to one entry per cluster of the longer chain.
5. **Size.** Matched to the curated NMR homo- and heterodimers by optimal one-to-one matching on log total length. Larger heterodimers that pass every other rule are added without matching.

| Step | Entries |
| --- | --- |
| Search: X-ray, ≤ 2.5 Å, released after 2021-09-30, two-chain protein assembly | 13,278 |
| No mutation, modification or covalent link between the chains | 7,963 |
| Chains ≥ 30 residues, interface reported, two distinct proteins (heterodimers) | 6,537 |
| Novel across all methods | 3,050 |
| Published | 2,166 |
| Non-redundant | 907 homo + 214 hetero |
| Natural proteins | 817 homo + 121 hetero |
| Selected | 39 homo + 53 hetero (8SJJ, a metastable coiled coil, later removed) |

## Scores: X-ray-hetero

![ipTM against DockQ input→output for the 53 X-ray heterodimers](figures/iptm_vs_dockq_hetero.png)

**All 53 heterodimers** (one model each; small points: one weight, large points: the dimer mean over the weights of a set; ρ is the Spearman correlation of the dimer means).

| Weight set | Mean ipTM | Mean DockQ input→output | ρ |
| --- | --- | --- | --- |
| AF2-ptm | 0.89 | 0.97 | 0.12 |
| AF2.1 multimer | 0.85 | 0.95 | 0.19 |
| AF2.2 multimer | 0.87 | 0.93 | 0.38 |
| AF2.3 multimer | 0.87 | 0.95 | 0.48 |

AlphaFold keeps every crystal interface (DockQ input→output 0.93–0.97 on average; lowest dimer means 0.87–0.91) with high ipTM. The larger complexes behave like the size-matched ones, with slightly higher ipTM. With so little spread in DockQ, the correlation is weak except for AF2.3 (ρ 0.48), unlike the NMR heterodimers. At full resolution the template keeps its sequence and side chains, so well-packed crystal interfaces may largely be reproduced; the NMR–X-ray difference reflects structure quality as well as interface recognition. Scores: `results/scores_xray_hetero_all_weights.csv`; figure: `scripts/plot_scores.py`.

## Scores in the AF2Rank paper's masked setting

The same 53 heterodimers and 17 weights, scored as in the AF2Rank paper (Roney & Ovchinnikov, [PRL 129, 238101, 2022](https://doi.org/10.1103/PhysRevLett.129.238101)):
- template sequence masked;
- side chains masked;
- 1 recycle, i.e. `num_recycles=1` and `iterations=1`.

These are the pipeline defaults (`--recycles 1`, without `--no-mask-sequence` / `--no-mask-sidechains`), and they match the authors' [AF2Rank notebook](https://github.com/sokrypton/ColabDesign/blob/main/af/examples/AF2Rank.ipynb). The question: do natives still score high when AlphaFold cannot simply copy the template's sequence and side chains?

![Native X-ray heterodimers, full resolution vs masked](figures/natives_masked_vs_fullres.png)

Medians over dimers × weights:

| Weight set | ipTM full → masked | Runs with ipTM ≥ 0.6 (masked) | DockQ input→output full → masked | Interface lost (masked DockQ_io < 0.23) |
| --- | --- | --- | --- | --- |
| AF2-ptm | 0.89 → 0.84 | 99% | 0.97 → 0.83 | 0 of 106 |
| AF2.1 multimer | 0.86 → 0.76 | 70% | 0.96 → 0.76 | 31 of 265 |
| AF2.2 multimer | 0.87 → 0.82 | 96% | 0.95 → 0.83 | 1 of 265 |
| AF2.3 multimer | 0.87 → 0.83 | 96% | 0.96 → 0.89 | 1 of 265 |

**With masking, AlphaFold no longer just copies the crystal interface.**
- DockQ input→output drops from about 0.96 to 0.83.
- With AF2-ptm, AF2.2 and AF2.3, natives stay confident and in the acceptable range: only 8PFC is lost, once each with AF2.2 and AF2.3.
- AF2.1 is the exception: 31 of 265 runs lose the native interface, worst for 7WZ6, 6ZH1, 8PFC, 8HFP and 9WLP.

At full resolution AlphaFold reproduces almost any template interface, right or wrong. The masked setting leaves room for wrong interfaces to move while native ones stay, which is the condition for using "confident and kept" to rank docking decoys.

![Masked DockQ input→output per dimer and weight](figures/masked_dockq_io_heatmap.png)

Scores: `results/scores_xray_hetero_masked_all_weights.csv`; figures: `scripts/plot_masked.py`.

## Structures

![X-ray-hetero structures](figures/structures_hetero.png)

**X-ray-hetero** (biological assembly 1; first chain teal, second amber). X-ray-homo is in [`figures/structures_homo.png`](figures/structures_homo.png).

## Comparison with Genz et al. 2025

The heterodimer benchmark of Genz et al. (2025) has 262 heterodimers released after the AlphaFold2 cutoff ([doi:10.1002/pro.70327](https://doi.org/10.1002/pro.70327)). Its selection does not require interface novelty. [`genz2025_overlap.csv`](genz2025_overlap.csv) traces each of its entries through our filters:
- 117 are outside our search (most at 2.5–3.0 Å resolution);
- 77 fail our quality filters;
- 19 involve proteins or pairs seen before the cutoff;
- 8 have no paper.

26 of the larger natural complexes in X-ray-hetero come from this overlap.

## Files

| Path | Contents |
| --- | --- |
| `CURATION.md` | criteria, funnel, size match, decision log, open points |
| `xray_hetero.csv`, `xray_homo.csv` | one row per dimer: chains, entities, UniProt, residues, resolution, release date, novelty reason, 30%/100% clusters, interface area, citation; `source` gives the size-matched or added group |
| `kept_xray_hetero.txt`, `kept_xray_homo.txt` | PDB IDs per set |
| `fasta_hetero/` | trimmed FASTA files of the 53 heterodimers (coordinates: biological assembly 1 of each PDB entry, chains as in the CSV) |
| `interface_ss_hetero.csv` | interface residues and their helix/strand fraction per chain (8PFC chain A: loops only) |
| `genz2025_overlap.csv` | where each Genz et al. 2025 heterodimer leaves our selection |
| `funnel.csv` | counts after each selection step |
| `results/scores_xray_hetero_all_weights.csv` | AF2Rank scores at full resolution, 0 recycles; one row per heterodimer and weight (53 × 17) |
| `results/scores_xray_hetero_masked_all_weights.csv` | AF2Rank scores in the paper's masked setting (sequence and side chains masked, 1 recycle); same layout |
| `scripts/` | `step1_fetch.py` and `step2_select.py` (RCSB search and selection), `rcsb_util.py`, `prep/` (input building, strict-cleaning check, interface secondary structure), `plot_scores.py` and `plot_masked.py` (figures) |

The selection uses the public RCSB APIs (metadata only) and reflects the PDB on 2026-11-21. The PDB archive data are [CC0](https://www.rcsb.org/pages/usage-policy).
