# Unseen X-ray dimer benchmark: curation record

Curated 2026-11-21. Two separate sets: **X-ray-homo** (`data/xray_homo.csv`, 39 homodimers) and **X-ray-hetero** (`data/xray_hetero.csv`, 53 heterodimers: 27 size-matched to the NMR heterodimers plus 26 larger natural complexes); IDs in `data/kept_xray_homo.txt` and `data/kept_xray_hetero.txt`. Selection scripts: `scripts/step1_fetch.py` and `scripts/step2_select.py` (RCSB metadata only). Intermediate tables named below (`xray_dimers*.csv`, the NMR curation table) are earlier working files and are not included.

## The logic

The NMR benchmark tests AF2Rank on interfaces AlphaFold2 never saw, but NMR models are flexible and sometimes not fully experimental. The X-ray set is a control: crystal structures held to the same novelty rule, the same quality rules we applied to the NMR set, and matched to it in size.

**Novel to every AlphaFold2 version.** Released after 2021-09-30, the PDB cutoff of AF2-multimer v3 (earlier versions have earlier cutoffs). Novelty is judged per protein (UniProt accession, plus its 100% sequence cluster so identical sequences under other accessions count):
- **Homodimers:** the protein had no structure released before the cutoff.
- **Heterodimers:** at least one partner had no earlier structure, or both had earlier structures but never together in one entry.

The rule follows the thesis' NMR curation but checks earlier structures of **any** method (X-ray, NMR, cryo-EM), not X-ray only.

**Same quality rules as the curated NMR set:**
- published (RCSB primary citation is a paper, not "To be published");
- no engineered mutations, no modified residues (selenomethionine allowed), no glycans, covalent cofactors or other modifications;
- no covalent link between the two chains;
- no peptides: every chain has at least 30 modelled residues;
- natural proteins: no synthetic constructs or nanobodies in either set.

**Clean, two-chain biology.** X-ray at 2.5 Å or better. Biological assembly 1 is exactly two protein chains (one entity ×2 or two entities ×1), no nucleic acid, and RCSB reports an interface between the chains. Same-protein fragment pairs (e.g. a protease and its own propeptide) are dropped; two constructs of one protein count as a homodimer.

**Non-redundant.** The best-resolution entry is kept per 30%-identity cluster (homodimers) or per pair of clusters (heterodimers). Heterodimers are further limited to one entry per cluster of the longer chain, so one domain family bound to different partners appears only once.

**Matched to the curated NMR benchmark.** The target is the 67 kept NMR dimers (39 homo, 28 hetero) as they stood before 6ES7 and 6DMP were excluded from the NMR set (now 65); the X-ray sets were not re-matched. Homodimers and heterodimers are selected separately by optimal one-to-one matching on log total length (fixed seed 20261121), and are kept as separate sets.

## Funnel

| Step | Entries |
| --- | --- |
| Search: X-ray, ≤ 2.5 Å, released after 2021-09-30, assembly 1 = 2 protein chains, no nucleic acid | 13,278 |
| No engineered mutation | 10,108 |
| No modified residue or modification | 8,184 |
| No covalent link between the chains | 7,963 |
| Every chain ≥ 30 residues, 69–900 in total | 6,805 |
| Interface reported | 6,792 |
| Two distinct proteins for heterodimers | 6,537 |
| Novel by X-ray history | 3,087 |
| Novel across all methods | 3,050 |
| Published | 2,166 |
| Non-redundant (907 homo + 263 hetero) | 1,170 |
| One heterodimer per long-chain cluster | 214 hetero |
| Natural proteins (UniProt on both chains) | 817 homo, 121 hetero |
| Selected | 39 homo + 28 hetero |

## Match to the curated NMR set (total residues)

| | X-ray median (range) | NMR median (range) | KS distance |
| --- | --- | --- | --- |
| Homodimers | 158 (83–865) | 158 (72–864) | 0.08 |
| Heterodimers, size-matched 28 | 174.5 (113–263) | 150 (93–263) | 0.29 |

The size-matched heterodimers run larger because few novel, published, natural, non-peptide X-ray heterodimers are small (the smallest in the pool has 113 residues). The 26 added heterodimers (161–798 residues) are deliberately not size-matched.

## Decision log (2026-11-21)

1. **v1:** first version. Every chain ≥ 30 residues, novelty judged by X-ray history only, matched to the original 99 NMR dimers.
2. **v2:** peptide partners admitted (shorter chain ≥ 8 residues) to mirror the NMR set's protein–peptide fraction. Novelty extended to all methods. Synthetic constructs and nanobodies flagged.
3. **v3:** heterodimers matched on total length plus shorter-chain length, one entry per long-chain cluster.
4. **v4:** the NMR set was curated (unpublished entries, docking models, a lipid-modified peptide, peptide complexes and a pseudocontact-shift-only complex removed). The X-ray set followed: published only, no peptides (every chain ≥ 30), matched to the 67 kept NMR dimers.
5. **Heterodimers without synthetic or nanobody chains:** the 10 heterodimers with a synthetic or nanobody chain (33ZJ, 6WMK, 7XLD, 8G4K, 9BK6, 9CC5, 9NKR, 9NLT, 9NZH, 9S83) were replaced by matched natural heterodimers (10VC, 7R49, 7VNI, 8PFC, 8SVG, 9E9G, 9L93, 9MLU, 9T3F, 9WLP). The previous table is `xray_dimers_v4_with_synthetic_hetero.csv`.

6. **Homodimers without synthetic chains:** the 9 homodimers with a synthetic-construct chain (8D03, 8D04, 8JVO, 8JVP, 8UPB, 8UTK, 9QDO, 9TLM, 9U9G) were each replaced by the unused natural homodimer from the pool closest in total length (7F9I, 8ZG9, 7YFV, 9IJX, 8CHW, 7WJT, 7RRM, 9N2G, 7Y86); the other 30 homodimers are unchanged. Previous table: `xray_dimers_v4b_with_synthetic_homo.csv`.
7. **Two sets:** homodimers and heterodimers are kept as separate benchmarks (X-ray-homo, X-ray-hetero), as for the NMR set.
8. **Larger heterodimers added:** X-ray can determine larger complexes than NMR, so size matching is not required for heterodimers. 38 heterodimers of the benchmark of Genz et al. 2025 (Protein Sci 34:e70327, doi:10.1002/pro.70327) passed all our filters but were not picked by size (161–798 residues). The 11 with a designed binder or a chain without UniProt were excluded and the 27 natural ones added; 8CM0 was then dropped because its chain B coordinates match the SEQRES sequence at only 89% identity, failing the pipeline's strict cleaning. X-ray-hetero is now 54 dimers (`source` column: size-matched or Genz et al. overlap). Three added entries share one partner family with a size-matched entry but with a different partner (8JTK/8PFC, 6ZH1/9MLU, 7V7L/7VNI), which the pair-level and longer-chain redundancy rules allow. The full comparison with their 262 heterodimers is in `data/genz2025_overlap.csv`.
9. **8SJJ removed:** the SEPT14–SEPT7 coiled coil is described by its authors as a metastable heterodimer, so its crystal interface is not a stable native complex. X-ray-hetero is now 53 dimers (27 size-matched, 26 larger).

## Open points

- **9ZDI** (homodimer, 52+52) is nearly all coil in the deposited model; worth checking before it is scored.
- **Symmetry copies:** 13 dimers are built from one deposited chain plus a crystal-symmetry copy and need the assembly file.
- **Entry to look at:** 7WZ6 (long crossed helices) may be a construct or packing arrangement.
- **8PFC:** chain A touches its partner only through loops (no helix or strand at the interface); kept.
- **Runs:** X-ray-hetero is scored in both settings (see README). X-ray-homo inputs are not built yet.
