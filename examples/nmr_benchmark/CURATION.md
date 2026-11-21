# NMR dimer benchmark: curation record

Curated 2026-11-21. Start: the 99 NMR dimers of B. Savaş's AF2Rank chapter. Kept: 65; excluded: 34.

## Primary curation criteria (B. Savaş, PhD thesis chapter 4.1.2)

The starting 99 dimers were selected in the chapter as follows (paraphrased):

1. NMR structures only, since no AF2 version was trained on NMR-derived structures; dimers only.
2. Screened 336 NMR dimers released between March 2016 and December 2025 (of 1,508 NMR dimers in the PDB).
3. Novel interface: at least one monomer had no X-ray structure before the AF2 training cutoff; if both monomers had earlier X-ray structures, they never appeared together in one X-ray complex before the cutoff (160 left).
4. No mutations, no post-translational modifications, no covalently bound interfaces (122 left).
5. Secondary structure at the interface required, to leave out intrinsically disordered interactions (99 left: 40 homodimers, 59 heterodimers).

Source: B. Savaş, PhD thesis (in preparation), chapter 4; per-target scores of the original 99 at github.com/busrasavas/AF2Rank-PPI.

## Kept, by sub-benchmark

**NMR-homo** (29): symmetric homodimers, no long disordered segments

2NDP 5JPW 5NF8 5NOC 5YAM 6C44 6D6S 6LUL 6LXG 6N2M 6Q2Z 6SAI 6TDM 6TV5 6U6P 6U6Q 6U6S 6UHW 6V88 7C4O 7DFE 7VU7 7WJ0 7YWQ 8DSB 8DSX 8PEK 8YAP 8ZU6

**NMR-hetero** (12): heterodimers, no long disordered segments

2N77 2NBW 5L85 5VF0 6B7G 6BUT 6DSL 6IVU 7CSQ 7QUU 7ZEY 8BA1

**NMR-flexible** (24): homo- and heterodimers with long disordered segments (curators' list)

2N73 2N9P 5N7Y 5TMX 5US5 5XV8 6E4H 6EVI 6FEG 6LQZ 6TDD 6TDN 6U19 6U4N 6XFL 6ZDB 7DTI 7JYN 7JYZ 7PKU 7QRI 8RAJ 8TT7 9R3Y

## Excluded

| PDB | Reason |
|---|---|
| 2N23 | no published paper (PDB primary citation: To be published) |
| 2N3K | peptide partner (chain shorter than 30 residues) |
| 2N4Q | no published paper (PDB primary citation: To be published) |
| 2N7K | no published paper (PDB primary citation: To be published) |
| 2NBV | assembled from earlier structures (1J8C, 5IRS) with HADDOCK |
| 2NCZ | peptide partner (chain shorter than 30 residues) |
| 2ND0 | peptide partner (chain shorter than 30 residues) |
| 5GWM | no published paper (PDB primary citation: To be published) |
| 5H7P | peptide partner (chain shorter than 30 residues) |
| 5VTO | no published paper (PDB primary citation: To be published) |
| 5XBO | complex built from lanthanide pseudocontact shifts only (no NOESY; methods paper) |
| 5ZAU | HADDOCK docking model (docking with NMR restraints), not an experimentally determined complex |
| 6BGG | peptide partner (chain shorter than 30 residues) |
| 6BNH | peptide partner (chain shorter than 30 residues) |
| 6CTB | no published paper (PDB primary citation: To be published) |
| 6DMP | de novo designed heterodimer (not a natural protein pair) |
| 6E5N | peptide partner (chain shorter than 30 residues) |
| 6ES7 | reconstructed ancestral complex of two intrinsically disordered proteins (low-affinity, dynamic) |
| 6JXU | HADDOCK complex from intermolecular NOEs only (SUMO1 coordinates likely from an earlier structure; no intramolecular restraints deposited) |
| 6JXW | HADDOCK complex from intermolecular NOEs only (SUMO2 coordinates likely from an earlier structure; no intramolecular restraints deposited) |
| 6M78 | assembled from earlier structures (1ZYO, 6LXF) with HADDOCK |
| 6S3W | peptide partner (chain shorter than 30 residues) |
| 7JQ8 | peptide partner (chain shorter than 30 residues) |
| 7K3S | no published paper (PDB primary citation: To be published) |
| 7NQC | lipid-modified peptide: farnesyl at the interface (PTM) |
| 7OVC | peptide partner (chain shorter than 30 residues) |
| 7QDW | peptide partner (chain shorter than 30 residues) |
| 7SFT | peptide partner (chain shorter than 30 residues) |
| 7X5C | peptide partner (chain shorter than 30 residues) |
| 8DGH | assembled from an earlier structure (calmodulin lobe from 2F3Y) with HADDOCK |
| 8DGK | assembled from an earlier structure (calmodulin lobe from 2F3Y) with HADDOCK |
| 8HEQ | peptide partner (chain shorter than 30 residues) |
| 8U2Y | peptide partner (chain shorter than 30 residues) |
| 9ATN | peptide partner (chain shorter than 30 residues) |

## Kept, with notes

- 5N7Y: interface ligand/ion: Zn ions coordinated mainly by one chain, 3 contacts to the partner
- 7QRI: interface ligand/ion: L-Phe bound between the two chains (2 copies)
- 7CSQ: low AF2.3 ipTM (< 0.6 with all five weights); deposited method is NOE-based, so kept.
- 7YWQ: HADDOCK listed only as a processing step; computed and refined with X-PLOR-NIH from 3D 13C/15N NOESY data.

## Decision log (2026-11-21)

1. Long disordered segments: 33 dimers named by the curators form the NMR-flexible sub-benchmark (7JNY in that list read as 7JYN).
2. 7NQC excluded: the RalA peptide carries a farnesyl group (ligand ULW on the C-terminal Cys) buried in the calmodulin pocket; cleaning removed it, and a lipid modification breaks the set's no-PTM rule.
3. 5ZAU excluded: RCSB refinement 'Docking with NMR restraints' (HADDOCK); FYN SH3 taken from an existing structure, interface mainly from titration chemical-shift perturbations.
4. Entries without a paper excluded (RCSB primary citation 'To be published'): 2N23, 2N4Q, 2N7K, 5GWM, 5VTO, 6CTB, 7K3S.
5. Assembled from earlier structures with HADDOCK, excluded: 8DGH and 8DGK (calmodulin lobe from 2F3Y), 2NBV (1J8C, 5IRS), 6M78 (1ZYO, 6LXF); 6CTB (rigid docking) is already out under rule 4.
6. Interface ligands checked in model 1 of every entry: 7NQC (excluded, see 2), 7QRI (L-Phe between the chains) and 5N7Y (zinc, mostly one chain) kept with a note.
7. 6JXU and 6JXW excluded (SUMO1/SUMO2 with an ICP0 SIM peptide): computed with HADDOCK; the deposited restraints are only 17 and 19 intermolecular NOEs (no intramolecular restraints), so the SUMO coordinates likely come from an earlier structure.
8. Peptide complexes excluded (a chain shorter than 30 modelled residues), as in the X-ray set: 2N3K, 2NCZ, 2ND0, 5H7P, 6BGG, 6BNH, 6E5N, 6S3W, 7JQ8, 7OVC, 7QDW, 7SFT, 7X5C, 8HEQ, 8U2Y, 9ATN.
9. Deposited methods checked for the low-ipTM heterodimers: 6ES7 (CYANA, 3D NOESY + RDCs) and 7CSQ (CYANA/Amber, 13C/15N-filtered and 4D NOESY) kept; 5XBO excluded (only HSQC and 49 lanthanide pseudocontact-shift restraints, no NOESY; a lanthanide-tagging methods paper).
10. Natural proteins only, as in the X-ray set: 6ES7 excluded (reconstructed ancestral CID/NCBD complex of two intrinsically disordered proteins; low-affinity and dynamic) and 6DMP excluded (de novo designed heterodimer).
