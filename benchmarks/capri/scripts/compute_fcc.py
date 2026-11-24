"""Decoy-vs-decoy interface similarity: fraction of common contacts (FCC), as used by HADDOCK's clustering
(Rodrigues et al., Proteins 80:1810, 2012).

Contacts are residue pairs across the two chains with any heavy atoms within 5 Å. FCC(i, j) = |Ci ∩ Cj| / |Ci|.
Input: one multi-model PDB per target (the Scorer decoys trimmed to the residues resolved in the native, MODEL
serials = scoreset model numbers), at <decoy dir>/<target>/models/*.pdb. The decoys are not included in this
repository (scoreset.org); the resulting matrices are, in results/fcc/fcc_<target>.npz (ids, fcc, ncontacts).
Usage (from benchmarks/capri): python scripts/compute_fcc.py <decoy dir>
"""
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parents[1]


def models(path):
    cur, serial = [], None
    for line in open(path):
        if line.startswith("MODEL"):
            serial, cur = int(line.split()[1]), []
        elif line.startswith(("ATOM", "HETATM")) and line[76:78].strip() != "H" and line[12:16].strip()[0] != "H":
            cur.append((line[21], line[22:27], float(line[30:38]), float(line[38:46]), float(line[46:54])))
        elif line.startswith("ENDMDL"):
            yield serial, cur


def contacts(atoms):
    chains = sorted({a[0] for a in atoms})
    a_atoms = [a for a in atoms if a[0] == chains[0]]; b_atoms = [a for a in atoms if a[0] == chains[1]]
    pairs = cKDTree(np.array([a[2:] for a in a_atoms])).query_ball_tree(cKDTree(np.array([b[2:] for b in b_atoms])), 5.0)
    return {(a_atoms[i][1], b_atoms[j][1]) for i, js in enumerate(pairs) for j in js}


if __name__ == "__main__":
    out = HERE / "results/fcc"; out.mkdir(parents=True, exist_ok=True)
    for d in sorted(Path(sys.argv[1]).iterdir()):
        ids, contact_sets = [], []
        for serial, atoms in models(next((d / "models").glob("*.pdb"))):
            ids.append(f"{d.name}_scorers_model_{serial:04d}"); contact_sets.append(contacts(atoms))
        n = len(ids); fcc = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                fcc[i, j] = len(contact_sets[i] & contact_sets[j]) / len(contact_sets[i]) if contact_sets[i] else 0.0
        np.savez_compressed(out / f"fcc_{d.name}.npz", ids=np.array(ids), fcc=fcc,
                            ncontacts=np.array([len(c) for c in contact_sets]))
        print(d.name, n, "decoys")
