"""Buried surface area (BSA) of each X-ray heterodimer, from the exact two-chain inputs AF2Rank scored.

BSA = SASA(chain A) + SASA(chain B) - SASA(complex), summed over both sides (Å²), heavy atoms only,
Shrake-Rupley with Biopython defaults (probe 1.40 Å, 100 points per atom).
Usage: python scripts/prep/buried_surface.py <inputs dir with <id>/models/<id>_model1.pdb> [out csv]
Default output: interface_bsa_hetero.csv (only the ids in kept_xray_hetero.txt).
"""
import sys
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

HERE = Path(__file__).resolve().parents[2]
inputs = Path(sys.argv[1])
out = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "interface_bsa_hetero.csv"
kept = [x.strip().lower() for x in (HERE / "kept_xray_hetero.txt").read_text().split() if x.strip()]
sr, parser = ShrakeRupley(), PDBParser(QUIET=True)


def sasa(entity):
    sr.compute(entity, level="A")
    return sum(a.sasa for a in entity.get_atoms())


rows = ["pdb_id,chain_A,chain_B,sasa_A,sasa_B,sasa_complex,bsa_total,bsa_per_side"]
for pid in kept:
    model = parser.get_structure(pid, inputs / pid / "models" / f"{pid}_model1.pdb")[0]
    chains = list(model)
    if len(chains) != 2:
        raise SystemExit(f"{pid}: expected 2 chains, found {len(chains)}")
    a, b = (sasa(c.copy()) for c in chains)
    ab = sasa(model)
    bsa = a + b - ab
    rows.append(f"{pid},{chains[0].id},{chains[1].id},{a:.1f},{b:.1f},{ab:.1f},{bsa:.1f},{bsa / 2:.1f}")
    print(pid, round(bsa), flush=True)
out.write_text("\n".join(rows) + "\n")
print("wrote", out)
