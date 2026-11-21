"""Render E2A-HPr ensembles (examples/e2a_hpr) with PyMOL, headless: every dimer model as a cartoon, overlaid.

Usage: python render_ensembles.py <out_dir> examples/e2a_hpr/model_1ggr_nmr_ensemble.pdb examples/e2a_hpr/model_1hdn_ensemble_30.pdb
"""
import sys
import numpy as np
from pymol import cmd

out, ggr, hdn = sys.argv[1], sys.argv[2], sys.argv[3]
TEAL, AMBER = (0x0f, 0x9d, 0x8a), (0xe5, 0xa8, 0x3b)


def render(pdb, png):
    cmd.reinitialize()
    cmd.bg_color("white"); cmd.set("ray_opaque_background", 1); cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2); cmd.set("ambient", 0.45); cmd.set("spec_reflect", 0.2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.load(pdb, "ens")
    cmd.remove("hydro or hetatm or solvent")
    cmd.split_states("ens", prefix="m_"); cmd.delete("ens")
    models = sorted(n for n in cmd.get_names("objects") if n.startswith("m_"))
    cmd.hide("everything")
    for i, n in enumerate(models):
        cmd.dss(n)
        f = 0.75 + 0.25 * i / max(len(models) - 1, 1)  # slightly darker to lighter per model
        cmd.set_color(f"t{i}", [c / 255 * f for c in TEAL]); cmd.set_color(f"a{i}", [min(c / 255 * f + 0.05, 1) for c in AMBER])
        cmd.show("cartoon", n)
        cmd.color(f"t{i}", f"{n} and chain A"); cmd.color(f"a{i}", f"{n} and chain B")
    cmd.orient(models[0])
    ca = np.array(cmd.get_coords(f"{models[0]} and chain A and name CA", 1)).mean(0)
    cb = np.array(cmd.get_coords(f"{models[0]} and chain B and name CA", 1)).mean(0)
    R = np.array(cmd.get_view()[:9]).reshape(3, 3)
    v = (cb - ca) @ R
    cmd.turn("z", -np.degrees(np.arctan2(v[1], v[0])))
    cmd.zoom("m_*", 2)
    cmd.png(png, width=1400, height=1100, dpi=300, ray=1)
    print("saved", png, len(models), "dimer models")


render(ggr, f"{out}/e2a_hpr_1ggr_ensemble.png")
render(hdn, f"{out}/e2a_hpr_1hdn_ensemble.png")
