"""Render the E2A-HPr ensembles (examples/e2a_hpr) with PyMOL, headless: E2A once, every HPr model overlaid.

Usage: python render_ensembles.py <out_dir> examples/e2a_hpr/model_1ggr_nmr_ensemble.pdb examples/e2a_hpr/model_1hdn_ensemble_30.pdb
"""
import sys
from pymol import cmd

out, ggr, hdn = sys.argv[1], sys.argv[2], sys.argv[3]


def render(pdb, png, hpr_ribbon):
    cmd.reinitialize()
    cmd.bg_color("white"); cmd.set("ray_opaque_background", 1); cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2); cmd.set("ambient", 0.45); cmd.set("spec_reflect", 0.2)
    cmd.set_color("teal_c", [0x0f / 255, 0x9d / 255, 0x8a / 255]); cmd.set_color("amber_c", [0xe5 / 255, 0xa8 / 255, 0x3b / 255])
    cmd.load(pdb, "ens")
    cmd.remove("hydro or hetatm or solvent")
    cmd.dss("ens", state=0)
    # E2A from state 1 only; HPr from every state.
    cmd.create("e2a", "ens and chain A", 1, 1)
    cmd.split_states("ens", prefix="hpr_")
    cmd.delete("ens")
    hprs = [n for n in cmd.get_names("objects") if n.startswith("hpr_")]
    for n in hprs:
        cmd.remove(f"{n} and not chain B")
    cmd.hide("everything")
    cmd.show("cartoon", "e2a"); cmd.color("teal_c", "e2a")
    for n in hprs:
        if hpr_ribbon:
            cmd.show("ribbon", n); cmd.set("ribbon_width", 2.5, n)
        else:
            cmd.show("cartoon", n)
        cmd.color("amber_c", n)
    cmd.orient("e2a or hpr_*")
    import numpy as np
    ca = np.array(cmd.get_coords("e2a and name CA", 1)).mean(0)
    cb = np.array(cmd.get_coords(f"{hprs[0]} and name CA", 1)).mean(0)
    R = np.array(cmd.get_view()[:9]).reshape(3, 3)
    v = (cb - ca) @ R
    cmd.turn("z", -np.degrees(np.arctan2(v[1], v[0])))
    cmd.zoom("e2a or hpr_*", 2)
    cmd.png(png, width=1400, height=1100, dpi=300, ray=1)
    print("saved", png, len(hprs), "HPr models")




def render_all_states(pdb, png):
    """Every model of both chains as a thin ribbon, one shade per model, so the ensemble shows."""
    cmd.reinitialize()
    cmd.bg_color("white"); cmd.set("ray_opaque_background", 1); cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2); cmd.set("ambient", 0.5)
    cmd.load(pdb, "ens")
    cmd.remove("hydro or hetatm or solvent")
    cmd.split_states("ens", prefix="m_"); cmd.delete("ens")
    models = sorted(n for n in cmd.get_names("objects") if n.startswith("m_"))
    teal, amber = (0x0f, 0x9d, 0x8a), (0xe5, 0xa8, 0x3b)
    cmd.hide("everything")
    for i, n in enumerate(models):
        f = 0.55 + 0.45 * i / max(len(models) - 1, 1)  # darker to lighter per model
        cmd.set_color(f"t{i}", [c / 255 * f + (1 - f) * 0.15 for c in teal])
        cmd.set_color(f"a{i}", [c / 255 * f + (1 - f) * 0.15 for c in amber])
        cmd.show("ribbon", n); cmd.set("ribbon_width", 3.0, n)
        cmd.color(f"t{i}", f"{n} and chain A"); cmd.color(f"a{i}", f"{n} and chain B")
    cmd.orient(models[0])
    import numpy as np
    ca = np.array(cmd.get_coords(f"{models[0]} and chain A and name CA", 1)).mean(0)
    cb = np.array(cmd.get_coords(f"{models[0]} and chain B and name CA", 1)).mean(0)
    R = np.array(cmd.get_view()[:9]).reshape(3, 3)
    v = (cb - ca) @ R
    cmd.turn("z", -np.degrees(np.arctan2(v[1], v[0])))
    cmd.zoom("m_*", 2)
    cmd.png(png, width=1400, height=1100, dpi=300, ray=1)
    print("saved", png, len(models), "models")


render_all_states(ggr, f"{out}/e2a_hpr_1ggr_ensemble.png")
render(hdn, f"{out}/e2a_hpr_1hdn_ensemble.png", hpr_ribbon=True)
