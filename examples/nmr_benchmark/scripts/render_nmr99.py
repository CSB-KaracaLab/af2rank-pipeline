"""Render model 1 of every NMR benchmark dimer as a cartoon thumbnail (PyMOL, headless).

Usage: python render_nmr99.py <nmr99_inputs dir> <out dir>
First chain teal, second amber, as in the report figures; the inter-chain axis is laid along screen x.
"""
import sys
from pathlib import Path
import numpy as np
from pymol import cmd

inp, out = Path(sys.argv[1]), Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
for pdb in sorted(inp.glob("*/models/*_model1.pdb")):
    t = pdb.parent.parent.name
    cmd.reinitialize()
    cmd.bg_color("white"); cmd.set("ray_opaque_background", 1); cmd.set("ray_shadows", 0); cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1); cmd.set("spec_reflect", 0.2); cmd.set("ambient", 0.45)
    cmd.set_color("teal_c", [0x0f/255, 0x9d/255, 0x8a/255]); cmd.set_color("amber_c", [0xe5/255, 0xa8/255, 0x3b/255])
    cmd.load(str(pdb), t); cmd.remove("hydro or hetatm or solvent"); cmd.dss(t)
    cmd.hide("everything"); cmd.show("cartoon")
    cmd.set("cartoon_tube_radius", 0.6); cmd.set("cartoon_loop_radius", 0.35)
    a, b = cmd.get_chains(t)[:2]  # file order: first chain teal, second amber
    cmd.color("teal_c", f"chain {a}"); cmd.color("amber_c", f"chain {b}")
    cmd.orient(t)
    ca = np.array(cmd.get_coords(f"chain {a} and name CA", 1)).mean(0)
    cb = np.array(cmd.get_coords(f"chain {b} and name CA", 1)).mean(0)
    R = np.array(cmd.get_view()[:9]).reshape(3, 3)
    v = (cb - ca) @ R
    cmd.turn("z", -np.degrees(np.arctan2(v[1], v[0])))
    cmd.zoom(t, 1.5, complete=1)
    cmd.png(str(out / f"{t}.png"), width=600, height=600, dpi=300, ray=1)
print("rendered", len(list(out.glob("*.png"))))
