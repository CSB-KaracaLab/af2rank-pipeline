"""Render the two input complexes for the report figures with PyMOL (headless).

Usage: python render_structures.py <out_dir> examples/2ptc/model.pdb examples/e2a_hpr/model_1ggr_nmr_ensemble.pdb
"""
import sys
from pymol import cmd
out, inputs = sys.argv[1], sys.argv[2:]
def setup():
    cmd.reinitialize()
    cmd.bg_color("white"); cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_transparency", 0); cmd.set("ray_shadows", 0); cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1); cmd.set("spec_reflect", 0.2); cmd.set("ambient", 0.45)
    cmd.set_color("teal_c", [0x0f/255, 0x9d/255, 0x8a/255]); cmd.set_color("amber_c", [0xe5/255, 0xa8/255, 0x3b/255])

def render(pdb, name, chains, png, states=False):
    setup()
    cmd.load(pdb, name)
    cmd.remove("hydro or hetatm or solvent")
    cmd.dss(name)
    cmd.hide("everything"); cmd.show("cartoon")
    if states:
        cmd.set("all_states", 1)
    a, b = chains
    cmd.color("teal_c", f"chain {a}"); cmd.color("amber_c", f"chain {b}")
    # interface side chains (within 4.5 A of the partner) as thin sticks
    cmd.select("iface", f"byres ((chain {a} within 4.5 of chain {b}) or (chain {b} within 4.5 of chain {a})) and not name N+C+O", state=1)
    cmd.show("sticks", "iface and not name N+C+O"); cmd.set("stick_radius", 0.18)
    cmd.util.cnc("iface")
    cmd.orient(name)
    # put chain A left, chain B right: rotate so the inter-chain vector lies along screen x
    import numpy as np
    ca = np.array(cmd.get_coords(f"chain {a} and name CA", 1)).mean(0)
    cb = np.array(cmd.get_coords(f"chain {b} and name CA", 1)).mean(0)
    view = cmd.get_view(); R = np.array(view[:9]).reshape(3, 3)
    v = (cb - ca) @ R  # model->camera
    ang = np.degrees(np.arctan2(v[1], v[0]))
    cmd.turn("z", -ang)
    cmd.zoom(name, 2)
    cmd.png(png, width=1400, height=1100, dpi=300, ray=1)
    print("saved", png)

render(inputs[0], "ptc", ("E", "I"), f"{out}/2ptc_structure.png")
render(inputs[1], "ggr", ("A", "B"), f"{out}/1ggr_structure.png", states=True)
