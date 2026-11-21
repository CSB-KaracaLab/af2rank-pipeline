"""Two-panel report figure: input structure (A) beside ipTM vs DockQ input->output for 17 AF2 weights (B).

Usage: python make_figure.py --structure 2ptc_structure.png --scores examples/2ptc/results/scores_all_weights.csv \
         --out figure_2ptc.png --title "..." --chain-a "Trypsin, chain E (223 aa)" --chain-b "BPTI, chain I (58 aa)"
"""
import argparse
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from PIL import Image, ImageChops

ap = argparse.ArgumentParser()
ap.add_argument("--structure"); ap.add_argument("--scores"); ap.add_argument("--out")
ap.add_argument("--title"); ap.add_argument("--chain-a"); ap.add_argument("--chain-b"); ap.add_argument("--note", default="")
args = ap.parse_args()

SETS = [("AF2-ptm", "ptm", "#52514e", "o"), ("AF2.1", "v1", "#2a78d6", "s"),
        ("AF2.2", "v2", "#eb6834", "^"), ("AF2.3", "v3", "#1baf7a", "D")]
t = pd.read_csv(args.scores)
t["set"] = t["Weight"].str.split().str[0]; t["n"] = t["Weight"].str.split().str[-1]
x, y = "DockQ input→output", "ipTM"
models = sorted(t["Model"].unique())
SHAPES = ["o", "s", "^"]

img = Image.open(args.structure).convert("RGB")
bbox = ImageChops.difference(img, Image.new("RGB", img.size, (255, 255, 255))).getbbox()
img = img.crop((max(bbox[0] - 30, 0), max(bbox[1] - 30, 0), min(bbox[2] + 30, img.width), min(bbox[3] + 30, img.height)))

fig = plt.figure(figsize=(13, 5.6), dpi=200)
fig.patch.set_facecolor("white")
axA = fig.add_axes([0.01, 0.13, 0.44, 0.74]); axB = fig.add_axes([0.53, 0.11, 0.44, 0.76])
axA.imshow(img); axA.axis("off")
axA.legend(handles=[Patch(color="#0f9d8a", label=args.chain_a), Patch(color="#e5a83b", label=args.chain_b)],
           loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=2, frameon=False, fontsize=9)
for ax, letter in ((axA, "A"), (axB, "B")):
    fig.text(ax.get_position().x0, 0.905, letter, fontsize=15, fontweight="bold", va="top")

for name, tag, colour, marker in SETS:
    for i, model in enumerate(models):
        s = t[(t["set"] == name) & (t["Model"] == model)]
        axB.scatter(s[x], s[y], s=48, c=colour, marker=marker if len(models) == 1 else SHAPES[i],
                    edgecolors="white", linewidths=1.0, zorder=3)
axB.set_xlim(0, 1); axB.set_ylim(0, 1)
axB.set_xlabel(x); axB.set_ylabel(y)
axB.grid(color="#e4e3df", linewidth=0.6, zorder=0)
for side in ("top", "right"): axB.spines[side].set_visible(False)
for side in ("left", "bottom"): axB.spines[side].set_color("#c3c2b7")

# one label per weight (mean of the models), placed where it clears labels and markers
fig.canvas.draw(); renderer = fig.canvas.get_renderer(); box_ax = axB.get_window_extent(renderer)
taken = [matplotlib.transforms.Bbox.from_extents(px - 5, py - 5, px + 5, py + 5)
         for px, py in axB.transData.transform(t[[x, y]].to_numpy())]
groups = t.groupby("Weight").agg(x=(x, "mean"), y=(y, "mean"), set=("set", "first"), n=("n", "first"))
tags = {s[0]: s[1] for s in SETS}
for _, g in groups.sort_values("y", ascending=False).iterrows():
    placed = False
    for dist in (8, 14, 21, 30, 40, 52):
        for dx, dy in ((1, 0), (1, 1), (1, -1), (-1, 0), (-1, 1), (-1, -1), (0, 1), (0, -1)):
            note = axB.annotate(f"{tags[g['set']]}·{g['n']}", (g["x"], g["y"]), xytext=(dx * dist, dy * dist * 0.8),
                                textcoords="offset points", fontsize=7.5, color="#52514e", va="center",
                                ha="left" if dx > 0 else "right" if dx < 0 else "center", zorder=4,
                                arrowprops=dict(arrowstyle="-", color="#c3c2b7", lw=0.6, shrinkA=0, shrinkB=4) if dist > 8 else None)
            b = note.get_window_extent(renderer).expanded(1.03, 1.08)
            if box_ax.contains(b.x0, b.y0) and box_ax.contains(b.x1, b.y1) and not any(b.overlaps(o) for o in taken):
                taken.append(b); placed = True; break
            note.remove()
        if placed: break

handles = [plt.Line2D([], [], ls="", marker=m, ms=7, mfc=c, mec=c, label=n + ("" if n == "AF2-ptm" else " multimer"))
           for n, _, c, m in SETS]
if len(models) > 1:
    handles = [plt.Line2D([], [], ls="", marker="o", ms=7, mfc=c, mec=c, label=n + ("" if n == "AF2-ptm" else " multimer"))
               for n, _, c, m in SETS]
    leg = axB.legend(handles=handles, title="Weight set (colour)", loc="upper left", frameon=False, fontsize=8, title_fontsize=8)
    axB.add_artist(leg)
    axB.legend(handles=[plt.Line2D([], [], ls="", marker=SHAPES[i], ms=7, mfc="#8a8984", mec="#8a8984", label=f"NMR model {i + 1}")
                        for i in range(len(models))], title="Model (shape)", loc="lower right", frameon=False, fontsize=8, title_fontsize=8)
else:
    axB.legend(handles=handles, loc="upper left", frameon=False, fontsize=8)
fig.suptitle(args.title, x=0.01, y=0.985, ha="left", fontsize=12)
if args.note:
    fig.text(0.53, 0.01, args.note, fontsize=7, color="#8a8984")
fig.savefig(args.out, facecolor="white")
print("saved", args.out)
