"""Score vs DockQ input->output for the X-ray heterodimers, full resolution (top row) vs the paper's masked setting
(bottom row), one column per weight set.

The inputs are the crystal structures, so DockQ input->output is the true DockQ of AF2Rank's output.
Usage (from benchmarks/xray): python scripts/plot_masked_scatter.py [score: actifptm | i_ptm]
Small points: each weight; large points: dimer means over the weights of a set.
Writes figures/<score>_vs_dockq_full_vs_masked.png
"""
import sys
from pathlib import Path

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
SCORE = sys.argv[1] if len(sys.argv) > 1 else "actifptm"
LABEL = {"actifptm": "actifpTM", "i_ptm": "ipTM"}[SCORE]
d = pd.concat([pd.read_csv(HERE / "results/scores_xray_hetero_all_weights.csv").assign(setting="full"),
               pd.read_csv(HERE / "results/scores_xray_hetero_masked_all_weights.csv").assign(setting="masked")])
SETS = [("ptm", "AF2-ptm", "#52514e", "o"), ("multimer", "AF2.1 multimer", "#2a78d6", "s"),
        ("multimer_v2", "AF2.2 multimer", "#eb6834", "^"), ("multimer_v3", "AF2.3 multimer", "#1baf7a", "D")]
ROWS = [("full", "full resolution, 0 recycles"), ("masked", "masked seq. + side chains, 1 recycle (paper)")]
d["set"] = d["weight"].str.split("_", n=2).str[2]
INK, MUTED = "#0b0b0b", "#898781"

fig, axes = plt.subplots(2, 4, figsize=(16, 8.8), dpi=200, sharex=True, sharey=True)
for i, (setting, rlab) in enumerate(ROWS):
    for j, (key, name, col, mk) in enumerate(SETS):
        ax = axes[i, j]
        s = d[(d.setting == setting) & (d.set == key)]
        m = s.groupby("target")[[SCORE, "dockq_io"]].mean()
        ax.axvline(0.23, color=MUTED, lw=0.9, ls=(0, (3, 3)), zorder=1)
        ax.scatter(s["dockq_io"], s[SCORE], s=6, c=col, marker=mk, alpha=0.2, edgecolors="none", zorder=2)
        ax.scatter(m["dockq_io"], m[SCORE], s=46, c=col, marker=mk, edgecolors="white", linewidths=0.8, zorder=4)
        for t, r in m.iterrows():  # label dimers that leave the main cluster
            if r["dockq_io"] < 0.6 or r[SCORE] < 0.6:
                ax.annotate(t.upper(), (r["dockq_io"], r[SCORE]), fontsize=6.5, xytext=(3, 3),
                            textcoords="offset points", zorder=5, color=INK)
        lost = int((s["dockq_io"] < 0.23).sum())
        ax.set_title(f"{name}\nmean {LABEL} {s[SCORE].mean():.2f} · DockQ {s['dockq_io'].mean():.2f} · "
                     f"lost {lost}/{len(s)}", loc="left", fontsize=8.5, color=INK)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
        ax.grid(color="#e4e3df", linewidth=0.6, zorder=0)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        if i == 1:
            ax.set_xlabel("DockQ input→output (= true DockQ of the output)")
        if j == 0:
            ax.set_ylabel(f"{rlab}\n{LABEL}")
fig.suptitle(f"X-ray heterodimers (53 natives × 17 weights): {LABEL} vs DockQ, full resolution vs masked · "
             "small points: each weight, large: dimer means · dashed line: DockQ 0.23 (acceptable); "
             "labelled: dimer means below 0.6", x=0.01, ha="left", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.96))
out = HERE / f"figures/{SCORE}_vs_dockq_full_vs_masked.png"
fig.savefig(out); plt.close(fig); print("wrote", out)
