"""ipTM vs DockQ input->output, one row per curated NMR sub-benchmark (heterodimers, homodimers, flexible), one column per weight set.

Usage (from examples/nmr_benchmark):
  python scripts/plot_subsets_grid.py results/scores_model1_all_weights.csv targets.csv figures/iptm_vs_dockq_by_subbenchmark.png
Only dimers with status "kept" in targets.csv are plotted.
"""
import sys

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

d = pd.read_csv(sys.argv[1]); k = pd.read_csv(sys.argv[2]); out = sys.argv[3]
k = k[k["status"] == "kept"].copy(); k["target"] = k["pdb_id"].str.lower(); d["target"] = d["pdb_id"].str.lower()
k["sub"] = k.apply(lambda r: "flexible" if r["long_disordered_segments"] else ("homo" if r["dimer"] == "homodimer" else "hetero"), axis=1)
d = d.merge(k[["target", "sub"]], on="target")
ROWS = [("hetero", "NMR heterodimers"), ("homo", "NMR homodimers (symmetric)"), ("flexible", "NMR-flexible (long disordered segments)")]
SETS = [("AF2-ptm", "#52514e", "o"), ("AF2.1 multimer", "#2a78d6", "s"),
        ("AF2.2 multimer", "#eb6834", "^"), ("AF2.3 multimer", "#1baf7a", "D")]


def rho(x):
    return x[["dockq_io", "iptm"]].corr(method="spearman").iloc[0, 1] if len(x) > 2 else float("nan")


fig, axes = plt.subplots(3, 4, figsize=(15, 13.2), dpi=200, sharex=True, sharey=True)
for i, (sub, label) in enumerate(ROWS):
    n = k[k["sub"] == sub]
    for j, (setname, col, mk) in enumerate(SETS):
        ax = axes[i, j]; s = d[(d["sub"] == sub) & (d["set"] == setname)]
        title = f"{setname} · ρ {rho(s):.2f}\nmean ipTM {s['iptm'].mean():.2f} · DockQ {s['dockq_io'].mean():.2f}"
        if sub == "flexible":  # mixed row: heterodimers filled, homodimers hollow
            homo = s["target"].isin(set(k.loc[k["dimer"] == "homodimer", "target"]))
            ax.scatter(s.loc[~homo, "dockq_io"], s.loc[~homo, "iptm"], s=13, c=col, marker=mk, alpha=0.55, edgecolors="none", zorder=3)
            ax.scatter(s.loc[homo, "dockq_io"], s.loc[homo, "iptm"], s=16, facecolors="none", edgecolors=col, marker=mk,
                       linewidths=0.8, zorder=4)
            title += f"\nhetero ρ {rho(s[~homo]):.2f} · homo ρ {rho(s[homo]):.2f}"
        else:
            ax.scatter(s["dockq_io"], s["iptm"], s=13, c=col, marker=mk, alpha=0.55, edgecolors="none", zorder=3)
        ax.set_title(title, loc="left", fontsize=8.5)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
        ax.grid(color="#e4e3df", linewidth=0.6, zorder=0)
        for side in ("top", "right"): ax.spines[side].set_visible(False)
        if i == 2: ax.set_xlabel("DockQ input→output")
    detail = f"{len(n)} dimers" + (f" ({int((n['dimer'] == 'homodimer').sum())} homo, {int((n['dimer'] != 'homodimer').sum())} hetero)" if sub == "flexible" else "")
    axes[i, 0].set_ylabel(f"{label}\n{detail}\n\nipTM", fontsize=9.5)
fig.suptitle(f"Curated NMR benchmark ({len(k)} dimers) · AF2Rank, full resolution, 0 recycles · one point per weight and dimer",
             x=0.01, ha="left", fontsize=11)
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([], [], ls="", marker="o", color="#52514e", alpha=0.6, label="heterodimer"),
                    Line2D([], [], ls="", marker="o", mfc="none", mec="#52514e", label="homodimer")],
           title="NMR-flexible row", loc="lower center", ncol=2, frameon=False, fontsize=9, title_fontsize=9)
fig.tight_layout(rect=(0, 0.035, 1, 0.97), h_pad=2.2); fig.savefig(out); plt.close(fig)
print("wrote", out)
