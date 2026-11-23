"""X-ray-hetero natives at full resolution vs the AF2Rank paper's masked setting.

Full resolution: sequence and side chains kept, 0 recycles (results/scores_xray_hetero_all_weights.csv).
Masked: sequence and side chains masked, 1 recycle, as in Roney & Ovchinnikov 2022 and the authors'
AF2Rank notebook (results/scores_xray_hetero_masked_all_weights.csv).
Usage (from examples/xray_benchmark): python scripts/plot_masked.py
Writes figures/natives_masked_vs_fullres.png and figures/masked_dockq_io_heatmap.png and prints a summary.
"""
from pathlib import Path

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "figures"
frames = []
for setting, f in [("fullres_rec0", "scores_xray_hetero_all_weights.csv"), ("masked_rec1", "scores_xray_hetero_masked_all_weights.csv")]:
    frames.append(pd.read_csv(HERE / "results" / f).assign(setting=setting))
d = pd.concat(frames, ignore_index=True)
d["family"] = d.weight.str.replace(r"model_\d_", "", regex=True)

INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SET = {"fullres_rec0": ("#2a78d6", "full resolution, 0 recycles (our runs)"),
       "masked_rec1": ("#eb6834", "masked sequence + side chains, 1 recycle (AF2Rank paper)")}
FAM = [("ptm", "AF2-ptm\n(2 models)"), ("multimer", "AF2.1\nmultimer"), ("multimer_v2", "AF2.2\nmultimer"), ("multimer_v3", "AF2.3\nmultimer")]
plt.rcParams.update({"font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": MUTED, "axes.titlecolor": INK, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})

# ---- figure 1: distributions by family ----
METRICS = [("i_ptm", "ipTM", (0, 1), [(0.6, "0.6")]),
           ("actifptm", "actifpTM", (0, 1), [(0.6, "0.6")]),
           ("dockq_io", "DockQ input → output", (0, 1), [(0.23, "acceptable"), (0.8, "high")]),
           ("pae_contacts_mean", "PAE at native contacts (Å)", (0, 25), [])]
fig, axes = plt.subplots(1, 4, figsize=(16, 4.6))
rng = np.random.default_rng(0)
for ax, (col, lab, ylim, refs) in zip(axes, METRICS):
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    for y, t in refs:
        ax.axhline(y, color=MUTED, lw=1, ls=(0, (3, 3)), zorder=1)
        ax.text(4.42, y, t, color=MUTED, fontsize=7.5, va="bottom", ha="right")
    for i, (fam, _) in enumerate(FAM):
        for k, (s, (c, _)) in enumerate(SET.items()):
            v = d[(d.family == fam) & (d.setting == s)][col].dropna().values
            x = i + 1 + (k - 0.5) * 0.38
            ax.scatter(x + rng.uniform(-0.12, 0.12, len(v)), v, s=6, color=c, alpha=0.35, lw=0, zorder=2)
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.plot([x - 0.15, x + 0.15], [med, med], color=INK, lw=2, zorder=4, solid_capstyle="round")
            ax.plot([x, x], [q1, q3], color=INK, lw=1, zorder=4)
    ax.set_xticks(range(1, 5), [f for _, f in FAM]); ax.set_xlim(0.45, 4.55); ax.set_ylim(*ylim)
    ax.set_title(lab, fontsize=10)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
axes[0].set_ylabel("one dot per target × weight")
fig.legend(handles=[Patch(color=c, label=l) for c, l in SET.values()] +
           [plt.Line2D([], [], color=INK, lw=2, label="median (bar: interquartile range)")],
           loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.99))
fig.suptitle("Native X-ray heterodimers (53 targets × 17 weights): natives stay confident and in place when masked, except with AF2.1 multimer",
             y=1.06, fontsize=11.5, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(OUT / "natives_masked_vs_fullres.png", dpi=170, bbox_inches="tight"); plt.close(fig)

# ---- figure 2: masked DockQ_io per target x weight ----
m = d[d.setting == "masked_rec1"]
ORDER = ["model_1_ptm", "model_2_ptm"] + [f"model_{i}_{f}" for f in ["multimer", "multimer_v2", "multimer_v3"] for i in range(1, 6)]
WL = {w: ("ptm " if "ptm" in w else {"multimer": "2.1 ", "multimer_v2": "2.2 ", "multimer_v3": "2.3 "}[w.split("_", 2)[2]]) + w.split("_")[1]
      for w in ORDER}
piv = m.pivot_table(index="target", columns="weight", values="dockq_io").reindex(columns=ORDER)
piv = piv.loc[piv.median(axis=1).sort_values(ascending=False).index]
cmap = LinearSegmentedColormap.from_list("blues", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
fig, ax = plt.subplots(figsize=(8.2, 13))
im = ax.imshow(piv.values, cmap=cmap, vmin=0, vmax=1, aspect="auto")
bad = np.argwhere(piv.values < 0.23)
ax.scatter(bad[:, 1], bad[:, 0], marker="x", s=28, color="#d03b3b", lw=1.6, zorder=3)
for x in [1.5, 6.5, 11.5]:
    ax.axvline(x, color="#fcfcfb", lw=3)
ax.set_xticks(range(len(ORDER)), [WL[w] for w in ORDER], rotation=90, fontsize=8)
ax.set_yticks(range(len(piv)), piv.index, fontsize=7.5)
for s in ax.spines.values():
    s.set_visible(False)
cb = fig.colorbar(im, ax=ax, shrink=0.35, pad=0.02); cb.set_label("DockQ input → output (masked)"); cb.outline.set_visible(False)
ax.set_title("Masked setting: how much of each native interface AF2 keeps\n"
             "(targets sorted by median; × = DockQ_io < 0.23, native no longer acceptable)", fontsize=10, color=INK)
fig.tight_layout(); fig.savefig(OUT / "masked_dockq_io_heatmap.png", dpi=170, bbox_inches="tight"); plt.close(fig)
m = d[d.setting == "masked_rec1"]
print(d.groupby(["family", "setting"])[["i_ptm", "dockq_io"]].median().round(2).to_string())
print(m.groupby("family").agg(runs=("dockq_io", "size"), iptm_ge_0_6=("i_ptm", lambda x: (x >= 0.6).mean()),
                              lost=("dockq_io", lambda x: int((x < 0.23).sum()))).round(2).to_string())
print("wrote", OUT)
