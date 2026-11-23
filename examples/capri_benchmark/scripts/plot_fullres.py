"""AF2Rank scores vs CAPRI's own DockQ and classes for the CAPRI Scorer decoys at full resolution (0 recycles).

Reads results/scores_capri_fullres_all_weights.csv. Writes figures/fullres_<score>_vs_dockq.png (one weight per
AF2 family x 9 targets), figures/fullres_ranking_summary.png (heatmaps of Spearman, AUC and top-10 hits for every
target x weight) and results/ranking_summary.csv.
Usage (from examples/capri_benchmark): python scripts/plot_fullres.py
"""
import sys
from pathlib import Path

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import HERE, load, targets  # noqa: E402

OUT = HERE / "figures"
d = load("fullres")
d["hit"] = d.good
sel = targets()
TARGETS = sorted(d.target.unique())
FAM = [("ptm", "AF2-ptm"), ("multimer", "AF2.1"), ("multimer_v2", "AF2.2"), ("multimer_v3", "AF2.3")]
WEIGHTS = ["model_1_ptm", "model_2_ptm"] + [f"model_{i}_{f}" for f, _ in FAM[1:] for i in range(1, 6)]
WLABEL = {w: f"{dict(FAM)[w.split('_', 2)[2]]} model {w.split('_')[1]}" for w in WEIGHTS}
SCATTER_WEIGHTS = ["model_1_ptm", "model_1_multimer", "model_1_multimer_v2", "model_1_multimer_v3"]  # one per family
# CAPRI quality is ordinal: incorrect recedes in grey, then a light -> dark blue ramp.
CLASSES = [("incorrect", "#c9c7c1", 10), ("acceptable", "#86b6ef", 16), ("medium", "#2a78d6", 16), ("high", "#0d366b", 16)]
INK, MUTED = "#2b2b29", "#76746e"

plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#b5b3ad", "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK})


def tlabel(t):
    tid = t.replace("_", ".")
    r = sel.loc[tid]
    return f"{tid} ({r.pdb})\n{r.desc.split(' ')[0][:22]}"


def scatter_grid(score, label):
    fig, axes = plt.subplots(len(SCATTER_WEIGHTS), len(TARGETS), figsize=(2.05 * len(TARGETS), 2.0 * len(SCATTER_WEIGHTS)),
                             sharex=True, sharey=True)
    for i, w in enumerate(SCATTER_WEIGHTS):
        for j, t in enumerate(TARGETS):
            ax = axes[i, j]
            g = d[(d.weight == w) & (d.target == t)]
            ax.grid(color="#ecebe7", lw=0.6); ax.set_axisbelow(True)
            if g.empty:
                ax.text(0.5, 0.5, "not run yet", ha="center", va="center", color=MUTED, transform=ax.transAxes)
            else:
                for cls, col, s in CLASSES:
                    h = g[g.main_classification == cls]
                    ax.scatter(h.dockq, h[score], s=s, c=col, lw=0.5 if cls != "incorrect" else 0,
                               edgecolors="#fcfcfb", zorder=3 if cls != "incorrect" else 2)
                rho = spearmanr(g.dockq, g[score])[0]
                ax.text(0.97, 0.05, f"ρ = {rho:.2f}", ha="right", va="bottom", transform=ax.transAxes, color=INK, fontsize=8)
            ax.set_xlim(-0.03, 1.0); ax.set_ylim(0, 1.02)
            if i == 0:
                ax.set_title(tlabel(t), fontsize=8.5)
            if j == 0:
                ax.set_ylabel(f"{WLABEL[w]}\n{label}")
            if i == len(SCATTER_WEIGHTS) - 1:
                ax.set_xlabel("CAPRI DockQ")
    handles = [plt.Line2D([], [], ls="", marker="o", ms=6, mfc=c, mec="none", label=n) for n, c, _ in CLASSES]
    fig.legend(handles=handles, title="CAPRI class", loc="upper center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, 1.0))
    fig.suptitle(f"AF2Rank {label} vs CAPRI DockQ at full resolution: CAPRI Scorer decoys, model 1 of each AF2 family",
                 y=1.025, fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(OUT / f"fullres_{score}_vs_dockq.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


# Ranking metrics per target x weight x score
rows = []
for (t, w), g in d.groupby(["target", "weight"]):
    n_hit = int(g.hit.sum())
    for score in ["i_ptm", "actifptm", "ptm", "plddt", "neg_pae_contacts"]:
        o = g.sort_values(score, ascending=False)
        rows.append({
            "target": t, "weight": w, "score": score, "n_decoys": len(g), "n_hits": n_hit,
            "spearman_dockq": spearmanr(g.dockq, g[score])[0],
            "auc_hit": roc_auc_score(g.hit, g[score]) if 0 < n_hit < len(g) else np.nan,
            "top1_class": o.main_classification.iloc[0],
            "top1_dockq": o.dockq.iloc[0],
            "hits_in_top10": int(o.hit.head(10).sum()),
            "random_top10": 10 * n_hit / len(g),
            "best_dockq_in_top10": o.dockq.head(10).max(),
            "best_dockq_overall": g.dockq.max(),
        })
summ = pd.DataFrame(rows)
summ.to_csv(HERE / "results/ranking_summary.csv", index=False)


def dark_cell(frac, cmap):
    # text on a cell: white only where the fill is dark enough
    return frac > 0.82 or frac < 0.15 if cmap == "RdBu_r" else frac > 0.55


def heat(ax, score, metric, title, vmin, vmax, cmap, fmt):
    m = summ[summ.score == score].pivot(index="weight", columns="target", values=metric).reindex(index=WEIGHTS, columns=TARGETS)
    im = ax.imshow(m.values, vmin=vmin, vmax=vmax, cmap=cmap, aspect="auto")
    for (i, j), v in np.ndenumerate(m.values):
        if np.isnan(v):
            ax.text(j, i, "–", ha="center", va="center", color=MUTED)
        else:
            frac = (v - vmin) / (vmax - vmin)
            ax.text(j, i, fmt(v), ha="center", va="center", fontsize=8,
                    color="white" if dark_cell(frac, cmap) else INK)
    ax.set_xticks(range(len(TARGETS)), [t.replace("_", ".") for t in TARGETS], rotation=0, fontsize=8)
    ax.set_yticks(range(len(WEIGHTS)), [WLABEL[w] for w in WEIGHTS])
    ax.set_title(title, fontsize=9.5)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


ROWS = [("i_ptm", "ipTM"), ("actifptm", "actifpTM"), ("neg_pae_contacts", "−PAE at input contacts")]
fig, axes = plt.subplots(len(ROWS), 3, figsize=(19, 6.2 * len(ROWS)))
for r, (score, lab) in enumerate(ROWS):
    heat(axes[r, 0], score, "spearman_dockq", f"{lab}: Spearman ρ with CAPRI DockQ", -1, 1, "RdBu_r", lambda v: f"{v:.2f}")
    heat(axes[r, 1], score, "auc_hit", f"{lab}: AUC, acceptable-or-better vs incorrect", 0, 1, "RdBu_r", lambda v: f"{v:.2f}")
    top = summ[summ.score == score].assign(lbl=lambda x: x.hits_in_top10.astype(str) + " / " + x.random_top10.round(1).astype(str))
    heat(axes[r, 2], score, "hits_in_top10", f"{lab}: acceptable-or-better in top 10 (random expectation)", 0, 10, "Blues",
         lambda v: f"{int(v)}")
    # add the random expectation under each count
    m = top.pivot(index="weight", columns="target", values="random_top10").reindex(index=WEIGHTS, columns=TARGETS)
    hits = top.pivot(index="weight", columns="target", values="hits_in_top10").reindex(index=WEIGHTS, columns=TARGETS)
    for (i, j), v in np.ndenumerate(m.values):
        if not np.isnan(v):
            axes[r, 2].text(j, i + 0.3, f"({v:.1f})", ha="center", va="center", fontsize=6.5,
                            color="#e6eefa" if dark_cell(hits.values[i, j] / 10, "Blues") else MUTED)
fig.suptitle("AF2Rank ranking of CAPRI Scorer decoys per target and weight (all 9 targets × 17 weights)",
             y=1.0, fontsize=12, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.985))
fig.savefig(OUT / "fullres_ranking_summary.png", dpi=170, bbox_inches="tight")
plt.close(fig)

scatter_grid("i_ptm", "ipTM")
scatter_grid("actifptm", "actifpTM")

# Console summary
for score in ["i_ptm", "actifptm"]:
    s = summ[summ.score == score]
    print(score, "median ρ by weight:", s.groupby("weight").spearman_dockq.median().round(2).to_dict())
    print(score, "median AUC by weight:", s.groupby("weight").auc_hit.median().round(2).to_dict())
    print(score, "top-1 acceptable+ :", (s.top1_class != "incorrect").sum(), "/", len(s))
print(summ[summ.score == "i_ptm"].pivot(index="target", columns="weight", values="spearman_dockq").round(2))
