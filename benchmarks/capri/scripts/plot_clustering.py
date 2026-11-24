"""Figures for the clustering analysis (run scripts/clustering.py first).

figures/clustering_per_target.png: top-1 and top-5 success per target (share of the 17 weights), actifpTM alone vs
  HADDOCK-style cluster ranking, with random picking.
figures/clustering_scatter_<weight>.png: masked actifpTM vs CAPRI DockQ per target for one weight; decoys in clusters of
  >= 4 solid, others faded; rings mark each method's top-1. Default: model 1 of each AlphaFold version.
Usage (from benchmarks/capri): python scripts/plot_clustering.py [weight ...]
"""
import sys
from pathlib import Path

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import HERE, load, targets  # noqa: E402
from clustering import CL  # noqa: E402  (clusters at the default cutoff)

INK, MUTED = "#0b0b0b", "#76746e"
plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": "#c3c2b7", "xtick.color": MUTED, "ytick.color": MUTED})
r = pd.read_csv(HERE / "results/clustering_per_run.csv")
pt = r.pivot_table(index="target", columns="method", values=["top1", "top5", "random_top1", "random_top5"])
fig, axes = plt.subplots(1, 2, figsize=(14, 4.6), dpi=200, sharey=True)
for ax, k in zip(axes, ["top1", "top5"]):
    x = np.arange(len(pt)); w = 0.38
    for off, m, col in [(-w / 2, "actifpTM alone", "#b5b3ad"), (w / 2, "clustered", "#1c5cab")]:
        lab = "actifpTM alone" if m == "actifpTM alone" else "HADDOCK-style clusters, ranked by actifpTM"
        ax.bar(x + off, pt[(k, m)], w * 0.92, color=col, zorder=2, label=f"{lab}: {r[r.method == m][k].mean():.0%} of runs")
    ax.scatter(x, pt[("random_" + k, "clustered")], marker="_", s=380, color=INK, lw=1.6, zorder=3,
               label=f"random picking: {r['random_' + k].mean():.0%}")
    ax.set_xticks(x, [t.replace("_1", "") for t in pt.index]); ax.set_ylim(0, 1.05)
    ax.grid(axis="y", color="#e4e3df", lw=0.6, zorder=0)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.set_title(f"{'top-1: good decoy ranked first' if k == 'top1' else 'top-5: good decoy in the first five'}",
                 loc="left", fontsize=10, color=INK)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.1), fontsize=8.5)
axes[0].set_ylabel("share of the 17 weights")
fig.suptitle("Masked CAPRI decoys: HADDOCK-style cluster ranking (FCC ≥ 0.6, clusters of ≥ 4, ranked by the mean "
             "actifpTM of their best 4)", x=0.01, ha="left", fontsize=10.5, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(HERE / "figures/clustering_per_target.png"); plt.close(fig)

d = load("masked"); sel = targets()
CLASSES = [("incorrect", "#c9c7c1", 14), ("acceptable", "#86b6ef", 22), ("medium", "#2a78d6", 22), ("high", "#0d366b", 26)]
for W in sys.argv[1:] or ["model_1_ptm", "model_1_multimer", "model_1_multimer_v2", "model_1_multimer_v3"]:
    fig, axes = plt.subplots(2, 5, figsize=(17, 7.6), dpi=200, sharex=True, sharey=True)
    for ax, t in zip(axes.flat, sorted(d.target.unique())):
        g = d[(d.target == t) & (d.weight == W)].copy(); g["cluster"] = g.model_id.map(CL[t]).values
        gc = g[g.cluster >= 0]
        rank = gc.groupby("cluster").actifptm.apply(lambda x: x.nlargest(4).mean()).sort_values(ascending=False)
        ax.grid(color="#e4e3df", lw=0.6, zorder=0)
        for c, col, s in CLASSES:
            h = g[g.main_classification == c]
            ax.scatter(h[h.cluster < 0].dockq, h[h.cluster < 0].actifptm, s=s, c=col, alpha=0.25, lw=0, zorder=2)
            ax.scatter(h[h.cluster >= 0].dockq, h[h.cluster >= 0].actifptm, s=s, c=col, alpha=0.95, lw=0, zorder=3)
        best = gc[gc.cluster == rank.index[0]]
        ax.scatter(best.dockq, best.actifptm, s=34, facecolors="none", edgecolors=INK, lw=0.7, zorder=4)
        a = g.loc[g.actifptm.idxmax()]; h1 = best.loc[best.actifptm.idxmax()]
        ax.scatter(a.dockq, a.actifptm, s=170, facecolors="none", edgecolors="#eb6834", lw=1.6, zorder=5)
        ax.scatter(h1.dockq, h1.actifptm, s=260, facecolors="none", edgecolors=INK, lw=1.6, zorder=5)
        rr = sel.loc[t.replace("_", ".")]
        ax.set_title(f"{t.replace('_', '.')} ({rr.pdb}) {rr.desc.split(' ')[0][:18]}\n{int(g.good.sum())} good / {len(g)} · "
                     f"{len(gc)} in {len(rank)} clusters", fontsize=8.5, color=INK)
        ax.set_xlim(-0.03, 1); ax.set_ylim(0, 1.02); ax.tick_params(labelsize=7.5)
        for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax = axes.flat[-1]; ax.axis("off")
    h = [plt.Line2D([], [], ls="", marker="o", ms=6, mfc=c, mec="none", label=n) for n, c, _ in CLASSES]
    h += [plt.Line2D([], [], ls="", marker="o", ms=6, mfc="#86b6ef", mec="none", alpha=0.25, label="not in a cluster of ≥ 4"),
          plt.Line2D([], [], ls="", marker="o", ms=5, mfc="none", mec=INK, mew=0.7, label="members of the top cluster"),
          plt.Line2D([], [], ls="", marker="o", ms=10, mfc="none", mec="#eb6834", mew=1.6, label="top-1, actifpTM alone"),
          plt.Line2D([], [], ls="", marker="o", ms=12, mfc="none", mec=INK, mew=1.6, label="top-1, HADDOCK-style")]
    ax.legend(handles=h, loc="center", frameon=False, fontsize=9)
    for a_ in axes[1, :4]: a_.set_xlabel("CAPRI DockQ", fontsize=8.5)
    for a_ in axes[:, 0]: a_.set_ylabel("actifpTM (masked)", fontsize=9)
    wl = (W.replace("model_", "model ").replace("_multimer_v2", ", AF2.2").replace("_multimer_v3", ", AF2.3")
          .replace("_multimer", ", AF2.1").replace("_ptm", ", AF2-ptm"))
    fig.suptitle(f"CAPRI decoys, {wl} (masked): actifpTM alone vs HADDOCK-style cluster ranking "
                 "(FCC ≥ 0.6, clusters of ≥ 4, ranked by the mean actifpTM of their best 4)", x=0.01, ha="left", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(HERE / f"figures/clustering_scatter_{W}.png"); plt.close(fig)
print("wrote figures/clustering_*.png")
