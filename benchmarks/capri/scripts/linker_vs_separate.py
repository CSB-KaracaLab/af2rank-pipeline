"""Masked multimer runs: 20-glycine linker vs separate chains, and their product.

For each run (one target scored with one multimer weight, 135 per mode), ranks the decoys by actifpTM with the linker,
actifpTM with separate chains, and their product (computed per decoy and weight). Reports AUC (good = CAPRI acceptable
or better vs incorrect), top-1 and top-5 per run, one ranking per version (score averaged over the version's 5 weights),
decoys pulled apart (DockQ input->output < 0.23) and the Spearman correlation between the two modes.
Writes results/linker_vs_separate_*.csv and two figures:
  figures/linker_vs_separate_topk_grid.png   class of the decoy ranked first, per target and weight
  figures/linker_vs_separate_by_target.png   score (mean over weights) vs CAPRI DockQ per target
Usage (from benchmarks/capri): python scripts/linker_vs_separate.py
"""
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import HERE, load_pair, targets  # noqa: E402

d = load_pair()
MODES = {"linker": "actifptm_linker", "separate": "actifptm_separate", "product": "product"}
NAMES = {"linker": "20-glycine linker", "separate": "separate chains", "product": "linker × separate"}
VERSIONS = ["AF2.1", "AF2.2", "AF2.3"]
tg = d.drop_duplicates(["target", "model_id"]).groupby("target").good.agg(["sum", "size"])
RAND1 = tg["sum"] / tg["size"]
RAND5 = 1 - tg.apply(lambda r: math.comb(int(r["size"] - r["sum"]), 5) / math.comb(int(r["size"]), 5), axis=1)

# ---- per run ----
rows = []
for (t, w), g in d.groupby(["target", "weight"]):
    rho = spearmanr(g.actifptm_linker, g.actifptm_separate)[0]
    for m, c in MODES.items():
        o = g.sort_values(c, ascending=False)
        rows.append(dict(target=t, weight=w, version=g.version.iloc[0], mode=m, first_class=o.main_classification.iloc[0],
                         top1=bool(o.good.iloc[0]), top5=bool(o.good.head(5).any()),
                         auc=roc_auc_score(g.good, g[c]) if 0 < g.good.sum() < len(g) else np.nan,
                         spearman_linker_separate=rho))
pr = pd.DataFrame(rows)
pr.to_csv(HERE / "results/linker_vs_separate_per_run.csv", index=False)


def one_ranking(sub, c):
    """Score averaged over the weights in sub, ranked once per target: (top-1, top-5) counts."""
    m = sub.groupby(["target", "model_id"]).agg(s=(c, "mean"), good=("good", "first")).reset_index()
    o = m.sort_values("s", ascending=False).groupby("target")
    return int(o.good.first().sum()), int(o.good.apply(lambda x: x.head(5).any()).sum())


ver = []
for v, sub in [(v, d[d.version == v]) for v in VERSIONS] + [("all multimer", d)]:
    p = pr[pr.version.isin(sub.version.unique())]
    for m, c in MODES.items():
        q = p[p["mode"] == m]; a1, a5 = one_ranking(sub, c)
        ver.append(dict(version=v, mode=m, runs=len(q), auc=q.auc.median(), top1=q.top1.mean(), top5=q.top5.mean(),
                        one_ranking_top1=a1, one_ranking_top5=a5, random_top1=q.target.map(RAND1).mean(),
                        random_top5=q.target.map(RAND5).mean()))
ver = pd.DataFrame(ver).round(3)
ver.to_csv(HERE / "results/linker_vs_separate_by_version.csv", index=False)
agg = dict(runs=("top1", "size"), top1=("top1", "mean"), top5=("top5", "mean"), auc=("auc", "median"))
bt = pr.groupby(["target", "mode"]).agg(**agg).unstack("mode").round(3)
bt.to_csv(HERE / "results/linker_vs_separate_by_target.csv")
bw = pr.groupby(["version", "weight", "mode"]).agg(**agg).unstack("mode").round(3)
bw.to_csv(HERE / "results/linker_vs_separate_by_weight.csv")

# ---- pulled apart and agreement between the modes ----
ap = []
for v, sub in [(v, d[d.version == v]) for v in VERSIONS] + [("all multimer", d)]:
    for good, lab in [(True, "good"), (False, "incorrect")]:
        s = sub[sub.good == good]
        ap.append(dict(version=v, decoys=lab, n=len(s), pulled_apart_linker=(s.dockq_io_linker < 0.23).mean(),
                       pulled_apart_separate=(s.dockq_io_separate < 0.23).mean()))
ap = pd.DataFrame(ap).round(3)
ap.to_csv(HERE / "results/linker_vs_separate_pulled_apart.csv", index=False)
rho = pr.drop_duplicates(["target", "weight"])

pd.set_option("display.width", 200)
print("== by version (per run: median AUC, top-1, top-5; one ranking per version: targets of 9) ==")
print(ver.to_string(index=False))
print("\n== per target, top-1 per run ==\n", bt["top1"].to_string())
print("\n== decoys pulled apart (DockQ input->output < 0.23) ==\n", ap.to_string(index=False))
print("\n== median Spearman between linker and separate actifpTM, per run ==")
print(rho.groupby("version").spearman_linker_separate.median().round(2).to_string(),
      "\nall", round(rho.spearman_linker_separate.median(), 2))

# ---- figure 1: class of the decoy ranked first, per target x weight ----
COL = {"incorrect": "#dddbd6", "acceptable": "#86b6ef", "medium": "#2a78d6", "high": "#0d366b"}
INK, MUTED = "#0b0b0b", "#5f5e5a"
tlist = sorted(d.target.unique())
wlist = sorted(d.weight.unique(), key=lambda w: (VERSIONS.index(d[d.weight == w].version.iloc[0]), w))
q = pr.set_index(["mode", "weight", "target"])
fig, axes = plt.subplots(1, 3, figsize=(23, 8.4), dpi=170)
for ax, m in zip(axes, MODES):
    for i, w in enumerate(wlist):
        y = len(wlist) - 1 - i
        for j, t in enumerate(tlist):
            c = q.loc[(m, w, t)]
            ax.add_patch(plt.Rectangle((j, y), 1, 1, fc=COL[c.first_class], ec="white", lw=1.5))
            if c.top5:
                ax.plot(j + 0.5, y + 0.5, "o", ms=4, color="white" if c.first_class in ("medium", "high") else INK)
        r = pr[(pr["mode"] == m) & (pr.weight == w)]
        ax.text(len(tlist) + 0.15, y + 0.5, f"{int(r.top1.sum())}/9 · {int(r.top5.sum())}/9", va="center", fontsize=8.5, color=INK)
    for j, t in enumerate(tlist):
        r = pr[(pr["mode"] == m) & (pr.target == t)]
        ax.text(j + 0.5, -0.35, f"{int(r.top1.sum())}/15\n{int(r.top5.sum())}/15", ha="center", va="top", fontsize=8, color=INK)
    for k in (5, 10):
        ax.axhline(len(wlist) - k, color=INK, lw=1.2)
    ax.set_xlim(0, len(tlist)); ax.set_ylim(-1.3, len(wlist)); ax.set_aspect("equal")
    ax.set_xticks(np.arange(len(tlist)) + 0.5, [t[:4] for t in tlist]); ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(wlist)) + 0.5, [f"{'AF2.3' if w.endswith('v3') else 'AF2.2' if w.endswith('v2') else 'AF2.1'} model {w[6]}"
                                                for w in wlist][::-1])
    ax.tick_params(length=0, labelsize=9)
    for s in ax.spines.values(): s.set_visible(False)
    p = pr[pr["mode"] == m]
    ax.set_title(f"{NAMES[m]}\ntop-1 {p.top1.mean():.0%} · top-5 {p.top5.mean():.0%} of {len(p)} runs", fontsize=11.5,
                 fontweight="bold" if m == "product" else "normal", pad=26)
    ax.text(len(tlist) + 0.15, len(wlist) + 0.2, "top-1 · top-5", fontsize=8, color=MUTED)
    ax.text(-0.1, -0.35, "top-1\ntop-5", ha="right", va="top", fontsize=8, color=MUTED)
h = [plt.Line2D([], [], ls="", marker="s", ms=10, mfc=c, mec="none", label=f"ranked first: {k}") for k, c in COL.items()]
h += [plt.Line2D([], [], ls="", marker="o", ms=5, mfc=INK, mec="none", label="good decoy in top 5")]
fig.legend(handles=h, loc="lower center", ncol=5, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.0))
fig.suptitle("Masked CAPRI, actifpTM: decoy ranked first per target and multimer weight — linker, separate chains, "
             "and their product", fontsize=13)
fig.subplots_adjust(left=0.05, right=0.97, top=0.86, bottom=0.08, wspace=0.42)
fig.savefig(HERE / "figures/linker_vs_separate_topk_grid.png")

# ---- figure 2: score (mean over weights) vs CAPRI DockQ, per target ----
sel = targets()
m = d.groupby(["target", "model_id"]).agg(**{c: (c, "mean") for c in MODES.values()}, dockq=("dockq", "first"),
                                          cls=("main_classification", "first"), good=("good", "first")).reset_index()
CLASSES = [("incorrect", "#c9c7c1", 9, 0.6), ("acceptable", "#86b6ef", 16, 0.9), ("medium", "#2a78d6", 16, 0.9),
           ("high", "#0d366b", 18, 0.95)]
fig, ax_all = plt.subplots(3, 11, figsize=(24.5, 9.6), dpi=170, sharex=True,
                           gridspec_kw={"width_ratios": [1, 1, 1, 0.22, 1, 1, 1, 0.22, 1, 1, 1]})
for a in ax_all[:, [3, 7]].flat: a.set_visible(False)
axes = ax_all[:, [0, 1, 2, 4, 5, 6, 8, 9, 10]]
for k, t in enumerate(tlist):
    r, c0 = divmod(k, 3); s = m[m.target == t]
    for i, (mode, c) in enumerate(MODES.items()):
        ax = axes[r, 3 * c0 + i]; ax.grid(color="#e4e3df", lw=0.6, zorder=0)
        for cl, col, sz, a in CLASSES:
            h_ = s[s.cls == cl]
            ax.scatter(h_.dockq, h_[c], s=sz, c=col, alpha=a, lw=0, zorder=2 if cl == "incorrect" else 3)
        top = s.loc[s[c].idxmax()]
        ax.scatter([top.dockq], [top[c]], s=110, facecolors="none", edgecolors="#d1495b", lw=1.6, zorder=4)
        p = pr[(pr["mode"] == mode) & (pr.target == t)]
        ax.text(0.04, 0.97, f"AUC {roc_auc_score(s.good, s[c]):.2f} · first: {top.cls}\nper run: AUC {p.auc.median():.2f}, "
                f"top-1 {p.top1.mean():.0%}", ha="left", va="top", transform=ax.transAxes, fontsize=7.3, color=INK,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.88), zorder=5)
        ax.set_xlim(-0.03, 1); ax.set_ylim(0, 1.0)
        if mode == "product": ax.set_facecolor("#f6f5f2")
        for side in ("top", "right"): ax.spines[side].set_visible(False)
        ax.tick_params(labelsize=7, colors=MUTED)
        if i > 0: ax.set_yticklabels([])
        ax.set_title(NAMES[mode], fontsize=8.5, color=INK if mode == "product" else MUTED, pad=3,
                     fontweight="bold" if mode == "product" else "normal")
        if i == 0:
            info = sel.loc[t.replace("_", ".")]; desc = info.desc.split("(")[0].strip()
            desc = desc if len(desc) <= 34 else desc[:34].rsplit(" ", 1)[0] + "…"
            ax.annotate(f"{t[:4]} · {info.pdb} {desc}\n{int(s.good.sum())}/{len(s)} good · 15 weights", xy=(0, 1.11),
                        xycoords="axes fraction", fontsize=9, color=INK, fontweight="bold", annotation_clip=False)
            if c0 == 0: ax.set_ylabel("actifpTM, mean over weights", fontsize=8.5)
        if r == 2: ax.set_xlabel("CAPRI DockQ", fontsize=8)
h = [plt.Line2D([], [], ls="", marker="o", ms=6, mfc=c, mec="none", label=nm) for nm, c, _, _ in CLASSES]
h.append(plt.Line2D([], [], ls="", marker="o", ms=9, mfc="none", mec="#d1495b", mew=1.6, label="ranked first"))
fig.legend(handles=h, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.955), fontsize=9)
fig.suptitle("Masked CAPRI per target: actifpTM with the linker, with separate chains, and their product vs CAPRI DockQ "
             "(one point per decoy, mean over the 15 multimer weights)", fontsize=11.5, color=INK, y=0.99)
fig.subplots_adjust(left=0.035, right=0.995, top=0.855, bottom=0.06, hspace=0.55, wspace=0.08)
fig.savefig(HERE / "figures/linker_vs_separate_by_target.png")
print("\nwrote results/linker_vs_separate_*.csv, figures/linker_vs_separate_topk_grid.png, "
      "figures/linker_vs_separate_by_target.png")
