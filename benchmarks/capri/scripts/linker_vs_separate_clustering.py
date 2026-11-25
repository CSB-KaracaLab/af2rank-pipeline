"""HADDOCK-style clustering for the masked multimer runs in both input modes: actifpTM with the 20-glycine linker,
with separate chains, and their product, each ranked alone and with clustering; plus mixed top-5 lists.

Clustering as in clustering.py (FCC >= 0.6 between input decoys, greedy, clusters of >= 4, ranked by the mean score of
their best 4); the clusters come from the input decoys, so they are the same in every mode.
Top-5 lists: 'alone' (top 5 by score), 'clustered' (best decoy of the 5 best clusters), and mixed:
'C4+A1' (4 best clusters, then the top decoy by score), 'C3+A2', 'interleave' (cluster 1, score 1, cluster 2, ...).
Picks are de-duplicated; a list that runs short is filled from the next clusters, then by score.
Top-1 rules: alone, clustered, and 'lone-if-ahead x' (the top decoy by score instead of the best cluster's pick when it
is in no cluster and scores at least x times that cluster's score).
Writes results/linker_vs_separate_clustering_{per_run,one_ranking}.csv and two figures:
  figures/linker_vs_separate_clustering_per_target.png   top-1 and top-5 per target, alone vs clustered
  figures/linker_vs_separate_clustering_scatter.png      score vs CAPRI DockQ per target, alone and clustered
Usage (from benchmarks/capri): python scripts/linker_vs_separate_clustering.py [cutoff]
"""
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import HERE, load_pair, targets  # noqa: E402

CUTOFF = float(sys.argv[1]) if len(sys.argv) > 1 else 0.6
d = load_pair()
MODES = {"linker": "actifptm_linker", "separate": "actifptm_separate", "product": "product"}
NAMES = {"linker": "20-glycine linker", "separate": "separate chains", "product": "linker × separate"}
VERSIONS = ["AF2.1", "AF2.2", "AF2.3"]
RATIOS = [1.1, 1.2, 1.3, 1.5]


def clusters(target, cutoff=CUTOFF, min_size=4):
    """Cluster label per decoy id; -1 = not in a cluster of at least min_size (same as clustering.py)."""
    z = np.load(HERE / f"results/fcc/fcc_{target}.npz"); f = z["fcc"]
    nb = np.minimum(f, f.T) >= cutoff; np.fill_diagonal(nb, False)
    empty = z["ncontacts"] == 0; nb[empty] = False; nb[:, empty] = False
    lab = -np.ones(len(f), int); k = 0
    while (lab < 0).any():
        free = lab < 0; cnt = (nb & free[None, :]).sum(1); cnt[~free] = -1
        c = int(np.argmax(cnt)); lab[c] = k; lab[np.where(nb[c] & free)[0]] = k; k += 1
    s = pd.Series(lab, index=z["ids"])
    return s.where(s.map(s.value_counts()) >= min_size, -1)


CL = {t: clusters(t) for t in d.target.unique()}


def lists(g, score):
    """Top-1 picks and top-5 lists for one target's decoys."""
    g = g.assign(cluster=g.model_id.map(CL[g.target.iloc[0]]).values); gc = g[g.cluster >= 0]
    rank = gc.groupby("cluster")[score].apply(lambda x: x.nlargest(4).mean()).sort_values(ascending=False)
    C = [gc[gc.cluster == c].sort_values(score, ascending=False).iloc[0] for c in rank.index]
    A = [r for _, r in g.sort_values(score, ascending=False).iterrows()]

    def fill(seq):
        out, seen = [], set()
        for r in seq + C + A:
            if r.model_id not in seen:
                out.append(r); seen.add(r.model_id)
            if len(out) == 5:
                return out
        return out
    top5 = {"alone": A[:5], "clustered": fill(C[:5]), "C4+A1": fill(C[:4] + A[:1]), "C3+A2": fill(C[:3] + A[:2]),
            "interleave": fill([x for pair in zip(C[:3], A[:3]) for x in pair][:5])}
    top1 = {"alone": A[0], "clustered": C[0]}
    for x in RATIOS:
        top1[f"lone-if-ahead {x}"] = A[0] if (A[0].cluster < 0 and A[0][score] >= x * rank.iloc[0]) else C[0]
    return top1, top5


rows = []
for (t, w), g in d.groupby(["target", "weight"]):
    for m, c in MODES.items():
        t1, t5 = lists(g, c)
        for k, r in t1.items():
            rows.append(dict(target=t, weight=w, version=g.version.iloc[0], mode=m, list="top1", strategy=k,
                             hit=bool(r.good), dockq=r.dockq, first_class=r.main_classification))
        for k, lst in t5.items():
            rows.append(dict(target=t, weight=w, version=g.version.iloc[0], mode=m, list="top5", strategy=k,
                             hit=any(r.good for r in lst), dockq=max(r.dockq for r in lst), first_class=lst[0].main_classification))
r = pd.DataFrame(rows); r.to_csv(HERE / "results/linker_vs_separate_clustering_per_run.csv", index=False)

one = []
for v, sub in [(v, d[d.version == v]) for v in VERSIONS] + [("all multimer", d)]:
    for t, g in sub.groupby("target"):
        for m, c in MODES.items():
            av = g.groupby("model_id").agg(s=(c, "mean"), good=("good", "first"), dockq=("dockq", "first"),
                                           main_classification=("main_classification", "first")).reset_index().assign(target=t)
            t1, t5 = lists(av, "s")
            for k in ["alone", "clustered"]:
                one.append(dict(version=v, target=t, mode=m, strategy=k, top1=bool(t1[k].good),
                                top5=any(x.good for x in t5[k]), first_class=t1[k].main_classification))
one = pd.DataFrame(one); one.to_csv(HERE / "results/linker_vs_separate_clustering_one_ranking.csv", index=False)

pd.set_option("display.width", 220)
for lst in ["top1", "top5"]:
    q = r[r.list == lst]; order = list(dict.fromkeys(q.strategy))
    s = q.groupby(["strategy", "mode"]).agg(hit=("hit", "mean"), dq=("dockq", "mean"))
    print(f"\n== {lst} per run: % of {q.groupby(['target', 'weight']).ngroups} runs "
          f"(mean CAPRI DockQ of {'the pick' if lst == 'top1' else 'the best of the five'}) ==")
    print(s.apply(lambda x: f"{100 * x.hit:3.0f}% ({x.dq:.2f})", axis=1).unstack("mode").reindex(order)[list(MODES)].to_string())
    print("-- by version, product --")
    p = q[q["mode"] == "product"].groupby(["strategy", "version"]).hit.mean().unstack("version").reindex(order)
    print((100 * p).round(0).astype(int).to_string())
q = r[(r.list == "top1") & r.strategy.isin(["alone", "clustered"])]
print("\n== top-1 by version, alone -> clustered ==")
print((100 * q.groupby(["version", "mode", "strategy"]).hit.mean()).round(0).astype(int).unstack(["mode", "strategy"]).to_string())
print("\n== top-1 per target (% of 15 runs) ==")
print((100 * q.groupby(["target", "mode", "strategy"]).hit.mean()).round(0).astype(int).unstack(["mode", "strategy"]).to_string())
nt = q[q.target != "T131_1"].groupby(["mode", "strategy"]).hit.agg(["sum", "size"])
print("\n== top-1 without T131 ==\n", nt.to_string())
print("\n== one ranking per version (targets of 9, top-1/top-5) ==")
oc = one.groupby(["version", "mode", "strategy"])[["top1", "top5"]].sum()
print((oc.top1.astype(str) + "/" + oc.top5.astype(str)).unstack(["mode", "strategy"]).to_string())

# ---- figure 1: top-1 and top-5 per target ----
tlist = sorted(d.target.unique())
tg = d.drop_duplicates(["target", "model_id"]).groupby("target").good.agg(["sum", "size"])
rnd = {"top1": tg["sum"] / tg["size"],
       "top5": 1 - tg.apply(lambda x: math.comb(int(x["size"] - x["sum"]), 5) / math.comb(int(x["size"]), 5), axis=1)}
SHADES = {"linker": ("#a9c8ef", "#1c5cab"), "separate": ("#f3c79a", "#b85f12"), "product": ("#cdc3ef", "#5a3fb0")}
INK, MUTED, RING = "#0b0b0b", "#5f5e5a", "#d1495b"
plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": "#c3c2b7", "xtick.color": MUTED, "ytick.color": MUTED})
fig, axes = plt.subplots(2, 1, figsize=(14, 8.2), dpi=200, sharex=True)
for ax, k, lab in [(axes[0], "top1", "top-1: good decoy ranked first"), (axes[1], "top5", "top-5: good decoy in the first five")]:
    x = np.arange(len(tlist)); w = 0.135
    for i, m in enumerate(MODES):
        for j, meth in enumerate(["alone", "clustered"]):
            q = r[(r.list == k) & (r["mode"] == m) & (r.strategy == meth)]
            ax.bar(x + (2 * i + j - 2.5) * w, q.groupby("target").hit.mean().reindex(tlist), w * 0.9, color=SHADES[m][j],
                   zorder=2, label=f"{NAMES[m]}, {meth}: {q.hit.mean():.0%} of runs")
    ax.scatter(x, [rnd[k][t] for t in tlist], marker="_", s=900, color=INK, lw=1.6, zorder=3,
               label=f"random picking: {d.drop_duplicates(['target', 'weight']).target.map(rnd[k]).mean():.0%}")
    ax.set_ylim(0, 1.05); ax.set_ylabel("share of the 15 multimer weights")
    ax.grid(axis="y", color="#e4e3df", lw=0.6, zorder=0)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.set_title(lab, loc="left", fontsize=10, color=INK)
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8.3)
axes[1].set_xticks(np.arange(len(tlist)), [t[:4] for t in tlist])
axes[0].text(tlist.index("T131_1"), 0.93, "6 good decoys,\nno cluster of 4", ha="center", fontsize=8, color=MUTED)
fig.suptitle("Masked CAPRI, 15 multimer weights: score alone (light) vs HADDOCK-style cluster ranking (dark; "
             f"FCC ≥ {CUTOFF}, clusters of ≥ 4, ranked by the mean score of their best 4)", x=0.01, ha="left", fontsize=10.5, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(HERE / "figures/linker_vs_separate_clustering_per_target.png", bbox_inches="tight")
plt.rcdefaults()

# ---- figure 2: score (mean over weights) vs CAPRI DockQ per target, alone and clustered ----
pr1 = r[r.list == "top1"].groupby(["mode", "strategy", "target"]).hit.mean()
m_ = d.groupby(["target", "model_id"]).agg(**{c: (c, "mean") for c in MODES.values()}, dockq=("dockq", "first"),
                                           cls=("main_classification", "first"), good=("good", "first")).reset_index()
sel = targets()
CLASSES = [("incorrect", "#c9c7c1", 9, 0.6), ("acceptable", "#86b6ef", 15, 0.9), ("medium", "#2a78d6", 15, 0.9),
           ("high", "#0d366b", 17, 0.95)]
fig, ax_all = plt.subplots(len(tlist), 8, figsize=(19, 2.45 * len(tlist) + 1.4), dpi=160, sharex=True,
                           gridspec_kw={"width_ratios": [1, 1, 0.16, 1, 1, 0.16, 1, 1]})
for a in ax_all[:, [2, 5]].flat: a.set_visible(False)
cols = [0, 1, 3, 4, 6, 7]
for ri, t in enumerate(tlist):
    s = m_[m_.target == t].copy(); s["cluster"] = s.model_id.map(CL[t]).values
    for k, (mode, c) in enumerate(MODES.items()):
        for j, meth in enumerate(["alone", "clustered"]):
            ax = ax_all[ri, cols[2 * k + j]]; ax.grid(color="#e4e3df", lw=0.6, zorder=0)
            for cl, col, sz, a in CLASSES:
                h_ = s[s.cls == cl]
                if meth == "alone":
                    ax.scatter(h_.dockq, h_[c], s=sz, c=col, alpha=a, lw=0, zorder=2 if cl == "incorrect" else 3)
                else:
                    ax.scatter(h_[h_.cluster < 0].dockq, h_[h_.cluster < 0][c], s=sz, c=col, alpha=0.18, lw=0, zorder=2)
                    ax.scatter(h_[h_.cluster >= 0].dockq, h_[h_.cluster >= 0][c], s=sz, c=col, alpha=a, lw=0, zorder=3)
            if meth == "alone":
                top = s.loc[s[c].idxmax()]
            else:
                gc = s[s.cluster >= 0]
                rank = gc.groupby("cluster")[c].apply(lambda x: x.nlargest(4).mean()).sort_values(ascending=False)
                best = gc[gc.cluster == rank.index[0]]
                ax.scatter(best.dockq, best[c], s=26, facecolors="none", edgecolors=INK, lw=0.6, zorder=4)
                top = best.loc[best[c].idxmax()]
                ax.set_facecolor("#f6f5f2")
            ax.scatter([top.dockq], [top[c]], s=120, facecolors="none", edgecolors=RING, lw=1.7, zorder=5)
            ax.text(0.04, 0.97, f"first: {top.cls}\nper run: top-1 {pr1[(mode, meth, t)]:.0%}", ha="left", va="top",
                    transform=ax.transAxes, fontsize=7.2, color=INK,
                    bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="none", alpha=0.88), zorder=6)
            ax.set_xlim(-0.03, 1); ax.set_ylim(0, 1.0)
            for side in ("top", "right"): ax.spines[side].set_visible(False)
            ax.tick_params(labelsize=6.8, colors=MUTED)
            if cols[2 * k + j] != 0: ax.set_yticklabels([])
            if ri == 0:
                ax.set_title(f"{NAMES[mode]}\n{'score alone' if meth == 'alone' else 'clustered'}", fontsize=9, color=INK,
                             fontweight="bold" if meth == "clustered" else "normal")
            if ri == len(tlist) - 1: ax.set_xlabel("CAPRI DockQ", fontsize=8)
    info = sel.loc[t.replace("_", ".")]; desc = info.desc.split("(")[0].strip()
    desc = desc if len(desc) <= 30 else desc[:30].rsplit(" ", 1)[0] + "…"
    ax_all[ri, 0].set_ylabel(f"{t[:4]} · {info.pdb}\n{desc}\n{int(s.good.sum())}/{len(s)} good · "
                             f"{int((s.cluster >= 0).sum())} in {s[s.cluster >= 0].cluster.nunique()} clusters\n\n"
                             "score, mean of 15 weights", fontsize=7.8, color=INK)
h = [plt.Line2D([], [], ls="", marker="o", ms=6, mfc=c, mec="none", label=nm) for nm, c, _, _ in CLASSES]
h += [plt.Line2D([], [], ls="", marker="o", ms=6, mfc="#86b6ef", mec="none", alpha=0.25, label="not in a cluster of ≥ 4 (clustered panels)"),
      plt.Line2D([], [], ls="", marker="o", ms=5, mfc="none", mec=INK, mew=0.7, label="members of the top cluster"),
      plt.Line2D([], [], ls="", marker="o", ms=10, mfc="none", mec=RING, mew=1.7, label="ranked first")]
fig.legend(handles=h, loc="upper center", ncol=7, frameon=False, bbox_to_anchor=(0.5, 0.985), fontsize=8.8)
fig.suptitle("Masked CAPRI, 15 multimer weights: actifpTM with the linker, with separate chains, and their product, ranked "
             f"alone and with HADDOCK-style clustering (FCC ≥ {CUTOFF}, clusters of ≥ 4, ranked by the mean score of their best 4)",
             fontsize=11, color=INK, y=0.997)
fig.subplots_adjust(left=0.085, right=0.995, top=0.945, bottom=0.03, hspace=0.18, wspace=0.07)
fig.savefig(HERE / "figures/linker_vs_separate_clustering_scatter.png")
print("\nwrote results/linker_vs_separate_clustering_*.csv and two figures")
