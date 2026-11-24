"""Masked (AF2Rank paper setting) vs full resolution on the same CAPRI runs.

Masked = the AF2Rank paper's setting (sequence and side chains masked, 1 recycle); full resolution = 0 recycles.
All 153 runs (9 targets x 17 weights) in both settings.
Writes results/masked_vs_fullres_summary.csv (per-run medians of AUC, top-1, good in top 10 for each score),
results/masked_vs_fullres_by_family.csv (the same per weight family), results/masked_vs_fullres_consensus.csv
(one ranking per target from each decoy's rank averaged over a set of weights), results/masked_vs_fullres_per_target.csv
(median AUC per target and family), and figures/masked_vs_fullres_<score>_<weight>.png (score vs CAPRI DockQ, one row per setting, one column per target).
Usage (from benchmarks/capri): python scripts/masked_vs_fullres.py [score] [weight ...]
"""
import math
import sys
from pathlib import Path

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import HERE, load, targets  # noqa: E402

SCORE = sys.argv[1] if len(sys.argv) > 1 else "actifptm"
WEIGHTS = sys.argv[2:] or ["model_1_ptm", "model_1_multimer", "model_1_multimer_v2", "model_1_multimer_v3"]
FAM = {"ptm": "AF2-ptm", "multimer": "AF2.1", "multimer_v2": "AF2.2", "multimer_v3": "AF2.3"}
SCORES = ["actifptm", "i_ptm", "neg_pae_contacts", "dockq_io", "af2rank", "composite"]
LABEL = {"actifptm": "actifpTM", "i_ptm": "ipTM", "af2rank": "pLDDT × pTM × DockQ_io"}.get(SCORE, SCORE)
m = load("masked").assign(setting="masked")
f = load("fullres").assign(setting="full")
runs = set(zip(m.weight, m.target))
f = f[[(w, t) in runs for w, t in zip(f.weight, f.target)]]
D = pd.concat([m, f], ignore_index=True)
D["af2rank"] = D.plddt * D.ptm * D.dockq_io
z = lambda c: D.groupby(["setting", "target", "weight"])[c].transform(lambda x: (x - x.mean()) / (x.std() + 1e-9))
D["composite"] = z("neg_pae_contacts") + z("actifptm") + z("fnat_io")

# per-run summary
rows = []
for (s, t, w), g in D.groupby(["setting", "target", "weight"]):
    if not 0 < g.good.sum() < len(g):
        continue
    for k in SCORES:
        o = g.sort_values(k, ascending=False)
        rows.append(dict(setting=s, target=t, weight=w, family=FAM[g.family.iloc[0]], score=k, auc=roc_auc_score(g.good, g[k]),
                         top1=o.good.iloc[0], top5=o.good.head(5).any(), top10=o.good.head(10).sum(),
                         random_top1=g.good.mean(), random_top10=10 * g.good.mean(),
                         random_top5=1 - math.comb(int((~g.good).sum()), 5) / math.comb(len(g), 5)))
r = pd.DataFrame(rows)
summ = r.groupby(["score", "setting"]).agg(median_auc=("auc", "median"), top1_success=("top1", "mean"),
                                           top5_success=("top5", "mean"), good_in_top10=("top10", "mean"),
                                           random_top1=("random_top1", "mean"), random_top5=("random_top5", "mean"),
                                           random_top10=("random_top10", "mean"), runs=("auc", "size")).round(3)
summ.to_csv(HERE / "results/masked_vs_fullres_summary.csv")
print(len(runs), "runs in both settings\n", summ.unstack("setting").to_string())
D["fam"] = D.family.map(FAM)
print("\nmedian DockQ_io by CAPRI class:\n",
      D.pivot_table(index="main_classification", columns=["fam", "setting"], values="dockq_io", aggfunc="median")
       .reindex(["incorrect", "acceptable", "medium", "high"]).round(2).to_string())

# per weight family
fam = r.groupby(["family", "score", "setting"]).agg(median_auc=("auc", "median"), top1_success=("top1", "mean"),
                                                     top5_success=("top5", "mean"), good_in_top10=("top10", "mean"),
                                                     runs=("auc", "size")).round(3)
fam.to_csv(HERE / "results/masked_vs_fullres_by_family.csv")
print("\nby family (median AUC / top-1 / top-5):\n",
      fam.unstack("setting")[["median_auc", "top1_success", "top5_success"]].to_string())

# one ranking per target: each decoy's percentile rank averaged over a set of weights
SETS = {"all 17 weights": list(FAM.values()), **{v: [v] for v in FAM.values()}}
crow = []
for name, fams in SETS.items():
    for s in ["full", "masked"]:
        g0 = D[(D.setting == s) & D.fam.isin(fams)]
        for k in ["actifptm", "composite", "neg_pae_contacts", "i_ptm"]:
            rk = g0.assign(rk=g0.groupby(["target", "weight"])[k].rank(pct=True))
            c = rk.groupby(["target", "model_id"]).agg(rk=("rk", "mean"), good=("good", "first"), dockq=("dockq", "first")).reset_index()
            for t, g in c.groupby("target"):
                o = g.sort_values("rk", ascending=False)
                crow.append(dict(weights=name, setting=s, score=k, target=t, auc=roc_auc_score(g.good, g.rk),
                                 top1=bool(o.good.iloc[0]), top5=bool(o.good.head(5).any()), top1_dockq=o.dockq.iloc[0],
                                 top10=int(o.good.head(10).sum())))
cons = pd.DataFrame(crow); cons.round(3).to_csv(HERE / "results/masked_vs_fullres_consensus.csv", index=False)
print("\nconsensus (targets with a good decoy first, median AUC):\n",
      cons.groupby(["weights", "score", "setting"]).agg(top1=("top1", "sum"), top5=("top5", "sum"),
                                                        median_auc=("auc", "median"),
                                                        good_top10=("top10", "mean")).unstack("setting").round(2).to_string())
pt = (r[r.score == SCORE].groupby(["target", "family", "setting"]).auc.median().unstack(["family", "setting"]).round(2))
pt.to_csv(HERE / "results/masked_vs_fullres_per_target.csv")
print(f"\nper target, median {SCORE} AUC:\n", pt.to_string())
sel = targets()
for WEIGHT in WEIGHTS:
    d = D[D.weight == WEIGHT]
    TARGETS = sorted(d.target.unique())
    SETTINGS = [("full", "full resolution\n0 recycles"), ("masked", "masked seq. + side chains\n1 recycle (paper)")]
    CLASSES = [("incorrect", "#c9c7c1", 10), ("acceptable", "#86b6ef", 16), ("medium", "#2a78d6", 16), ("high", "#0d366b", 16)]
    INK, MUTED = "#0b0b0b", "#898781"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c3c2b7", "xtick.color": MUTED, "ytick.color": MUTED,
                         "axes.labelcolor": "#52514e", "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})

    fig, axes = plt.subplots(2, len(TARGETS), figsize=(2.05 * len(TARGETS), 4.6), sharex=True, sharey=True)
    for i, (s, slab) in enumerate(SETTINGS):
        for j, t in enumerate(TARGETS):
            ax = axes[i, j]; g = d[(d.setting == s) & (d.target == t)]
            ax.grid(color="#e1e0d9", lw=0.6); ax.set_axisbelow(True)
            for cls, col, sz in CLASSES:
                h = g[g.main_classification == cls]
                ax.scatter(h.dockq, h[SCORE], s=sz, c=col, lw=0.5 if cls != "incorrect" else 0, edgecolors="#fcfcfb",
                           zorder=3 if cls != "incorrect" else 2)
            rho = spearmanr(g.dockq, g[SCORE])[0]
            auc = roc_auc_score(g.good, g[SCORE]) if 0 < g.good.sum() < len(g) else float("nan")
            ax.text(0.97, 0.05, f"ρ {rho:.2f}\nAUC {auc:.2f}", ha="right", va="bottom", transform=ax.transAxes,
                    fontsize=7.5, color=INK, bbox=dict(boxstyle="round,pad=0.2", fc="#fcfcfb", ec="none", alpha=0.85), zorder=5)
            ax.set_xlim(-0.03, 1.0); ax.set_ylim(0, 1.02)
            if i == 0:
                r = sel.loc[t.replace("_", ".")]
                ax.set_title(f"{t.replace('_', '.')} ({r.pdb})\n{r.desc.split(' ')[0][:20]}", fontsize=8.5, color=INK)
            if j == 0:
                ax.set_ylabel(f"{slab}\n{LABEL}")
            if i == 1:
                ax.set_xlabel("CAPRI DockQ")
    handles = [plt.Line2D([], [], ls="", marker="o", ms=6, mfc=c, mec="none", label=n) for n, c, _ in CLASSES]
    fig.legend(handles=handles, title="CAPRI class", loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.0))
    wl = WEIGHT.replace("model_", "").replace("_ptm", " pTM").replace("_multimer_v2", " AF2.2 multimer").replace("_multimer_v3", " AF2.3 multimer").replace("_multimer", " AF2.1 multimer")
    fig.suptitle(f"AF2Rank {LABEL} vs CAPRI DockQ, CAPRI Scorer decoys, weight: model {wl}", y=1.06, fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = HERE / f"figures/masked_vs_fullres_{SCORE}_{WEIGHT}.png"
    fig.savefig(out, dpi=170, bbox_inches="tight"); print("wrote", out)
