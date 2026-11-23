"""Can AF2Rank outputs separate good CAPRI decoys from bad ones? Full resolution, all 153 runs.

Reads results/scores_capri_fullres_all_weights.csv. "Good" = CAPRI acceptable or better.
Scores (higher = better; PAE and iRMSD are negated): pLDDT, pTM, ipTM, actifpTM, PAE, PAE at the input contacts,
DockQ / fnat / iRMSD input->output, ipSAE, LIS, pDockQ, min inter-chain PAE, and AF2Rank-style products of
confidence x how well AlphaFold kept the decoy's interface.
Writes results/per_run_metrics.csv, results/summary_by_score.csv, results/pooled_by_score.csv,
and results/loto_logistic.csv.
Usage (from examples/capri_benchmark): python scripts/analyze_scores.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import HERE, load  # noqa: E402

OUT = HERE / "results"
d = load("fullres")
HAVE_IFACE = True
S = {  # name -> (column or expression, label)
    "plddt": ("plddt", "pLDDT"), "ptm": ("ptm", "pTM"), "i_ptm": ("i_ptm", "ipTM"), "actifptm": ("actifptm", "actifpTM"),
    "neg_pae": (lambda x: -x.pae, "-mean PAE"), "neg_pae_contacts": (lambda x: -x.pae_contacts_mean, "-PAE at input contacts"),
    "dockq_io": ("dockq_io", "DockQ in→out"), "fnat_io": ("fnat_io", "fnat in→out"), "neg_irmsd_io": (lambda x: -x.irmsd_io, "-iRMSD in→out"),
    "af2rank_cplx": (lambda x: x.plddt * x.ptm * x.dockq_io, "pLDDT×pTM×DockQ_io"),
    "iptm_x_dockq": (lambda x: x.i_ptm * x.dockq_io, "ipTM×DockQ_io"),
    "actif_x_dockq": (lambda x: x.actifptm * x.dockq_io, "actifpTM×DockQ_io"),
    "actif_x_fnat": (lambda x: x.actifptm * x.fnat_io, "actifpTM×fnat_io"),
}
if HAVE_IFACE:
    S.update({
        "ipsae": ("ipsae", "ipSAE"), "ipsae_min": ("ipsae_min", "ipSAE (min dir.)"), "lis": ("lis", "LIS"),
        "pdockq": ("pdockq", "pDockQ"), "neg_pae_inter_min": (lambda x: -x.pae_inter_min, "-min inter-chain PAE"),
        "ipsae_x_dockq": (lambda x: x.ipsae * x.dockq_io, "ipSAE×DockQ_io"),
        "lis_x_dockq": (lambda x: x.lis * x.dockq_io, "LIS×DockQ_io"),
    })
for k, (e, _) in S.items():
    d[k] = d[e] if isinstance(e, str) else e(d)
LABEL = {k: v[1] for k, v in S.items()}

# ---- per target x weight ranking ----
rows = []
for (t, w), g in d.groupby(["target", "weight"]):
    for k in S:
        x = g[k]
        if x.isna().all():
            continue
        o = g.sort_values(k, ascending=False)
        rows.append({"target": t, "weight": w, "score": k, "n": len(g), "n_good": int(g.good.sum()),
                     "auc": roc_auc_score(g.good, x.fillna(x.min())) if 0 < g.good.sum() < len(g) else np.nan,
                     "spearman": spearmanr(g.dockq, x, nan_policy="omit")[0],
                     "top1_good": bool(o.good.iloc[0]), "top5_any_good": bool(o.good.head(5).any()),
                     "good_in_top10": int(o.good.head(10).sum()), "random_top10": 10 * g.good.mean(),
                     "top1_dockq": o.dockq.iloc[0], "best_dockq": g.dockq.max()})
pr = pd.DataFrame(rows); pr.to_csv(OUT / "per_run_metrics.csv", index=False)
summ = (pr.groupby("score").agg(median_auc=("auc", "median"), mean_auc=("auc", "mean"), median_rho=("spearman", "median"),
                                top1_success=("top1_good", "mean"), top5_success=("top5_any_good", "mean"),
                                mean_good_top10=("good_in_top10", "mean"), random_top10=("random_top10", "mean"),
                                runs=("auc", "size"))
        .sort_values("median_auc", ascending=False))
summ.insert(0, "label", summ.index.map(LABEL))
summ.round(3).to_csv(OUT / "summary_by_score.csv")
print("\n== per target x weight ranking (median over runs) ==\n", summ.round(3).to_string())

# ---- one global threshold across targets (pooled), per weight ----
pool = []
for w, g in d.groupby("weight"):
    for k in S:
        x = g[k].fillna(g[k].min())
        fpr, tpr, thr = roc_curve(g.good, x)
        j = np.argmax(tpr - fpr)
        sel = x >= thr[j]
        pool.append({"weight": w, "score": k, "pooled_auc": roc_auc_score(g.good, x), "youden_threshold": thr[j],
                     "sensitivity": tpr[j], "specificity": 1 - fpr[j], "precision": g.good[sel].mean(),
                     "base_rate": g.good.mean()})
pool = pd.DataFrame(pool); pool.round(3).to_csv(OUT / "pooled_by_score.csv", index=False)
print("\n== pooled over targets: AUC with one threshold for all targets ==\n",
      pool.pivot(index="score", columns="weight", values="pooled_auc").round(2).to_string())

# ---- does combining help? logistic regression, leave one target out, per weight ----
FEATS = ["plddt", "ptm", "i_ptm", "actifptm", "neg_pae", "neg_pae_contacts", "dockq_io", "fnat_io", "neg_irmsd_io"]
if HAVE_IFACE:
    FEATS += ["ipsae", "lis", "pdockq", "neg_pae_inter_min"]
lo = []
for w, g in d.groupby("weight"):
    g = g.dropna(subset=FEATS)
    for t in sorted(g.target.unique()):
        tr, te = g[g.target != t], g[g.target == t]
        if te.good.nunique() < 2:
            continue
        sc = StandardScaler().fit(tr[FEATS])
        m = LogisticRegression(max_iter=2000, C=0.5).fit(sc.transform(tr[FEATS]), tr.good)
        p = m.predict_proba(sc.transform(te[FEATS]))[:, 1]
        o = te.assign(p=p).sort_values("p", ascending=False)
        lo.append({"weight": w, "target": t, "auc": roc_auc_score(te.good, p), "top1_good": bool(o.good.iloc[0]),
                   "good_in_top10": int(o.good.head(10).sum()),
                   **{f"coef_{f}": c for f, c in zip(FEATS, m.coef_[0])}})
lo = pd.DataFrame(lo); lo.round(3).to_csv(OUT / "loto_logistic.csv", index=False)
print("\n== logistic regression on all scores, trained on the other targets ==")
print(lo.groupby("weight")[["auc", "top1_good", "good_in_top10"]].mean().round(2).to_string())
print("median AUC", round(lo.auc.median(), 3), " top-1 success", round(lo.top1_good.mean(), 3))
print("mean standardized coefficients:\n", lo.filter(like="coef_").mean().sort_values().round(2).to_string())

print("\nwrote", OUT)
