"""Composite scores, gates and consensus over the AF2 weights (full resolution, all 153 runs).

Per-target z-score sums (Gaudreault et al. 2023 Sci Rep), Bennett et al. 2023 style gates (confident AND kept),
AF2Rank-style confidence x input->output agreement, and rank averaging over weights (all, and within each family).
Writes results/strategies_per_run.csv, results/strategies_consensus.csv, results/strategies_by_family.csv.
Usage (from examples/capri_benchmark): python scripts/strategies.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import HERE, load  # noqa: E402

d = load("fullres")
d["actif_x_dockq"] = d.actifptm * d.dockq_io
HAVE_IFACE = "ipsae" in d
z = lambda s: d.groupby(["target", "weight"])[s].transform(lambda x: (x - x.mean()) / (x.std() + 1e-9))
for s in ["neg_pae_contacts", "actifptm", "dockq_io", "plddt", "ptm", "fnat_io"] + (["ipsae", "lis"] if HAVE_IFACE else []):
    d["z_" + s] = z(s)
C = {
    "PAE@contacts": d.neg_pae_contacts,
    "actifpTM×DockQ_io": d.actif_x_dockq,
    "Gaudreault z(pLDDT)+z(pTM)": d.z_plddt + d.z_ptm,
    "z(PAEc)+z(actifpTM)": d.z_neg_pae_contacts + d.z_actifptm,
    "z(PAEc)+z(actifpTM)+z(DockQ_io)": d.z_neg_pae_contacts + d.z_actifptm + d.z_dockq_io,
    "z(PAEc)+z(DockQ_io)": d.z_neg_pae_contacts + d.z_dockq_io,
    "z(PAEc)+z(actifpTM)+z(fnat_io)": d.z_neg_pae_contacts + d.z_actifptm + d.z_fnat_io,
    "gate DockQ_io≥0.5 → PAEc": np.where(d.dockq_io >= 0.5, d.neg_pae_contacts, d.neg_pae_contacts - 100),
    "gate iRMSD_io<2 → PAEc": np.where(d.irmsd_io < 2, d.neg_pae_contacts, d.neg_pae_contacts - 100),
}
if HAVE_IFACE:
    C.update({"ipSAE": d.ipsae, "LIS": d.lis, "pDockQ": d.pdockq,
              "z(ipSAE)+z(actifpTM)+z(fnat_io)": d.z_ipsae + d.z_actifptm + d.z_fnat_io,
              "z(PAEc)+z(ipSAE)+z(fnat_io)": d.z_neg_pae_contacts + d.z_ipsae + d.z_fnat_io})
for k, v in C.items():
    d[k] = v


def evaluate(frame, cols, by):
    out = []
    for k in cols:
        r = []
        for _, g in frame.groupby(by):
            if not (0 < g.good.sum() < len(g)):
                continue
            o = g.sort_values(k, ascending=False)
            r.append((roc_auc_score(g.good, g[k]), average_precision_score(g.good, g[k]) / g.good.mean(),
                      o.good.iloc[0], o.good.head(5).any(), o.good.head(10).sum(), 10 * g.good.mean()))
        r = np.array(r, dtype=float)
        out.append((k, np.median(r[:, 0]), np.median(r[:, 1]), r[:, 2].mean(), r[:, 3].mean(), r[:, 4].mean(), r[:, 5].mean(), len(r)))
    return pd.DataFrame(out, columns=["score", "median AUC", "PR-AUC/base", "top1", "top5", "good in top10",
                                      "random top10", "n"]).set_index("score").round(3)


per_run = evaluate(d, list(C), ["target", "weight"]).sort_values("top1", ascending=False)
for k in C:
    d["r_" + k] = d.groupby(["target", "weight"])[k].rank(pct=True)
cons = d.groupby(["target", "model_id"]).agg(good=("good", "first"), **{k: ("r_" + k, "mean") for k in C}).reset_index()
consensus = evaluate(cons, list(C), "target").sort_values("top1", ascending=False)
per_run.to_csv(HERE / "results/strategies_per_run.csv"); consensus.to_csv(HERE / "results/strategies_consensus.csv")
print("== per target x weight ==\n", per_run.to_string())
print("\n== consensus over the weights (one ranking per target) ==\n", consensus.to_string())

# by weight family: per-run results, and consensus within each family
KEY = ["PAE@contacts", "actifpTM×DockQ_io", "z(PAEc)+z(actifpTM)+z(fnat_io)", "Gaudreault z(pLDDT)+z(pTM)"]
fam_rows = []
for f, g in d.groupby("family"):
    pr = evaluate(g, KEY, ["target", "weight"]).assign(family=f, level="per run")
    cg = g.groupby(["target", "model_id"]).agg(good=("good", "first"), **{k: ("r_" + k, "mean") for k in KEY}).reset_index()
    cs = evaluate(cg, KEY, "target").assign(family=f, level="family consensus")
    fam_rows += [pr, cs]
fam = pd.concat(fam_rows).reset_index().set_index(["family", "level", "score"])
fam.to_csv(HERE / "results/strategies_by_family.csv")
print("\n== by weight family ==\n", fam[["median AUC", "top1", "top5", "good in top10", "random top10", "n"]].to_string())
