"""HADDOCK-style cluster ranking of the CAPRI decoys with masked actifpTM in place of the HADDOCK score.

Clustering as in HADDOCK: decoys are neighbours when they share enough interface contacts, min(FCC(i,j),
FCC(j,i)) >= cutoff (default 0.6); greedy clustering takes the decoy with most unassigned neighbours as a centre
together with those neighbours, and repeats; clusters with fewer than 4 members are dropped and their decoys are not
ranked. Clusters are ranked by the mean actifpTM of their best 4 members; cluster size is not used. Top-1 = best
decoy of the best cluster; top-5 = best decoy of each of the 5 best clusters. Compared with actifpTM alone.
Writes results/clustering_per_run.csv, results/clustering_by_version.csv, results/clustering_one_ranking.csv.
Usage (from benchmarks/capri): python scripts/clustering.py [cutoff]
"""
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import HERE, load  # noqa: E402

CUTOFF = float(sys.argv[1]) if len(sys.argv) > 1 else 0.6
VERSION = {"ptm": "AF2-ptm", "multimer": "AF2.1", "multimer_v2": "AF2.2", "multimer_v3": "AF2.3"}
d = load("masked"); d["version"] = d.family.map(VERSION)
MEDIUM_PLUS = {"medium", "high"}


def clusters(target, cutoff=CUTOFF, min_size=4):
    """Cluster label per decoy id; -1 = not in a cluster of at least min_size."""
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


def picks(g, score="actifptm"):
    """Top-5 lists: HADDOCK-style (best decoy of each of the 5 best clusters) and by score alone."""
    g = g.assign(cluster=g.model_id.map(CL[g.target.iloc[0]]).values); gc = g[g.cluster >= 0]
    rank = gc.groupby("cluster")[score].apply(lambda x: x.nlargest(4).mean()).sort_values(ascending=False)
    clustered = [gc[gc.cluster == c].sort_values(score, ascending=False).iloc[0] for c in rank.index[:5]]
    alone = [r for _, r in g.sort_values(score, ascending=False).head(5).iterrows()]
    return {"clustered": clustered, "actifpTM alone": alone}


def row(p, **kw):
    return dict(kw, top1=bool(p[0].good), top5=any(x.good for x in p), top1_class=p[0].main_classification,
                top1_medium_or_high=p[0].main_classification in MEDIUM_PLUS, top1_dockq=p[0].dockq)


# per run: one target scored with one weight
runs = []
for (t, w), g in d.groupby(["target", "weight"]):
    nbad = int((~g.good).sum())
    for m, p in picks(g).items():
        runs.append(row(p, target=t, weight=w, version=g.version.iloc[0], method=m, random_top1=g.good.mean(),
                        random_top5=1 - math.comb(nbad, 5) / math.comb(len(g), 5)))
runs = pd.DataFrame(runs); runs.to_csv(HERE / "results/clustering_per_run.csv", index=False)
byv = (runs.groupby(["version", "method"]).agg(runs=("top1", "size"), top1=("top1", "mean"), top5=("top5", "mean"),
                                               top1_medium_or_high=("top1_medium_or_high", "mean"),
                                               mean_dockq_top1=("top1_dockq", "mean")).round(3))
allruns = runs.groupby("method").agg(runs=("top1", "size"), top1=("top1", "mean"), top5=("top5", "mean"),
                                     top1_medium_or_high=("top1_medium_or_high", "mean"),
                                     mean_dockq_top1=("top1_dockq", "mean")).round(3)
byv = pd.concat([byv, allruns.assign(version="all 17").set_index("version", append=True).swaplevel()])
byv.to_csv(HERE / "results/clustering_by_version.csv")

# one ranking per version: actifpTM averaged over the version's models, then ranked once per target
one = []
for name, versions in {**{v: [v] for v in VERSION.values()}, "all 17": list(VERSION.values())}.items():
    for t, g in d[d.version.isin(versions)].groupby("target"):
        m = g.groupby("model_id").agg(actifptm=("actifptm", "mean"), good=("good", "first"), dockq=("dockq", "first"),
                                      main_classification=("main_classification", "first")).reset_index().assign(target=t)
        for meth, p in picks(m).items():
            one.append(row(p, version=name, method=meth, target=t))
one = pd.DataFrame(one); one.to_csv(HERE / "results/clustering_one_ranking.csv", index=False)

print(f"FCC cutoff {CUTOFF}, clusters of >= 4; random picking: top-1 {runs.random_top1.mean():.2f}, "
      f"top-5 {runs.random_top5.mean():.2f}")
print("\n== per run ==\n", byv.unstack("method").to_string())
print("\n== one ranking per version (targets of 9) ==\n",
      one.groupby(["version", "method"]).agg(top1=("top1", "sum"), top5=("top5", "sum"),
                                             mean_dockq_top1=("top1_dockq", "mean")).round(2).unstack("method").to_string())
print("\n== per target, share of the 17 weights (top-1) ==\n",
      runs.pivot_table(index="target", columns="method", values="top1").round(2).to_string())
