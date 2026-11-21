"""ipTM vs DockQ input->output over ALL conformers of one curated NMR sub-benchmark, one panel per weight set.

Usage: python plot_ens_subset.py <collected ens csv> <title> <out png> [expected target runs]
Small faint points: every conformer x weight; large points: per-dimer mean over conformers and weights of the set.
"""
import sys

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

d = pd.read_csv(sys.argv[1]); title, out = sys.argv[2], sys.argv[3]
expected = int(sys.argv[4]) if len(sys.argv) > 4 else None
d = d.rename(columns={"i_ptm": "iptm"})
for c in ("iptm", "dockq_io"):
    d[c] = pd.to_numeric(d[c], errors="coerce")
SETS = [("ptm", "AF2-ptm", "#52514e", "o"), ("multimer", "AF2.1 multimer", "#2a78d6", "s"),
        ("multimer_v2", "AF2.2 multimer", "#eb6834", "^"), ("multimer_v3", "AF2.3 multimer", "#1baf7a", "D")]
d["set"] = d["weight"].str.split("_", n=2).str[2].map({s[0]: s[1] for s in SETS})


def rho(x):
    return x[["dockq_io", "iptm"]].corr(method="spearman").iloc[0, 1] if len(x) > 2 else float("nan")


runs = d.groupby(["target", "weight"]).ngroups
fig, axes = plt.subplots(1, 4, figsize=(16, 5.6), dpi=200, sharex=True, sharey=True)
for ax, (_, name, col, mk) in zip(axes, SETS):
    s = d[d["set"] == name]
    m = s.groupby("target")[["iptm", "dockq_io"]].mean()
    ax.scatter(s["dockq_io"], s["iptm"], s=5, c=col, marker=mk, alpha=0.18, edgecolors="none", zorder=2)
    ax.scatter(m["dockq_io"], m["iptm"], s=46, c=col, marker=mk, edgecolors="white", linewidths=0.8, zorder=4)
    crowded = len(m) > 15  # label only dimers outside the main cluster when there are many
    for t, r in m.iterrows():
        if crowded and r["iptm"] >= 0.8 and r["dockq_io"] >= 0.85:
            continue
        ax.annotate(t.upper(), (r["dockq_io"], r["iptm"]), fontsize=6, xytext=(3, 3), textcoords="offset points", zorder=5)
    ax.set_title(f"{name} · {s.groupby(['target', 'weight']).ngroups} runs, {len(s)} conformer scores\n"
                 f"ρ all conformers {rho(s):.2f} · ρ dimer means {rho(m):.2f}\n"
                 f"mean ipTM {s['iptm'].mean():.2f} · DockQ {s['dockq_io'].mean():.2f}", loc="left", fontsize=8.5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal"); ax.set_xlabel("DockQ input→output")
    ax.grid(color="#e4e3df", linewidth=0.6, zorder=0)
    for side in ("top", "right"): ax.spines[side].set_visible(False)
axes[0].set_ylabel("ipTM")
done = f"{runs} of {expected} weight×dimer runs done" if expected else f"{runs} weight×dimer runs"
single = d.groupby(["target", "weight"]).size().max() == 1
what = "one model per dimer · small points: each weight, large: dimer means" if single else "all conformers · small points: conformers, large: dimer means"
fig.suptitle(f"{title} · {d['target'].nunique()} dimers · {done} · {what}",
             x=0.01, ha="left", fontsize=10.5)
fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(out); plt.close(fig)
print("wrote", out)
