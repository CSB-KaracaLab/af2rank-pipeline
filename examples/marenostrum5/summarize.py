"""Collect AF2Rank scores from one folder per weight; write a CSV and an ipTM vs DockQ plot."""
import sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

out = Path(sys.argv[1])
RELEASES = [("AF2-ptm", "_ptm", "#52514e", "o"), ("AF2.1 multimer", "_multimer", "#2a78d6", "s"),
            ("AF2.2 multimer", "_multimer_v2", "#eb6834", "^"), ("AF2.3 multimer", "_multimer_v3", "#1baf7a", "D")]
COLUMNS = {"model_id": "Model", "plddt": "pLDDT", "ptm": "pTM", "i_ptm": "ipTM", "actifptm": "actifpTM",
           "dockq_io": "DockQ input→output", "pae_contacts_mean": "Contact PAE (Å)"}

def release(name):
    for label, suffix, colour, marker in reversed(RELEASES):
        if name.endswith(suffix):
            return label, colour, marker

rows = []
for weight_dir in sorted(p for p in out.iterdir() if (p / "manifests/final_models.csv").is_file()):
    t = pd.read_csv(weight_dir / "manifests/final_models.csv", dtype={"model_id": str})
    label, colour, marker = release(weight_dir.name)
    t = t[[c for c in COLUMNS if c in t]].rename(columns=COLUMNS)
    rows.append(t.assign(Weight=f"{label.replace(' multimer', '')} {weight_dir.name.split('_')[1]}",
                         Release=label, colour=colour, marker=marker))
d = pd.concat(rows, ignore_index=True)
table = d.drop(columns=["colour", "marker", "Release"]).sort_values(["Model", "ipTM"], ascending=[True, False])
table.round(3).to_csv(out / "scores_all_weights.csv", index=False)
print(table.round(3).to_string(index=False))

models = sorted(d["Model"].unique())
fig, axes = plt.subplots(1, len(models), figsize=(5 * len(models), 4.8), dpi=200, squeeze=False)
fig.patch.set_facecolor("#fcfcfb")
for ax, model in zip(axes[0], models):
    ax.set_facecolor("#fcfcfb")
    sub = d[d["Model"] == model]
    for label, _, colour, marker in RELEASES:
        s = sub[sub["Release"] == label]
        if not s.empty:
            ax.scatter(s["DockQ input→output"], s["ipTM"], s=60, c=colour, marker=marker,
                       edgecolors="#fcfcfb", linewidths=1.2, label=label, zorder=3)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_title(model, loc="left", fontsize=10, color="#52514e")
    ax.set_xlabel("DockQ input→output"); ax.set_ylabel("ipTM")
    ax.grid(color="#e4e3df", linewidth=0.6, zorder=0)
    for side in ("top", "right"): ax.spines[side].set_visible(False)
axes[0][0].legend(loc="upper left", frameon=False, fontsize=8)
fig.suptitle(f"{out.name.rsplit('-', 1)[0]} · ipTM vs DockQ input→output · {d['Weight'].nunique()} AF2 weights",
             x=0.01, ha="left", fontsize=12)
fig.tight_layout(); fig.savefig(out / "iptm_vs_dockq.png", facecolor=fig.get_facecolor())
