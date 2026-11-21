#!/usr/bin/env python3
"""Interface secondary structure of the prepared X-ray dimers, and the per-target summary.

Secondary structure: PyMOL `dss` (DSSP-like; H = helix, S = strand, L = other) on the dimer as written
(both chains loaded, so intermolecular sheets count). Interface residue: any heavy atom < 5.0 A from a heavy atom
of the partner chain (same definition as interface_residues_model1 in nmr99_conformers.csv).

Usage: ~/opt/anaconda3/envs/pymol/bin/python -I interface_ss.py <xray_dimers.csv> <xray99_inputs dir> <clean dir>
Writes <inputs>/xray99_interface_ss.csv and <inputs>/xray99_summary.csv.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from pymol import cmd

dimers_csv, inputs, clean_dir = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
rows = list(csv.DictReader(open(dimers_csv)))
report = json.load(open(inputs / "prep_report.json"))
ss_rows, summary = [], []
for r in rows:
    pid = r["pdb_id"].upper()
    pdb = inputs / pid.lower() / "models" / f"{pid.lower()}_model1.pdb"
    cmd.reinitialize()
    cmd.load(str(pdb), "m")
    cmd.dss("m")
    atoms = {"A": [], "B": []}
    ss = {}
    model = cmd.get_model("m")
    for a in model.atom:
        key = (a.chain, a.resi)
        atoms[a.chain].append((key, a.coord))
        if a.name == "CA":
            ss[key] = a.ss or "L"
    keys = {c: [k for k, _ in atoms[c]] for c in "AB"}
    xyz = {c: np.array([x for _, x in atoms[c]]) for c in "AB"}
    d = np.sqrt(((xyz["A"][:, None, :] - xyz["B"][None, :, :]) ** 2).sum(-1))
    close = d < 5.0
    iface = {"A": {keys["A"][i] for i in np.where(close.any(1))[0]},
             "B": {keys["B"][j] for j in np.where(close.any(0))[0]}}
    out = {"pdb_id": pid, "type": r["type"]}
    flags = []
    for c in "AB":
        res = sorted(iface[c], key=lambda k: (int("".join(ch for ch in k[1] if ch.isdigit() or ch == "-")), k[1]))
        n = len(res)
        h = sum(ss.get(k) == "H" for k in res)
        s = sum(ss.get(k) == "S" for k in res)
        allc = [v for k, v in ss.items() if k[0] == c]
        out.update({f"iface_res_{c}": n, f"frac_helix_{c}": round(h / n, 3) if n else "",
                    f"frac_strand_{c}": round(s / n, 3) if n else "", f"frac_helix_or_strand_{c}": round((h + s) / n, 3) if n else "",
                    f"chain_helix_{c}": round(sum(v == "H" for v in allc) / len(allc), 3),
                    f"chain_strand_{c}": round(sum(v == "S" for v in allc) / len(allc), 3)})
        if n and h + s == 0:
            flags.append(f"no_ss_at_interface_{c}")
        if n == 0:
            flags.append(f"no_interface_{c}")
    out["n_modelled_A"], out["n_modelled_B"] = report[pid]["modelled_residues"]
    out["flag"] = ";".join(flags)
    ss_rows.append(out)

    # pipeline strict clean result
    rec = [json.loads(l) for l in open(clean_dir / pid.lower() / "manifests" / "cleaning.jsonl")]
    clean_ok = len(rec) == 1 and rec[0]["status"] == "ok" and rec[0]["cleaned_residue_count"] == sum(report[pid]["modelled_residues"])
    sflags = list(flags)
    if not clean_ok:
        sflags.append("strict_clean_failed")
    if r["symmetry_mate"] == "yes":
        sflags.append("symmetry_mate")
    if any("symmetry copy" in n for n in report[pid]["notes"]):
        sflags.append("chain_from_symmetry_copy")
    if r["type"] == "homo" and not report[pid]["homodimer"]:
        sflags.append("homo_with_different_sequences")
    if report[pid]["nonstandard"]:
        sflags.append("nonstandard_replaced")
    if r["synthetic_or_nanobody"]:
        sflags.append("synthetic_or_nanobody")
    summary.append({"pdb_id": pid, "type": r["type"],
                    "residues_scored": f"{report[pid]['modelled_residues'][0]}+{report[pid]['modelled_residues'][1]}",
                    "fasta_lengths": f"{report[pid]['lengths'][0]}+{report[pid]['lengths'][1]}",
                    "interface_residues": out["iface_res_A"] + out["iface_res_B"],
                    "iface_frac_helix_or_strand_A": out["frac_helix_or_strand_A"],
                    "iface_frac_helix_or_strand_B": out["frac_helix_or_strand_B"],
                    "strict_clean": "pass" if clean_ok else "FAIL", "flags": ";".join(sflags)})

ss_csv = inputs / "xray99_interface_ss.csv"
with open(ss_csv, "w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=list(ss_rows[0]))
    w.writeheader(); w.writerows(ss_rows)
with open(inputs / "xray99_summary.csv", "w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=list(summary[0]))
    w.writeheader(); w.writerows(summary)
fl = [x for x in ss_rows if x["flag"]]
print(len(ss_rows), "targets;", len(fl), "flagged:")
for x in fl:
    print(" ", x["pdb_id"], x["type"], x["flag"], "iface", x["iface_res_A"], x["iface_res_B"],
          "len", x["n_modelled_A"], x["n_modelled_B"])
print("passed strict clean:", sum(s["strict_clean"] == "pass" for s in summary))
