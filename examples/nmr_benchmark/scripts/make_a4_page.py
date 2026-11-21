"""Lay out dimer thumbnails on one A4 page (PDF + PNG): homodimers first, then heterodimers.

Usage: python make_a4_page.py <png dir> <table.csv> <id col> <type col> <residues col> <title> <out stem> [flag col] [flag label] [flag2 col] [flag2 label]
<type col> values containing "homo" count as homodimers; <residues col> is "nA+nB".
Optional [flag col] (boolean column): flagged panels get a red frame and a red label; [flag label] explains it.
Optional [flag2 col]: a dashed purple frame and a dagger after the label. Pass "-" as [flag col] to skip the first flag.
"""
import sys
from pathlib import Path

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageChops

png_dir, table, idc, typec, resc, title, out = sys.argv[1:8]
flagc = sys.argv[8] if len(sys.argv) > 8 and sys.argv[8] != "-" else None
flag_label = sys.argv[9] if len(sys.argv) > 9 else "flagged"
flag2c = sys.argv[10] if len(sys.argv) > 10 else None
flag2_label = sys.argv[11] if len(sys.argv) > 11 else "flagged"
t = pd.read_csv(table)
if "excluded" in t:  # dimers excluded from the benchmark are left off the page
    t = t[t["excluded"].fillna("") == ""]
t["id"] = t[idc].str.lower()
t["homo"] = t[typec].str.contains("homo")
t = t.sort_values(["homo", "id"], ascending=[False, True]).reset_index(drop=True)
TEAL, AMBER, INK, MUTED, FLAG = "#0f9d8a", "#e5a83b", "#1f1e1b", "#77756e", "#c8322f"
t["flag"] = t[flagc].astype(str).str.lower().isin(["true", "1", "yes"]) if flagc else False
t["flag2"] = t[flag2c].astype(str).str.lower().isin(["true", "1", "yes"]) if flag2c else False
FLAG2 = "#7a3db8"


def trimmed(path):
    im = Image.open(path).convert("RGB")
    box = ImageChops.difference(im, Image.new("RGB", im.size, "white")).getbbox()
    im = im.crop(box) if box else im
    w, h = im.size; s = max(w, h)
    sq = Image.new("RGB", (s, s), "white"); sq.paste(im, ((s - w) // 2, (s - h) // 2))
    return np.asarray(sq)


n = len(t)
cols = max(1, int(np.ceil(np.sqrt(n * (0.94 * 8.27) / (0.90 * 11.69)))))  # roughly square cells on A4
rows = int(np.ceil(n / cols)); fs = min(8.5, 5.6 * 9 / cols)
fig = plt.figure(figsize=(8.27, 11.69), dpi=300)
top, bottom, left, right = (0.928 if (flagc or flag2c) else 0.935), 0.025, 0.03, 0.97
cw, rh = (right - left) / cols, (top - bottom) / rows
fig.text(left, 0.968, title, fontsize=11, color=INK, weight="bold", va="center")
fig.text(left, 0.950, f"{t['homo'].sum()} homodimers (first), {(~t['homo']).sum()} heterodimers · "
         "first chain teal, second chain amber · labels: PDB ID, residues per chain (homodimers in bold)", fontsize=6.5, color=MUTED, va="center")
notes = []
if flagc:
    notes.append((f"Red frame and label: {flag_label} ({int(t['flag'].sum())} dimers)", FLAG))
if flag2c:
    notes.append((f"Dashed purple frame and †: {flag2_label} ({int(t['flag2'].sum())} dimers)", FLAG2))
for j, (txt, colr) in enumerate(notes):
    fig.text(left + j * 0.47, 0.938, txt, fontsize=6.5, color=colr, va="center")
for i, r in t.iterrows():
    row, col = divmod(i, cols)
    x0, y0 = left + col * cw, top - (row + 1) * rh
    lab_h = 0.016 * fs / 5.6
    ax = fig.add_axes([x0 + 0.003, y0 + lab_h, cw - 0.006, rh - lab_h - 0.003])
    f = Path(png_dir) / f"{r['id']}.png"
    if f.is_file():
        ax.imshow(trimmed(f), interpolation="lanczos")
    ax.set_axis_off()
    if r["flag"]:
        fig.patches.append(plt.Rectangle((x0 + 0.002, y0 + 0.002), cw - 0.004, rh - 0.004, transform=fig.transFigure,
                                         fill=False, edgecolor=FLAG, linewidth=0.8))
    if r["flag2"]:
        fig.patches.append(plt.Rectangle((x0 + 0.0045, y0 + 0.0035), cw - 0.009, rh - 0.007, transform=fig.transFigure,
                                         fill=False, edgecolor=FLAG2, linewidth=0.8, linestyle=(0, (3, 2))))
    fig.text(x0 + cw / 2, y0 + lab_h * 0.55, f"{r['id'].upper()}  {r[resc]}" + (" †" if r["flag2"] else ""), fontsize=fs, ha="center", va="center",
             color=FLAG if r["flag"] else INK, weight="bold" if r["homo"] else "normal")
fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.png", dpi=200)
print("wrote", f"{out}.pdf", f"{out}.png", len(t), "panels")
