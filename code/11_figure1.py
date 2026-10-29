#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 1：两库队列流程图"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "results")

MAIN = "#1f3b52"
EXCL = "#8a8a8a"
BOX_MAIN = "#eef2f6"
BOX_EXCL = "#f4f4f4"

def box(ax, x, y, w, h, text, fc, ec, fs=8.2, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.012",
                                linewidth=0.9, facecolor=fc, edgecolor=ec))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color="#111111", fontweight="bold" if bold else "normal", linespacing=1.4)

def arrow(ax, x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=9,
                                 linewidth=0.9, color="#444444", shrinkA=0, shrinkB=0))

def column(ax, x0, title, rows):
    w_main, w_excl = 0.205, 0.170
    h, gap = 0.058, 0.049
    y = 0.80
    ax.text(x0 + w_main / 2, 0.905, title, ha="center", va="center",
            fontsize=10, fontweight="bold", color=MAIN)
    last_main_y = None
    for text, is_excl in rows:
        if not is_excl:
            box(ax, x0, y, w_main, h, text, BOX_MAIN, MAIN)
            if last_main_y is not None:
                arrow(ax, x0 + w_main / 2, last_main_y, x0 + w_main / 2, y + h)
            last_main_y = y
            y -= (h + gap)
        else:
            ex = x0 + w_main + 0.020
            box(ax, ex, last_main_y + 0.006, w_excl, h * 0.78, text, BOX_EXCL, EXCL, fs=6.9)
            arrow(ax, x0 + w_main, last_main_y + h * 0.52, ex, last_main_y + h * 0.45)

fig, ax = plt.subplots(figsize=(9.6, 5.0), dpi=300)
ax.set_xlim(0, 1); ax.set_ylim(0.145, 0.95); ax.axis("off")

import json as _json, os as _os
_k = _json.load(open(_os.path.join(OUT, "key_numbers.json"), encoding="utf-8"))
F = _k["flow"]
M, E = _k["main"], _k["eicu"]

column(ax, 0.040, "MIMIC-IV (single centre)", [
    (f"{F['mimic_ne']:,} ICU stays receiving\nnorepinephrine", False),
    (f"{F['mimic_excl_no_hour']:,} excluded\n(no hour meeting\neligibility criteria)", True),
    (f"{F['mimic_eligible']:,} stays with ≥1 hour\nmeeting criteria", False),
    (f"{F['mimic_excl_no_sepsis']:,} excluded\n(no sepsis diagnosis code)", True),
    (f"{F['mimic_sepsis']:,} stays with sepsis", False),
    (f"{F['mimic_excl_filters']:,} excluded\n(not first ICU stay;\n<6 h prior norepinephrine)", True),
    (f"{F['mimic_after_filters']:,} analysed cohort\nat eligibility", False),
    (f"{F['mimic_excl_landmark']:,} excluded\n(before landmark)", True),
    (f"{F['mimic_analytic']:,} landmark cohort\n{M['early']:,} early / {M['delayed']:,} delayed", False),
])

column(ax, 0.530, "eICU-CRD (208 hospitals)", [
    (f"{F['eicu_ne']:,} ICU stays receiving\nnorepinephrine", False),
    (f"{F['eicu_ne'] - F['eicu_convertible']:,} excluded\n(volume-only records,\ndose not convertible)", True),
    (f"{F['eicu_convertible']:,} stays with a\nconvertible dose record", False),
    (f"{F['eicu_excl_no_elig']:,} excluded\n(never met eligibility\ncriteria)", True),
    (f"{F['eicu_eligible']:,} stays meeting\ncriteria at ≥1 hour", False),
    (f"{F['eicu_excl_no_sepsis']:,} excluded\n(no sepsis code;\nnot first stay)", True),
    (f"{F['eicu_sepsis']:,} stays with sepsis", False),
    (f"{F['eicu_excl_ne6h']:,} excluded\n(<6 h prior\nnorepinephrine)", True),
    (f"{F['eicu_after_ne6h']:,} analysed cohort\nat eligibility", False),
    (f"{F['eicu_excl_landmark']:,} excluded\n(before landmark)", True),
    (f"{E['n']:,} landmark cohort\n{E['early']:,} early / {E['n']-E['early']:,} delayed", False),
])



plt.tight_layout()
for ext in ["png", "pdf"]:
    plt.savefig(os.path.join(OUT, f"figure1_flow.{ext}"), bbox_inches="tight", facecolor="white")
print("saved", os.path.join(OUT, "figure1_flow.png"))
