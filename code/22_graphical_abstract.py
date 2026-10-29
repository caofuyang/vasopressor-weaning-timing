#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Graphical abstract for ACCPM submission.
Requirement: >=1328 px wide, longest side horizontal, readable at 13 cm.
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "results")
K = json.load(open(os.path.join(OUT, "key_numbers.json"), encoding="utf-8"))
M, E = K["main"], K["eicu"]

NAVY = "#1f3b57"
RED = "#a63d3d"
GREY = "#6b7280"
LIGHT = "#eef2f6"
WHITE = "#ffffff"

fig, ax = plt.subplots(figsize=(13.3, 4.6), dpi=170)   # 2261 x 782 px
ax.set_xlim(0, 100); ax.set_ylim(0, 34); ax.axis("off")

def box(x, y, w, h, fc, ec, lw=1.2, r=0.9):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle=f"round,pad=0.35,rounding_size={r}",
                 linewidth=lw, facecolor=fc, edgecolor=ec, zorder=2))

def arrow(x1, y1, x2, y2, color=GREY, lw=1.5):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                 mutation_scale=15, linewidth=lw, color=color, zorder=3,
                 shrinkA=0, shrinkB=0))

# ---------- title strip ----------
ax.text(1, 31.6, "Timing of the first vasopressor weaning attempt and 28-day mortality in sepsis",
        fontsize=14.5, fontweight="bold", color=NAVY, va="center")
ax.text(1, 28.9, "Target trial emulation in MIMIC-IV, with external validation in eICU-CRD",
        fontsize=9.5, color=GREY, va="center")

# ---------- panel 1: population & eligibility ----------
box(1, 12.6, 27, 13.2, LIGHT, NAVY)
ax.text(2.6, 23.6, "Population", fontsize=10.5, fontweight="bold", color=NAVY)
ax.text(2.6, 21.0, "Adults with sepsis receiving\nnorepinephrine", fontsize=9, color="#111")
ax.text(2.6, 17.4, "Eligibility (t₀), first hour with:", fontsize=8.6, color=GREY)
ax.text(2.6, 15.7, "NE ≤ 0.10 µg/kg/min  ·  MAP ≥ 65 mmHg  ·  no prior reduction",
        fontsize=7.8, color="#111")
ax.text(2.6, 14.0, f"Landmark = t₀ + 6 h      n = {M['n']:,}", fontsize=8.6, color="#111")

arrow(29.2, 19.2, 33.0, 19.2)

# ---------- panel 2: exposure split ----------
box(33.0, 20.8, 27, 7.4, WHITE, NAVY)
ax.text(34.6, 26.4, "Early weaning", fontsize=10, fontweight="bold", color=NAVY)
ax.text(34.6, 24.6, "≥ 25% dose reduction or", fontsize=8.4, color="#111")
ax.text(34.6, 23.2, "discontinuation within 6 h", fontsize=8.4, color="#111")
ax.text(34.6, 21.7, f"n = {M['early']:,}  ({M['pct_early']}%)", fontsize=9, fontweight="bold", color=NAVY)

box(33.0, 10.4, 27, 7.4, WHITE, RED)
ax.text(34.6, 16.0, "Delayed weaning", fontsize=10, fontweight="bold", color=RED)
ax.text(34.6, 14.2, "no reduction ≥ 25% and", fontsize=8.4, color="#111")
ax.text(34.6, 12.8, "no discontinuation within 6 h", fontsize=8.4, color="#111")
ax.text(34.6, 11.3, f"n = {M['delayed']:,}  ({100-M['pct_early']:.1f}%)", fontsize=9, fontweight="bold", color=RED)

arrow(60.2, 24.5, 64.0, 24.5, NAVY)
arrow(60.2, 14.1, 64.0, 14.1, RED)

# ---------- panel 3: outcome ----------
box(64.0, 10.4, 34.5, 17.8, "#f7f9fb", NAVY)
ax.text(65.8, 26.2, "28-day mortality", fontsize=10.5, fontweight="bold", color=NAVY)
ax.text(65.8, 23.2, "IPTW-weighted Cox model", fontsize=8.4, color=GREY)

ax.text(65.8, 20.3, f"HR {M['hr']:.2f}  (95% CI {M['lo']:.2f}–{M['hi']:.2f})",
        fontsize=15, fontweight="bold", color=NAVY)
ax.text(65.8, 17.6, f"28-day survival  {M['surv_early']}% vs {M['surv_delayed']}%", fontsize=8.6, color="#111")
ax.text(65.8, 15.7, f"E-value {M['evalue']:.2f}   ·   consistent across 6 sensitivity analyses", fontsize=8.2, color=GREY)
ax.text(65.8, 13.4, f"External validation (eICU-CRD, n = {E['n']}):  HR {E['hr']:.2f} ({E['lo']:.2f}–{E['hi']:.2f})",
        fontsize=8.2, color=GREY, style="italic")
ax.text(65.8, 11.6, "directionally consistent, not significant", fontsize=8.2, color=GREY, style="italic")

# ---------- conclusion strip ----------
box(1, 2.2, 97.5, 6.4, NAVY, NAVY)
ax.text(50.0, 5.4, "Attempting weaning within 6 hours of eligibility was associated with a 28% lower hazard of death",
        fontsize=10.6, color="white", ha="center", va="center", fontweight="bold")

plt.tight_layout(pad=0.2)
for ext in ["png", "pdf"]:
    plt.savefig(os.path.join(OUT, f"graphical_abstract.{ext}"),
                bbox_inches="tight", facecolor="white")
print("saved", os.path.join(OUT, "graphical_abstract.png"))
from PIL import Image
im = Image.open(os.path.join(OUT, "graphical_abstract.png"))
print("size:", im.size)
