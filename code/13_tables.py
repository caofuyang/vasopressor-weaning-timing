#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把三张结果表转成英文 markdown，供稿件直接使用"""
import os
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(BASE, "results")

MAP = {
    "age": "Age, years", "dose_t0": "Norepinephrine dose at eligibility, µg/kg/min",
    "map_t0": "MAP at eligibility, mmHg", "map_mean_6h": "Mean MAP, preceding 6 h, mmHg",
    "map_min_6h": "Minimum MAP, preceding 6 h, mmHg",
    "frac_map_lt65": "Fraction of preceding 6 h with MAP <65 mmHg",
    "hours_ne_before_t0": "Norepinephrine duration before eligibility, h",
    "hours_icu_before_t0": "ICU stay before eligibility, h",
    "lactate": "Lactate, mmol/L", "creatinine": "Creatinine, mg/dL",
    "bilirubin": "Total bilirubin, mg/dL", "platelet": "Platelet count, ×10⁹/L",
    "wbc": "White cell count, ×10⁹/L", "hemoglobin": "Haemoglobin, g/dL",
    "bun": "Blood urea nitrogen, mg/dL", "gcs_min": "Minimum GCS",
    "urine_6h": "Urine output, preceding 6 h, mL",
    "input_6h": "Fluid input, preceding 6 h, mL",
    "fluid_balance_6h": "Fluid balance, preceding 6 h, mL",
    "sofa": "SOFA score",
    "on_vent": "Invasive ventilation", "on_rrt": "Renal replacement therapy",
    "on_vasopressin": "Vasopressin co-administration",
}

def label(v, lab):
    if v in MAP:
        return MAP[v]
    if v.startswith("gender_"):
        return "Sex: " + {"F": "female", "M": "male"}.get(v.split("_", 1)[1], v)
    if v.startswith("adm_"):
        return "Admission type: " + v.split("_", 1)[1].capitalize()
    if v.startswith("unit_"):
        u = v.split("_", 1)[1]
        return "ICU type: " + ("Other" if u == "其他/未分类" else u)
    if v.endswith("_miss"):
        return "Missing: " + MAP.get(v[:-5], v[:-5]).split(",")[0].lower()
    return lab

def fmt_num(x):
    try:
        v = float(x)
    except Exception:
        return x
    return f"{v:.0f}" if abs(v - round(v)) < 1e-9 else f"{v:.2f}".rstrip("0").rstrip(".")

def fmt_cell(txt):
    t = str(txt)
    if " (" in t and t.endswith(")"):
        head, rest = t.split(" (", 1)
        body = rest[:-1]
        if "–" in body:
            a, b = body.split("–", 1)
            return f"{fmt_num(head)} ({fmt_num(a)}–{fmt_num(b)})"
    return t

t1 = pd.read_csv(os.path.join(R, "table1.csv"))
lines = ["**Table 1. Baseline characteristics by exposure group (observed data).**",
         "",
         f"| Characteristic | Weaning within 6 h (n={t1.columns[2].split('=')[1][:-1]}) | Weaning delayed (n={t1.columns[3].split('=')[1][:-1]}) | SMD |",
         "|---|---|---|---|"]
for _, r in t1.iterrows():
    if r["var"] == "n":
        continue
    smd = "" if pd.isna(r["SMD"]) else f"{r['SMD']:.3f}"
    a = fmt_cell(r[t1.columns[2]]); b = fmt_cell(r[t1.columns[3]])
    lines.append(f"| {label(r['var'], r['label'])} | {a} | {b} | {smd} |")
lines.append("")
lines.append("Values are median (IQR) or n (%). SMD, standardised mean difference.")
lines.append("")

t2 = pd.read_csv(os.path.join(R, "table_sensitivity.csv"))
NAME = {"early_wean_6h": "Primary: weaning within 6 h",
        "early_wean_2h": "Exposure window 2 h", "early_wean_1h": "Exposure window 1 h",
        "stop6h": "Complete discontinuation only",
        "truncate_2.5_97.5": "Weight truncation 2.5/97.5",
        "common_support_0.05_0.95": "Restricted to common support (PS 0.05–0.95)",
        "complete_case": "Complete-case analysis (including bilirubin)"}
lines += ["**Table 2. Sensitivity analyses for 28-day mortality.**", "",
          "| Analysis | HR (95% CI) | p |", "|---|---|---|"]
for _, r in t2.iterrows():
    lines.append(f"| {NAME.get(r['analysis'], r['analysis'])} | "
                 f"{r['HR']:.2f} ({r['lo']:.2f}–{r['hi']:.2f}) | "
                 f"{'<0.001' if r['p'] < 0.001 else f'{r[chr(112)]:.3f}'} |")
lines.append("")
lines.append("All models weighted by the propensity score; estimates pooled across 20 imputations.")
lines.append("")

t3 = pd.read_csv(os.path.join(R, "table3_subgroups.csv"))
lines += ["**Table 3. Subgroup analyses.**", "",
          "| Subgroup | n | Events | HR (95% CI) | p |", "|---|---|---|---|---|"]
for _, r in t3.iterrows():
    lines.append(f"| {r['subgroup']} | {int(r['n'])} | {int(r['events'])} | "
                 f"{r['HR']:.2f} ({r['lo']:.2f}–{r['hi']:.2f}) | "
                 f"{'<0.001' if r['p'] < 0.001 else f'{r[chr(112)]:.3f}'} |")
lines.append("")
lines.append("Propensity scores re-estimated within each subgroup; estimates pooled across 10 imputations.")

out = os.path.join(R, "tables_en.md")
open(out, "w", encoding="utf-8").write("\n".join(lines))
print("\n".join(lines[:8]))
print("...")
print("written:", out)
