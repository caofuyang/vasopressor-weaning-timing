#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成 Table 1（按暴露分组的基线特征）与投稿用 markdown 表"""
import os
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "results")
os.makedirs(OUT, exist_ok=True)

d = pd.read_csv(os.path.join(BASE, "analytic_v2.csv"))
d["age"] = d["age"].clip(upper=90)
d["fluid_balance_6h"] = d["input_6h"] - d["urine_6h"]

CONT = {
    "age": "年龄（岁）", "dose_t0": "入组时 NE 剂量（µg/kg/min）",
    "map_t0": "入组时 MAP（mmHg）", "map_mean_6h": "前 6h 平均 MAP（mmHg）",
    "map_min_6h": "前 6h 最低 MAP（mmHg）", "frac_map_lt65": "前 6h MAP<65 时间占比",
    "hours_ne_before_t0": "NE 已用时长（h）", "hours_icu_before_t0": "ICU 已住时长（h）",
    "lactate": "乳酸（mmol/L）", "creatinine": "肌酐（mg/dL）",
    "bilirubin": "总胆红素（mg/dL）", "platelet": "血小板（×10⁹/L）",
    "wbc": "白细胞（×10⁹/L）", "hemoglobin": "血红蛋白（g/dL）",
    "bun": "BUN（mg/dL）", "gcs_min": "GCS 最低分",
    "urine_6h": "6h 尿量（mL）", "input_6h": "6h 入量（mL）",
    "fluid_balance_6h": "6h 液体平衡（mL）",
    "sofa": "SOFA 评分",
}
BIN = {
    "on_vent": "机械通气", "on_rrt": "肾脏替代治疗",
    "on_vasopressin": "合用血管加压素",
}

# 稀有 ICU 类型合并为 Other
_keep = d["careunit_t0"].value_counts()
_keep = _keep[_keep >= 50].index
d["careunit_t0"] = np.where(d["careunit_t0"].isin(_keep), d["careunit_t0"], "其他/未分类")

def smd(a, b):
    a, b = a.dropna(), b.dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    return (a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)

g1 = d[d.early_wean_6h == 1]
g0 = d[d.early_wean_6h == 0]
rows = []

rows.append(("n", "", f"{len(g1)}", f"{len(g0)}", np.nan))
for c, lab in CONT.items():
    med1, q1a, q1b = g1[c].median(), g1[c].quantile(.25), g1[c].quantile(.75)
    med0, q0a, q0b = g0[c].median(), g0[c].quantile(.25), g0[c].quantile(.75)
    rows.append((c, lab,
                 f"{med1:.2f} ({q1a:.2f}–{q1b:.2f})",
                 f"{med0:.2f} ({q0a:.2f}–{q0b:.2f})",
                 smd(g1[c], g0[c])))
for c, lab in BIN.items():
    rows.append((c, lab, f"{int(g1[c].sum())} ({100*g1[c].mean():.1f}%)",
                 f"{int(g0[c].sum())} ({100*g0[c].mean():.1f}%)", smd(g1[c], g0[c])))

# 分类变量
def cat_row(col, lab, levels=None):
    tabs = []
    for lv in (levels or sorted(d[col].dropna().unique())):
        p1 = 100 * (g1[col] == lv).mean(); p0 = 100 * (g0[col] == lv).mean()
        tabs.append((lv, f"{int((g1[col]==lv).sum())} ({p1:.1f}%)",
                     f"{int((g0[col]==lv).sum())} ({p0:.1f}%)",
                     smd((g1[col] == lv).astype(float), (g0[col] == lv).astype(float))))
    return tabs

for lv, a, b, s in cat_row("gender", "性别"):
    rows.append(("gender_" + str(lv), f"  性别：{lv}", a, b, s))
for lv, a, b, s in cat_row("admission_type", "入院类型"):
    rows.append(("adm_" + str(lv), f"  入院类型：{lv}", a, b, s))
for lv, a, b, s in cat_row("careunit_t0", "ICU 类型"):
    rows.append(("unit_" + str(lv), f"  ICU：{lv}", a, b, s))

# 缺失
for c, lab in [("lactate", "乳酸"), ("bilirubin", "胆红素"), ("urine_6h", "尿量")]:
    m1 = 100 * g1[c].isna().mean(); m0 = 100 * g0[c].isna().mean()
    rows.append((c + "_miss", f"  缺失率：{lab}", f"{m1:.1f}%", f"{m0:.1f}%", np.nan))

C1, C0 = f"early(N={len(g1)})", f"late(N={len(g0)})"
t1 = pd.DataFrame(rows, columns=["var", "label", C1, C0, "SMD"])
t1.to_csv(os.path.join(OUT, "table1.csv"), index=False, float_format="%.3f")

# markdown
md = [f"| 变量 | 早撤组 (n={len(g1)}) | 延迟组 (n={len(g0)}) | SMD |", "|---|---|---|---|"]
for _, r in t1.iterrows():
    if r["var"] == "n":
        continue
    s = "" if pd.isna(r["SMD"]) else f"{r['SMD']:.3f}"
    md.append(f"| {r['label']} | {r[C1]} | {r[C0]} | {s} |")
open(os.path.join(OUT, "table1.md"), "w", encoding="utf-8").write("\n".join(md))
print("\n".join(md[:14]))
print("\n最大 |SMD| = %.3f" % t1.SMD.abs().max())
print("|SMD|>0.1 的变量：", int((t1.SMD.abs() > 0.1).sum()))
