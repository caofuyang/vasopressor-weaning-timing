#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 2：IPTW 加权 Kaplan-Meier 生存曲线（MIMIC 主队列）"""
import os, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression, BayesianRidge
from sklearn.preprocessing import StandardScaler
from lifelines import KaplanMeierFitter

warnings.filterwarnings("ignore")
BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "results"); os.makedirs(OUT, exist_ok=True)

d = pd.read_csv(os.path.join(BASE, "analytic_v2.csv"))
d["age"] = d["age"].clip(upper=90)
for c, n in [("lactate", "log_lact"), ("creatinine", "log_creat"), ("urine_6h", "log_urine")]:
    d[n] = np.log1p(d[c])

CONT = ["age", "dose_t0", "map_t0", "map_mean_6h", "map_min_6h", "frac_map_lt65",
        "hours_ne_before_t0", "hours_icu_before_t0", "log_lact", "log_creat",
        "platelet", "wbc", "hemoglobin", "bun", "gcs_min", "log_urine", "input_6h"]
BIN = ["on_vent", "on_rrt", "on_vasopressin", "urine_missing", "input_missing"]
CAT = ["gender", "admission_type", "insurance", "careunit_t0"]

imp = IterativeImputer(estimator=BayesianRidge(), sample_posterior=True,
                       max_iter=10, random_state=1, skip_complete=True)
df = pd.DataFrame(imp.fit_transform(d[CONT + BIN]), columns=CONT + BIN)
for c in CAT:
    df[c] = d[c]
X = pd.concat([df[CONT + BIN], pd.get_dummies(df[CAT], drop_first=True, dtype=float)], axis=1)
X = X.loc[:, X.std() > 0]
X = pd.DataFrame(StandardScaler().fit_transform(X), columns=X.columns)
y = d["early_wean_6h"].values
ps = LogisticRegression(max_iter=3000).fit(X, y).predict_proba(X)[:, 1]
p = y.mean()
w = np.clip(np.where(y == 1, p / ps, (1 - p) / (1 - ps)), *np.quantile(
    np.where(y == 1, p / ps, (1 - p) / (1 - ps)), (0.01, 0.99)))

fig, ax = plt.subplots(figsize=(6.2, 4.4), dpi=300)
n1, n0 = int((y == 1).sum()), int((y == 0).sum())
from lifelines import CoxPHFitter
_m = pd.DataFrame({"time": d["time_days"].values, "event": d["event_28d"].values,
                   "expo": y, "w": w})
_c = CoxPHFitter(); _c.fit(_m, "time", "event", weights_col="w", robust=True)
_s = _c.summary.loc["expo"]
hr, lo, hi = _s["exp(coef)"], _s["exp(coef) lower 95%"], _s["exp(coef) upper 95%"]
# 单一数据源
import json as _json, os as _os
_kp = _os.path.join(OUT, "key_numbers.json")
if _os.path.exists(_kp):
    _k = _json.load(open(_kp, encoding="utf-8"))["main"]
    hr, lo, hi = _k["hr"], _k["lo"], _k["hi"]
    n1, n0 = _k["early"], _k["delayed"]
for g, lab, col in [(1, f"Weaning attempted within 6 h (n={n1:,})", "#1f4e79"),
                    (0, f"Weaning delayed (n={n0:,})", "#b02418")]:
    m = y == g
    km = KaplanMeierFitter()
    km.fit(d.loc[m, "time_days"], d.loc[m, "event_28d"], weights=w[m])
    km.plot_survival_function(ax=ax, ci_show=True, color=col, linewidth=2, label=lab)
    ax.fill_between(km.survival_function_.index,
                    km.confidence_interval_.iloc[:, 0],
                    km.confidence_interval_.iloc[:, 1], color=col, alpha=0.12, linewidth=0)

ax.set_xlim(0, 28); ax.set_ylim(0.55, 1.0)
ax.set_xlabel("Days from landmark (eligibility + 6 h)", fontsize=10)
ax.set_ylabel("Survival probability", fontsize=10)
ax.set_xticks(range(0, 29, 7))
ax.grid(alpha=0.25, linewidth=0.6)
ax.legend(frameon=False, fontsize=8.5, loc="lower left")
ax.text(0.98, 0.96, f"IPTW-weighted\nHR {hr:.2f} (95% CI {lo:.2f}–{hi:.2f})",
        transform=ax.transAxes, fontsize=8.5, color="#333333", ha="right", va="top")
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
for ext in ["png", "pdf"]:
    plt.savefig(os.path.join(OUT, f"figure2_km.{ext}"), bbox_inches="tight")
print("saved:", os.path.join(OUT, "figure2_km.png"))
print("加权 28 天生存率：早撤 %.3f  延迟 %.3f" % (
    KaplanMeierFitter().fit(d.loc[y == 1, "time_days"], d.loc[y == 1, "event_28d"],
                            weights=w[y == 1]).predict(28),
    KaplanMeierFitter().fit(d.loc[y == 0, "time_days"], d.loc[y == 0, "event_28d"],
                            weights=w[y == 0]).predict(28)))
