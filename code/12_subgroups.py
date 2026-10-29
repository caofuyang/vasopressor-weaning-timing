#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Table 3：预设亚组分析（多重填补 + 亚组内重新估计倾向评分）"""
import os, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression, BayesianRidge
from sklearn.preprocessing import StandardScaler
from lifelines import CoxPHFitter
from scipy import stats

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "results"); os.makedirs(OUT, exist_ok=True)
M = 10

d = pd.read_csv(os.path.join(BASE, "analytic_v2.csv"))
d["age"] = d["age"].clip(upper=90)
for c, n in [("lactate", "log_lact"), ("creatinine", "log_creat"), ("urine_6h", "log_urine")]:
    d[n] = np.log1p(d[c])

CONT = ["age", "dose_t0", "map_t0", "map_mean_6h", "map_min_6h", "frac_map_lt65",
        "hours_ne_before_t0", "hours_icu_before_t0", "log_lact", "log_creat",
        "platelet", "wbc", "hemoglobin", "bun", "gcs_min", "log_urine", "input_6h"]
BIN = ["on_vent", "on_rrt", "on_vasopressin", "urine_missing", "input_missing"]
CAT = ["gender", "admission_type", "insurance", "careunit_t0"]

sets = []
for m in range(M):
    imp = IterativeImputer(estimator=BayesianRidge(), sample_posterior=True,
                           max_iter=10, random_state=m, skip_complete=True)
    df = pd.DataFrame(imp.fit_transform(d[CONT + BIN]), columns=CONT + BIN)
    for c in CAT:
        df[c] = d[c].values
    for c in ["early_wean_6h", "time_days", "event_28d"]:
        df[c] = d[c].values
    sets.append(df)

def fit(df):
    X = pd.concat([df[CONT + BIN], pd.get_dummies(df[CAT], drop_first=True, dtype=float)], axis=1)
    X = X.loc[:, X.std() > 0]
    X = pd.DataFrame(StandardScaler().fit_transform(X), columns=X.columns)
    y = df["early_wean_6h"].values.astype(int)
    ps = LogisticRegression(max_iter=3000).fit(X, y).predict_proba(X)[:, 1]
    p = y.mean(); w = np.where(y == 1, p / ps, (1 - p) / (1 - ps))
    w = np.clip(w, *np.quantile(w, (0.01, 0.99)))
    mdf = pd.DataFrame({"time": df["time_days"].values, "event": df["event_28d"].values,
                        "expo": y, "w": w})
    c = CoxPHFitter(); c.fit(mdf, "time", "event", weights_col="w", robust=True)
    s = c.summary.loc["expo"]
    return s["coef"], s["se(coef)"], int(mdf.event.sum())

def pool(fits):
    Q = np.array([f[0] for f in fits]); U = np.array([f[1] ** 2 for f in fits])
    m = len(fits); Qb = Q.mean(); T = U.mean() + (1 + 1 / m) * Q.var(ddof=1); se = np.sqrt(T)
    return np.exp(Qb), np.exp(Qb - 1.96 * se), np.exp(Qb + 1.96 * se), 2 * stats.norm.sf(abs(Qb / se))

subs = {
    "Mechanical ventilation — yes": d.on_vent == 1,
    "Mechanical ventilation — no": d.on_vent == 0,
    "Renal replacement — yes": d.on_rrt == 1,
    "Renal replacement — no": d.on_rrt == 0,
    "Vasopressin — yes": d.on_vasopressin == 1,
    "Vasopressin — no": d.on_vasopressin == 0,
    "Lactate >4 mmol/L": d.lactate.fillna(d.lactate.median()) > 4,
    "Lactate ≤4 mmol/L": d.lactate.fillna(d.lactate.median()) <= 4,
    "Prior NE >24 h": d.hours_ne_before_t0 > 24,
    "Prior NE ≤24 h": d.hours_ne_before_t0 <= 24,
}

rows = []
for lab, mask in subs.items():
    idx = np.where(mask.values)[0]
    fits = []
    for df in sets:
        sub = df.iloc[idx].reset_index(drop=True)
        try:
            fits.append(fit(sub))
        except Exception:
            continue
    if len(fits) < 3:
        continue
    hr, lo, hi, p = pool(fits)
    rows.append({"subgroup": lab, "n": len(idx), "events": fits[0][2],
                 "HR": hr, "lo": lo, "hi": hi, "p": p})

t3 = pd.DataFrame(rows)
t3.to_csv(os.path.join(OUT, "table3_subgroups.csv"), index=False, float_format="%.3f")
print(t3.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
