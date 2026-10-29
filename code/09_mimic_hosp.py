#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MIMIC 侧改用「院内死亡」结局重跑，与 eICU 可比"""
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
d = pd.read_csv(os.path.join(BASE, "analytic_v2.csv"))
d["age"] = d["age"].clip(upper=90)
for c in ["lactate", "creatinine", "bilirubin", "urine_6h"]:
    d["log_" + c.replace("bilirubin", "bili").replace("urine_6h", "urine")] = np.log1p(d[c])
d["log_lact"] = np.log1p(d["lactate"]); d["log_creat"] = np.log1p(d["creatinine"])
d["log_urine"] = np.log1p(d["urine_6h"])

dt = pd.to_datetime(d["death_time"], errors="coerce")
dc = pd.to_datetime(d["dischtime"], errors="coerce")
lm = pd.to_datetime(d["t_landmark"], errors="coerce")
end = pd.concat([dt, dc], axis=1).min(axis=1)
d["event_hosp"] = (dt.notna() & (dt <= dc)).astype(int)
d["time_hosp"] = (end - lm).dt.total_seconds() / 86400.0
d = d[d.time_hosp > 0].copy()

CONT = ["age", "dose_t0", "map_t0", "map_mean_6h", "map_min_6h", "frac_map_lt65",
        "hours_ne_before_t0", "hours_icu_before_t0", "log_lact", "log_creat",
        "platelet", "wbc", "hemoglobin", "bun", "gcs_min", "log_urine", "input_6h"]
BIN = ["on_vent", "on_rrt", "on_vasopressin", "urine_missing", "input_missing"]
CAT = ["gender", "admission_type", "insurance", "careunit_t0"]
M = 20
raw = d[CONT + BIN].copy()
sets = []
for m in range(M):
    imp = IterativeImputer(estimator=BayesianRidge(), sample_posterior=True,
                           max_iter=10, random_state=m, skip_complete=True)
    df = pd.DataFrame(imp.fit_transform(raw), columns=raw.columns)
    for c in CAT:
        df[c] = d[c].values
    for c in ["early_wean_6h", "early_wean_2h", "time_hosp", "event_hosp"]:
        df[c] = d[c].values
    sets.append(df)

def fit(df, expo="early_wean_6h"):
    X = pd.concat([df[CONT + BIN], pd.get_dummies(df[CAT], drop_first=True, dtype=float)], axis=1)
    X = X.loc[:, X.std() > 0]
    X = pd.DataFrame(StandardScaler().fit_transform(X), columns=X.columns)
    y = df[expo].values.astype(int)
    ps = LogisticRegression(max_iter=2000).fit(X, y).predict_proba(X)[:, 1]
    p = y.mean(); w = np.where(y == 1, p / ps, (1 - p) / (1 - ps))
    w = np.clip(w, *np.quantile(w, (0.01, 0.99)))
    mdf = pd.DataFrame({"time": np.clip(df["time_hosp"].values, 1e-3, None),
                        "event": df["event_hosp"].values.astype(int), "expo": y, "w": w})
    c = CoxPHFitter(); c.fit(mdf, "time", "event", weights_col="w", robust=True)
    s = c.summary.loc["expo"]
    return dict(coef=s["coef"], se=s["se(coef)"])

def pool(fits):
    Q = np.array([f["coef"] for f in fits]); U = np.array([f["se"] ** 2 for f in fits])
    m = len(fits); Qb = Q.mean(); T = U.mean() + (1 + 1 / m) * Q.var(ddof=1); se = np.sqrt(T)
    return dict(hr=np.exp(Qb), lo=np.exp(Qb - 1.96 * se), hi=np.exp(Qb + 1.96 * se),
                p=2 * stats.norm.sf(abs(Qb / se)))

r = pool([fit(df) for df in sets])
print("MIMIC 院内死亡结局: n=%d, 事件=%d, HR %.3f (%.3f-%.3f) p=%.4f"
      % (len(d), int(d.event_hosp.sum()), r["hr"], r["lo"], r["hi"], r["p"]))
open(os.path.join(BASE, "results", "mimic_hosp.txt"), "w", encoding="utf-8").write(
    "HR %.3f (%.3f-%.3f) n=%d events=%d" % (r["hr"], r["lo"], r["hi"], len(d), int(d.event_hosp.sum())))
print("按暴露分组的粗死亡率：")
print(d.groupby("early_wean_6h").event_hosp.agg(["count", "sum", "mean"]).to_string())
