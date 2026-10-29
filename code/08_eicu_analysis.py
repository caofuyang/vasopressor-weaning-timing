#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eICU 外部验证：与 MIMIC 同口径的 IPTW + 加权 Cox（结局：院内死亡）"""
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

d0 = pd.read_csv(os.path.join(BASE, "eicu_analytic.csv"))
if "event_hosp" in d0.columns:
    d0["event_28d"] = d0["event_hosp"]
d0["age"] = d0["age"].clip(upper=90)
d0["log_lact"] = np.log1p(d0["lactate"])
d0["log_creat"] = np.log1p(d0["creatinine"])
d0["log_urine"] = np.log1p(d0["urine_6h"])
d0["log_bili"] = np.log1p(d0["bilirubin"])

CONT = ["age", "weight", "dose_t0", "map_t0", "map_mean_6h", "map_min_6h",
        "frac_map_lt65", "hours_ne_before_t0", "log_lact", "log_creat",
        "wbc", "hemoglobin", "platelet", "bun", "log_urine"]
BIN = ["on_vent", "on_rrt", "on_vasopressin"]
CAT = ["gender", "unittype"]
M = 20

raw = d0[CONT + BIN].copy()
sets = []
for m in range(M):
    imp = IterativeImputer(estimator=BayesianRidge(), sample_posterior=True,
                           max_iter=10, random_state=m, skip_complete=True)
    df = pd.DataFrame(imp.fit_transform(raw), columns=raw.columns)
    for c in CAT:
        df[c] = d0[c].values
    for c in ["early_wean_6h", "early_wean_2h", "time_days", "event_28d", "hospitalid"]:
        df[c] = d0[c].values
    sets.append(df)

def fit(df, expo="early_wean_6h", weighted=True):
    X = pd.concat([df[CONT + BIN], pd.get_dummies(df[CAT], drop_first=True, dtype=float)], axis=1)
    X = X.loc[:, X.std() > 0]
    X = pd.DataFrame(StandardScaler().fit_transform(X), columns=X.columns)
    y = df[expo].values.astype(int)
    ps = LogisticRegression(max_iter=2000).fit(X, y).predict_proba(X)[:, 1]
    p = y.mean()
    w = np.where(y == 1, p / ps, (1 - p) / (1 - ps))
    w = np.clip(w, *np.quantile(w, (0.01, 0.99)))
    mdf = pd.DataFrame({"time": np.clip(df["time_days"].values, 1e-3, None),
                        "event": df["event_28d"].values.astype(int),
                        "expo": y, "w": w if weighted else np.ones(len(w)),
                        "hosp": df["hospitalid"].values})
    c = CoxPHFitter()
    c.fit(mdf, "time", "event", weights_col="w", cluster_col="hosp", robust=True)
    s = c.summary.loc["expo"]
    return dict(coef=s["coef"], se=s["se(coef)"], hr=s["exp(coef)"],
                lo=s["exp(coef) lower 95%"], hi=s["exp(coef) upper 95%"], p=s["p"],
                n=len(mdf), ev=int(mdf.event.sum()), wmax=float(w.max()))

def pool(fits):
    Q = np.array([f["coef"] for f in fits]); U = np.array([f["se"] ** 2 for f in fits])
    m = len(fits); Qb = Q.mean(); Ub = U.mean(); B = Q.var(ddof=1)
    T = Ub + (1 + 1 / m) * B; se = np.sqrt(T)
    lam = (1 + 1 / m) * B / T
    df_old = (m - 1) / lam ** 2 if lam > 0 else np.inf
    dfc = fits[0]["n"] - 2
    dfn = dfc * df_old / (dfc + df_old) if np.isfinite(df_old) else dfc
    t = stats.t.ppf(0.975, dfn)
    return dict(hr=np.exp(Qb), lo=np.exp(Qb - t * se), hi=np.exp(Qb + t * se),
                p=2 * stats.t.sf(abs(Qb / se), dfn))

main = pool([fit(df) for df in sets])
unadj = pool([fit(df, weighted=False) for df in sets])
print("eICU 主分析（IPTW, 院内死亡）: HR %.3f (%.3f-%.3f) p=%.4f | n=%d 事件=%d"
      % (main["hr"], main["lo"], main["hi"], main["p"], fit(sets[0])["n"], fit(sets[0])["ev"]))
print("未加权: HR %.3f (%.3f-%.3f) p=%.4f" % (unadj["hr"], unadj["lo"], unadj["hi"], unadj["p"]))
_f0 = fit(sets[0])
open(os.path.join(OUT, "eicu_summary.txt"), "w", encoding="utf-8").write(
    "HR %.3f (%.3f-%.3f) p=%.4f | n=%d 事件=%d" % (main["hr"], main["lo"], main["hi"], main["p"], _f0["n"], _f0["ev"]))

rows = []
for expo in ["early_wean_6h", "early_wean_2h"]:
    r = pool([fit(df, expo=expo) for df in sets])
    rows.append({"analysis": expo, "HR": r["hr"], "lo": r["lo"], "hi": r["hi"], "p": r["p"]})
sens = pd.DataFrame(rows)
sens.to_csv(os.path.join(OUT, "eicu_sensitivity.csv"), index=False, float_format="%.3f")
print(sens.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
