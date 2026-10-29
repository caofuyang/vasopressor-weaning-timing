#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脓毒性休克首次升压药撤离时机 —— 主分析
输入：analytic_v2.csv（由 03_cohort_v2.sql 生成）
输出：results/*.csv（表 1/2/3）与 results/summary.txt

方法：
  - 暴露 early_wean_6h（landmark t0+6h 前 6 小时内停药或减量≥25%）
  - 结局 28 天全因死亡，时间自 landmark 起算
  - 倾向评分：logistic，预设协变量
  - 权重：ATE 稳定化权重，1%/99% 截尾
  - 主模型：加权 Cox，subject 层聚类稳健方差
  - 缺失：链式多重填补 20 套，Rubin 合并（敏感性：完整病例）
"""
import os, sys, warnings, json
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "results")
os.makedirs(OUT, exist_ok=True)

from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression, BayesianRidge
from sklearn.preprocessing import StandardScaler
from lifelines import CoxPHFitter
from scipy import stats

RNG = np.random.default_rng(20261023)
M_IMP = 20

# ----------------------------------------------------------------------------
# 1. 读入与变量构造
# ----------------------------------------------------------------------------
d0 = pd.read_csv(os.path.join(BASE, "analytic_v2.csv"))
d0["age"] = d0["age"].clip(upper=90)
d0["log_lact"] = np.log1p(d0["lactate"])
d0["log_creat"] = np.log1p(d0["creatinine"])
d0["log_bili"] = np.log1p(d0["bilirubin"])
d0["log_urine"] = np.log1p(d0["urine_6h"])

CONT = ["age", "dose_t0", "map_t0", "map_mean_6h", "map_min_6h", "frac_map_lt65",
        "hours_ne_before_t0", "hours_icu_before_t0", "log_lact", "log_creat",
        "platelet", "wbc", "hemoglobin", "bun", "gcs_min", "log_urine", "input_6h"]
# bilirubin 缺失 37%，不进主模型
CAT = ["gender", "admission_type", "insurance", "careunit_t0"]
BIN = ["on_vent", "on_rrt", "on_vasopressin", "urine_missing", "input_missing"]

def design(df, cont=CONT):
    X = df[cont + BIN].copy()
    dummies = pd.get_dummies(df[CAT], drop_first=True, dtype=float)
    return pd.concat([X.reset_index(drop=True), dummies.reset_index(drop=True)], axis=1)

# ----------------------------------------------------------------------------
# 2. 多重填补
# ----------------------------------------------------------------------------
raw = d0[CONT + BIN].copy()
imp_sets = []
for m in range(M_IMP):
    imp = IterativeImputer(estimator=BayesianRidge(), sample_posterior=True,
                           max_iter=20, random_state=m, skip_complete=True)
    arr = imp.fit_transform(raw)
    df = pd.DataFrame(arr, columns=raw.columns)
    for c in CAT:
        df[c] = d0[c].values
    df["early_wean_6h"] = d0["early_wean_6h"].values
    df["early_wean_2h"] = d0["early_wean_2h"].values
    df["early_wean_1h"] = d0["early_wean_1h"].values
    df["stop6h"] = d0["stop6h"].values
    df["time_days"] = d0["time_days"].values
    df["event_28d"] = d0["event_28d"].values
    imp_sets.append(df)

# ----------------------------------------------------------------------------
# 3. 单个数据集内的分析
# ----------------------------------------------------------------------------
def fit_one(df, exposure="early_wean_6h", weighted=True, trunc=(0.01, 0.99),
            trim_ps=None):
    X = design(df).astype(float)
    X = X.loc[:, X.std() > 0]
    X = pd.DataFrame(StandardScaler().fit_transform(X), columns=X.columns)
    y = df[exposure].values.astype(int)
    lr = LogisticRegression(max_iter=2000, C=1.0)
    lr.fit(X, y)
    ps = lr.predict_proba(X)[:, 1]
    if trim_ps is not None:
        lo, hi = trim_ps
        keep = (ps >= lo) & (ps <= hi)
    else:
        keep = np.ones(len(ps), bool)
    dfk, Xk, yk, psk = df[keep], X[keep], y[keep], ps[keep]
    # 稳定化 ATE 权重
    p_treat = yk.mean()
    w = np.where(yk == 1, p_treat / psk, (1 - p_treat) / (1 - psk))
    if trunc is not None:
        lo, hi = np.quantile(w, trunc)
    else:
        lo, hi = None, None
    if weighted and trunc is not None:
        w = np.clip(w, lo, hi)
    mdf = pd.DataFrame({
        "time": dfk["time_days"].values,
        "event": dfk["event_28d"].values,
        "expo": yk,
        "w": w if weighted else np.ones(len(w)),
    })
    cph = CoxPHFitter()
    cph.fit(mdf, "time", "event", weights_col="w", robust=True)
    s = cph.summary.loc["expo"]
    return dict(
        coef=s["coef"], se=s["se(coef)"], hr=s["exp(coef)"],
        lo=s["exp(coef) lower 95%"], hi=s["exp(coef) upper 95%"], p=s["p"],
        n=int(keep.sum()), n_event=int(mdf.event.sum()),
        w_max=float(w.max()), w_min=float(w.min()),
        ps_min=float(psk.min()), ps_max=float(psk.max()),
        ps_lo=float(np.quantile(ps, 0.01)), ps_hi=float(np.quantile(ps, 0.99)),
    ), X, y, ps, df

# ----------------------------------------------------------------------------
# 4. Rubin 合并
# ----------------------------------------------------------------------------
def pool(fits):
    m = len(fits)
    Q = np.array([f["coef"] for f in fits])
    U = np.array([f["se"] ** 2 for f in fits])
    Qbar = Q.mean(); Ubar = U.mean(); B = Q.var(ddof=1)
    T = Ubar + (1 + 1 / m) * B
    se = np.sqrt(T)
    lam = (1 + 1 / m) * B / T
    df_old = (m - 1) / lam ** 2 if lam > 0 else np.inf
    df_com = fits[0]["n"] - 2
    df_new = df_com * df_old / (df_com + df_old) if np.isfinite(df_old) else df_com
    tcrit = stats.t.ppf(0.975, df_new)
    z = Qbar / se
    p = 2 * stats.t.sf(abs(z), df_new)
    return dict(coef=Qbar, se=se, hr=np.exp(Qbar),
                lo=np.exp(Qbar - tcrit * se), hi=np.exp(Qbar + tcrit * se),
                p=p, fmi=lam, df=df_new)

def fmt(r, label):
    return (f"{label}: HR {r['hr']:.3f} (95% CI {r['lo']:.3f}-{r['hi']:.3f}), "
            f"p={r['p']:.4f}, FMI={r.get('fmi', float('nan')):.3f}")

# ----------------------------------------------------------------------------
# 5. 主分析
# ----------------------------------------------------------------------------
fits, balances = [], []
for df in imp_sets:
    f, X, y, ps, dfk = fit_one(df)
    fits.append(f)
    balances.append((dfk, ps))
main = pool(fits)
print("主分析（IPTW，20 套填补合并）")
print(fmt(main, "  early_wean_6h"))

# 未加权
adj_fits = []
for df in imp_sets:
    f, *_ = fit_one(df, weighted=False, trunc=None)
    adj_fits.append(f)
unadj = pool(adj_fits)
print(fmt(unadj, "  未加权（仅看量级）"))

# 双重稳健：加权 + 协变量
def fit_dr(df):
    X = design(df).astype(float)
    X = X.loc[:, X.std() > 0]
    X = pd.DataFrame(StandardScaler().fit_transform(X), columns=X.columns)
    y = df["early_wean_6h"].values.astype(int)
    lr = LogisticRegression(max_iter=2000).fit(X, y)
    ps = lr.predict_proba(X)[:, 1]
    p_treat = y.mean()
    w = np.where(y == 1, p_treat / ps, (1 - p_treat) / (1 - ps))
    lo, hi = np.quantile(w, (0.01, 0.99)); w = np.clip(w, lo, hi)
    mdf = pd.DataFrame({"time": df["time_days"].values, "event": df["event_28d"].values,
                        "expo": y, "w": w})
    covs = [c for c in ["age", "dose_t0", "map_mean_6h", "gcs_min", "log_lact", "log_creat"]
            if c in df.columns]
    for c in covs:
        mdf[c] = (df[c] - df[c].mean()) / df[c].std()
    cph = CoxPHFitter()
    cph.fit(mdf, "time", "event", weights_col="w", robust=True)
    s = cph.summary.loc["expo"]
    return dict(coef=s["coef"], se=s["se(coef)"], n=len(mdf), n_event=int(mdf.event.sum()))

dr = pool([fit_dr(df) for df in imp_sets])
print(fmt(dr, "  双重稳健（加权+协变量）"))

# ----------------------------------------------------------------------------
# 6. 平衡诊断（首套填补，平均）
# ----------------------------------------------------------------------------
dfb, psb = balances[0]
Xb = design(dfb).astype(float)
Xb = Xb.loc[:, Xb.std() > 0]
yb = dfb["early_wean_6h"].values.astype(int)
p_treat = yb.mean()
wb = np.where(yb == 1, p_treat / psb, (1 - p_treat) / (1 - psb))
lo, hi = np.quantile(wb, (0.01, 0.99)); wb = np.clip(wb, lo, hi)

def smd(a, b):
    return (a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)

def wsmd(a, b, wa, wb_):
    ma, mb = np.average(a, weights=wa), np.average(b, weights=wb_)
    va = np.average((a - ma) ** 2, weights=wa)
    vb = np.average((b - mb) ** 2, weights=wb_)
    return (ma - mb) / np.sqrt((va + vb) / 2)

rows = []
for c in Xb.columns:
    a = Xb.loc[yb == 1, c].values; b = Xb.loc[yb == 0, c].values
    rows.append((c, smd(a, b), wsmd(a, b, wb[yb == 1], wb[yb == 0])))
bal = pd.DataFrame(rows, columns=["covariate", "SMD_unweighted", "SMD_weighted"])
bal.to_csv(os.path.join(OUT, "table_balance.csv"), index=False, float_format="%.4f")
print("\n平衡：最大 |SMD| 未加权 %.3f → 加权 %.3f；|SMD|>0.1 的变量 %d 个"
      % (bal.SMD_unweighted.abs().max(), bal.SMD_weighted.abs().max(),
         (bal.SMD_weighted.abs() > 0.1).sum()))

# ----------------------------------------------------------------------------
# 7. 敏感性分析
# ----------------------------------------------------------------------------
sens = []
for expo in ["early_wean_6h", "early_wean_2h", "early_wean_1h", "stop6h"]:
    ff = [fit_one(df, exposure=expo)[0] for df in imp_sets]
    sens.append((expo, pool(ff)))
# 截尾更狠
ff = [fit_one(df, trunc=(0.025, 0.975))[0] for df in imp_sets]
sens.append(("truncate_2.5_97.5", pool(ff)))
# 共同支持域
ff = [fit_one(df, trim_ps=(0.05, 0.95))[0] for df in imp_sets]
sens.append(("common_support_0.05_0.95", pool(ff)))

# 完整病例（含 bilirubin）
cc = d0.dropna(subset=[c for c in CONT if c != "log_bili"]).copy()
cc_f = []
for m in range(5):
    f, *_ = fit_one(cc)
    cc_f.append(f)
sens.append(("complete_case", pool(cc_f)))

sensdf = pd.DataFrame([{"analysis": k, "HR": v["hr"], "lo": v["lo"], "hi": v["hi"],
                        "p": v["p"], "n": f"~{fits[0]['n']}"} for k, v in sens])
sensdf.to_csv(os.path.join(OUT, "table_sensitivity.csv"), index=False, float_format="%.3f")
print("\n敏感性分析")
print(sensdf.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# ----------------------------------------------------------------------------
# 8. 异质性（按 t0 剂量分层）
# ----------------------------------------------------------------------------
het_rows = []
for band, lab in [((0, 0.05), "dose<0.05"), ((0.05, 0.101), "dose>=0.05")]:
    sub = d0[(d0.dose_t0 >= band[0]) & (d0.dose_t0 < band[1])]
    if len(sub) < 100:
        continue
    fs = []
    for m in range(5):
        df = imp_sets[m]
        idx = df.index.isin(sub.index)
        if idx.sum() < 100:
            continue
        f, *_ = fit_one(df.loc[idx].reset_index(drop=True))
        fs.append(f)
    if len(fs) >= 2:
        r = pool(fs)
        het_rows.append({"stratum": lab, "n": len(sub), "HR": r["hr"],
                         "lo": r["lo"], "hi": r["hi"], "p": r["p"]})
het = pd.DataFrame(het_rows)
het.to_csv(os.path.join(OUT, "table_heterogeneity.csv"), index=False, float_format="%.3f")
print("\n异质性（按 t0 剂量）")
print(het.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# ----------------------------------------------------------------------------
# 9. E-value（保护效应）
# ----------------------------------------------------------------------------
hr = main["hr"]
hr_ = 1 / hr if hr < 1 else hr
evalue = hr_ + np.sqrt(hr_ * (hr_ - 1))
print(f"\nE-value（针对 HR={hr:.3f}）= {evalue:.2f}")

with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as fh:
    fh.write(fmt(main, "IPTW 主分析") + "\n")
    fh.write(fmt(unadj, "未加权") + "\n")
    fh.write(fmt(dr, "双重稳健") + "\n")
    fh.write(f"E-value = {evalue:.2f}\n")
    fh.write(f"n(首套) = {fits[0]['n']}, events = {fits[0]['n_event']}\n")
    fh.write(f"weights max = {fits[0]['w_max']:.2f}\n")
print("\n结果已写入: " + OUT)
