#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""单一数据源：把所有要写进论文的数字算出来，存成 results/key_numbers.json
   之后图、Word、校验脚本都只读这个文件，不再手写数字。
"""
import os, json, re, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression, BayesianRidge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

BASE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(BASE, "results")
K = {}

# ---------- 1. 流程 ----------
fm = pd.read_csv(os.path.join(R, "flow_mimic.csv")).iloc[0]
fe = pd.read_csv(os.path.join(R, "flow_eicu.csv")).iloc[0]
K["flow"] = {
    "mimic_ne": int(fm.ne_stays), "mimic_eligible": int(fm.eligible_hour),
    "mimic_sepsis": int(fm.sepsis), "mimic_analytic": int(fm.landmark),
    "mimic_early": int(fm.early), "mimic_delayed": int(fm.delayed),
    "eicu_ne": int(fe.ne_stays), "eicu_convertible": int(fe.convertible),
    "eicu_eligible": int(fe.eligible), "eicu_sepsis": int(fe.sepsis_first),
    "eicu_after_ne6h": int(fe.after_ne6h), "eicu_analytic": int(fe.landmark),
    "eicu_early": int(fe.early), "eicu_delayed": int(fe.delayed),
}
K["flow"]["mimic_excl_no_hour"] = K["flow"]["mimic_ne"] - K["flow"]["mimic_eligible"]
K["flow"]["mimic_excl_no_sepsis"] = K["flow"]["mimic_eligible"] - K["flow"]["mimic_sepsis"]
K["flow"]["mimic_after_filters"] = int(fm.after_filters)
K["flow"]["mimic_excl_filters"] = K["flow"]["mimic_sepsis"] - K["flow"]["mimic_after_filters"]
K["flow"]["mimic_excl_landmark"] = K["flow"]["mimic_after_filters"] - K["flow"]["mimic_analytic"]
K["flow"]["eicu_excl_no_dose"] = K["flow"]["eicu_ne"] - K["flow"]["eicu_convertible"]
K["flow"]["eicu_excl_no_elig"] = K["flow"]["eicu_convertible"] - K["flow"]["eicu_eligible"]
K["flow"]["eicu_excl_no_sepsis"] = K["flow"]["eicu_eligible"] - K["flow"]["eicu_eligible"] + (K["flow"]["eicu_eligible"] - K["flow"]["eicu_sepsis"])
K["flow"]["eicu_excl_ne6h"] = K["flow"]["eicu_sepsis"] - K["flow"]["eicu_after_ne6h"]
K["flow"]["eicu_excl_landmark"] = K["flow"]["eicu_after_ne6h"] - K["flow"]["eicu_analytic"]

# ---------- 2. MIMIC 主分析 ----------
d = pd.read_csv(os.path.join(BASE, "analytic_v2.csv"))
d["age"] = d["age"].clip(upper=90)
for c, n in [("lactate", "log_lact"), ("creatinine", "log_creat"), ("urine_6h", "log_urine")]:
    d[n] = np.log1p(d[c])
CONT = ["age", "dose_t0", "map_t0", "map_mean_6h", "map_min_6h", "frac_map_lt65",
        "hours_ne_before_t0", "hours_icu_before_t0", "log_lact", "log_creat",
        "platelet", "wbc", "hemoglobin", "bun", "gcs_min", "log_urine", "input_6h"]
BIN = ["on_vent", "on_rrt", "on_vasopressin", "urine_missing", "input_missing"]
CAT = ["gender", "admission_type", "insurance", "careunit_t0"]
imp = IterativeImputer(estimator=BayesianRidge(), sample_posterior=True, max_iter=10,
                       random_state=1, skip_complete=True)
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
mdf = pd.DataFrame({"time": d["time_days"], "event": d["event_28d"], "expo": y, "w": w})
cph = CoxPHFitter(); cph.fit(mdf, "time", "event", weights_col="w", robust=True)
s = cph.summary.loc["expo"]
ph = proportional_hazard_test(cph, mdf, time_transform="rank").summary.loc["expo", "p"]
km = {}
for g in (1, 0):
    from lifelines import KaplanMeierFitter
    km[g] = float(KaplanMeierFitter().fit(d.loc[y == g, "time_days"], d.loc[y == g, "event_28d"],
                                          weights=w[y == g]).predict(28))
_sum = open(os.path.join(R, "summary.txt"), encoding="utf-8").read()
_pooled = re.search(r"IPTW 主分析: HR ([\d.]+) \(95% CI ([\d.]+)-([\d.]+)\)", _sum)
K["main"] = {
    "n": int(len(d)), "early": int((y == 1).sum()), "delayed": int((y == 0).sum()),
    "pct_early": round(100 * y.mean(), 1),
    "deaths": int(d["event_28d"].sum()), "pct_death": round(100 * d["event_28d"].mean(), 1),
    "hr": float(_pooled.group(1)),
    "lo": float(_pooled.group(2)),
    "hi": float(_pooled.group(3)),
    "cstat": round(float(roc_auc_score(y, ps)), 2),
    "wmax": round(float(w.max()), 2),
    "ess": int(round((w.sum() ** 2) / (w ** 2).sum())),
    "ph_p": round(float(ph), 2),
    "surv_early": round(100 * km[1], 1), "surv_delayed": round(100 * km[0], 1),
    "restart_early": round(100 * np.average(d.loc[y == 1, "restart_24h"], weights=w[y == 1]), 1),
    "restart_delayed": round(100 * np.average(d.loc[y == 0, "restart_24h"], weights=w[y == 0]), 1),
}
# 未加权与双重稳健从 summary.txt 读
txt = open(os.path.join(R, "summary.txt"), encoding="utf-8").read()
mm = re.findall(r"HR\s+([\d.]+)\s+\(95% CI\s+([\d.]+)-([\d.]+)\)", txt)
K["main"]["unadj_hr"], K["main"]["unadj_lo"], K["main"]["unadj_hi"] = map(float, mm[1])
K["main"]["dr_hr"], K["main"]["dr_lo"], K["main"]["dr_hi"] = map(float, mm[2])
K["main"]["evalue"] = float(re.search(r"E-value = ([\d.]+)", txt).group(1))
K["main"]["evalue_ci"] = round(1 / K["main"]["hi"] + np.sqrt(1 / K["main"]["hi"] * (1 / K["main"]["hi"] - 1)), 2)

# 平衡
bal = pd.read_csv(os.path.join(R, "table_balance.csv"))
K["main"]["max_smd"] = round(float(bal.SMD_weighted.abs().max()), 2)

# ---------- 3. 表 2 / 表 3 ----------
K["sensitivity"] = pd.read_csv(os.path.join(R, "table_sensitivity.csv")).round(3).to_dict("records")
K["subgroups"] = pd.read_csv(os.path.join(R, "table3_subgroups.csv")).round(3).to_dict("records")

# ---------- 4. 表 1 基本情况 ----------
t1 = pd.read_csv(os.path.join(R, "table1.csv"))
K["table1"] = {"n_early": int(y.sum()), "n_delayed": int((1 - y).sum()),
               "max_smd": round(float(t1.SMD.abs().max()), 2)}
sofa = t1[t1["var"] == "sofa"].iloc[0]
K["table1"]["sofa_med"] = round(float(str(sofa[t1.columns[2]]).split(" ")[0]), 0)
K["table1"]["sofa_iqr"] = str(sofa[t1.columns[2]]).split("(")[1].split(")")[0]

# ---------- 5. eICU ----------
eicu_sum = os.path.join(R, "eicu_summary.txt")
if os.path.exists(eicu_sum):
    t = open(eicu_sum, encoding="utf-8").read()
    m = re.search(r"HR ([\d.]+) \(([\d.]+)-([\d.]+)\) p=([\d.]+) \| n=(\d+) 事件=(\d+)", t)
    K["eicu"] = {"hr": float(m.group(1)), "lo": float(m.group(2)), "hi": float(m.group(3)),
                 "p": float(m.group(4)), "n": int(m.group(5)), "events": int(m.group(6))}
    ed = pd.read_csv(os.path.join(BASE, "eicu_analytic.csv"))
    K["eicu"]["early"] = int(ed["early_wean_6h"].sum())
    K["eicu"]["pct_early"] = round(100 * ed["early_wean_6h"].mean(), 1)
    K["eicu"]["deaths"] = int(ed["event_hosp"].sum()) if "event_hosp" in ed else int(ed["event_28d"].sum())
    K["eicu"]["pct_death"] = round(100 * ed["early_wean_6h"].mean() * 0 + 100 * (K["eicu"]["deaths"] / len(ed)), 1)

K["mimic_hosp"] = {}
h = open(os.path.join(R, "mimic_hosp.txt"), encoding="utf-8").read() if os.path.exists(os.path.join(R, "mimic_hosp.txt")) else ""
m = re.search(r"HR ([\d.]+) \(([\d.]+)-([\d.]+)\)", h)
if m:
    K["mimic_hosp"] = {"hr": float(m.group(1)), "lo": float(m.group(2)), "hi": float(m.group(3))}


# ---------- 6. eICU 严重度（APACHE IV） ----------
_eicu_base = "/sessions/funny-sweet-brahmagupta/mnt/eICU-CRD/data/eicu-collaborative-research-database-2.0"
_ap = os.path.join(_eicu_base, "apachePatientResult.csv.gz")
if os.path.exists(_ap) and "eicu" in K:
    import duckdb
    _c = duckdb.connect()
    _q = _c.sql(f"""
      with a as (
        select patientunitstayid, try_cast(apachescore as double) aps,
               try_cast(predictedhospitalmortality as double) phm,
               row_number() over (partition by patientunitstayid order by try_cast(apacheversion as varchar) desc) rn
        from read_csv_auto('{_ap}', header=true)
        where patientunitstayid in (select stay_id from read_csv_auto('{os.path.join(BASE, "eicu_analytic.csv")}', header=true))
      )
      select count(*) n, median(aps) aps_med, quantile_cont(aps,0.25) aps_q1, quantile_cont(aps,0.75) aps_q3,
             median(phm) phm_med, quantile_cont(phm,0.25) phm_q1, quantile_cont(phm,0.75) phm_q3
      from a where rn=1""").df().iloc[0]
    K["eicu"]["apache_n"] = int(_q.n)
    K["eicu"]["aps_med"] = round(float(_q.aps_med), 0)
    K["eicu"]["aps_q1"] = round(float(_q.aps_q1), 0)
    K["eicu"]["aps_q3"] = round(float(_q.aps_q3), 0)
    K["eicu"]["phm_med"] = round(100 * float(_q.phm_med), 1)
    K["eicu"]["phm_q1"] = round(100 * float(_q.phm_q1), 1)
    K["eicu"]["phm_q3"] = round(100 * float(_q.phm_q3), 1)

with open(os.path.join(R, "key_numbers.json"), "w", encoding="utf-8") as f:
    json.dump(K, f, ensure_ascii=False, indent=2)
print(json.dumps({k: K[k] for k in ["flow", "main"]}, ensure_ascii=False, indent=2)[:1200])
print("\nwritten:", os.path.join(R, "key_numbers.json"))
