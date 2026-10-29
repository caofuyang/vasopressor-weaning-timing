-- ============================================================================
-- eICU-CRD 外部验证 —— 队列构建（与 MIMIC 03_cohort_v2.sql 同口径）
-- 运行：& E:\数据库\duckdb\duckdb.exe "C:/Users/caofu/tte/eicu.duckdb" -c ".read 'C:/Users/caofu/tte/07_eicu_cohort.sql'"
-- 产出：C:/Users/caofu/tte/eicu_analytic.csv
--
-- 与 MIMIC 的差异（必须在论文里写清楚）：
--   1. eICU 无 28 天随访，主结局改为「院内死亡」，时间自 landmark 至出院/死亡
--   2. 剂量以 mcg/min 为主，用 admissionweight 换算为 µg/kg/min；
--      体重缺失者不纳入（MIMIC 中体重由 inputevents 提供）
--   3. 无小时级结构化输注，用 infusionoffset 的观测值前向携带，
--      间隔 >180 分钟视为中断
--   4. 无 GCS 时点值，改用插管/透析/血管加压素等二分类协变量
-- ============================================================================

CREATE OR REPLACE MACRO base() AS 'E:/数据库/eICU-CRD/data/eicu-collaborative-research-database-2.0';

SET threads = 6;
SET memory_limit = '10GB';
.timer on

-- ----------------------------------------------------------------------------
-- STEP 0  患者表
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE pat AS
SELECT patientunitstayid AS stay_id,
       patienthealthsystemstayid AS pt_id,
       gender,
       CASE WHEN age LIKE '>%' THEN 90 ELSE TRY_CAST(age AS INT) END AS age,
       unittype, hospitalid,
       TRY_CAST(admissionweight AS DOUBLE) AS weight,
       unitstaytype,
       CAST(unitdischargeoffset AS INT)     AS unit_out,
       CAST(hospitaldischargeoffset AS INT) AS hosp_out,
       hospitaldischargestatus AS hosp_status
FROM read_csv_auto(base() || '/patient.csv.gz', header = true);

-- ----------------------------------------------------------------------------
-- STEP 0b  脓毒症定义（eICU diagnosis 表）
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE sepsis_stay AS
SELECT DISTINCT patientunitstayid AS stay_id
FROM read_csv_auto(base() || '/diagnosis.csv.gz', header = true)
WHERE icd9code LIKE '038%'
   OR lower(diagnosisstring) LIKE '%sepsis%'
   OR lower(diagnosisstring) LIKE '%septic%';

-- ----------------------------------------------------------------------------
-- STEP 1  去甲肾上腺素记录 → µg/kg/min
-- 只取单位明确为 mcg/min 或 mcg/kg/min 的记录；ml/hr 无法换算，丢弃
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE ne_raw AS
SELECT i.patientunitstayid AS stay_id,
       CAST(i.infusionoffset AS INT) AS t_off,
       TRY_CAST(i.drugrate AS DOUBLE) AS rate,
       lower(i.drugname) AS dn
FROM read_csv_auto(base() || '/infusionDrug.csv.gz', header = true) i
WHERE lower(i.drugname) LIKE '%norepinephrine%'
  AND TRY_CAST(i.drugrate AS DOUBLE) IS NOT NULL
  AND (lower(i.drugname) LIKE '%mcg/min%' OR lower(i.drugname) LIKE '%mcg/kg/min%')
  AND lower(i.drugname) NOT LIKE '%ml/hr%';

CREATE OR REPLACE TABLE ne_dose AS
SELECT r.stay_id, r.t_off,
       CASE WHEN r.dn LIKE '%mcg/kg/min%' THEN r.rate
            WHEN p.weight > 20 THEN r.rate / p.weight END AS dose
FROM ne_raw r
JOIN pat p USING (stay_id)
WHERE (r.dn LIKE '%mcg/kg/min%' OR p.weight > 20)
  AND r.t_off >= 0;

-- ----------------------------------------------------------------------------
-- STEP 2  逐小时剂量（前向携带，间隔 >180 分钟视为缺失）
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE hrs AS
SELECT s.stay_id, t.hr, t.hr * 60 AS t_min
FROM (SELECT stay_id, max(t_off) AS mx, min(t_off) AS mn FROM ne_dose GROUP BY stay_id) s,
     unnest(range(0, CAST(ceil(s.mx / 60.0) AS INT) + 1)) AS t(hr);

CREATE OR REPLACE TABLE ne_hourly AS
SELECT h.stay_id, h.hr, h.t_min,
       CASE WHEN d.t_off IS NOT NULL AND h.t_min - d.t_off <= 180 THEN d.dose END AS dose,
       d.t_off AS last_obs
FROM hrs h
ASOF LEFT JOIN ne_dose d
  ON h.stay_id = d.stay_id AND h.t_min >= d.t_off;

-- ----------------------------------------------------------------------------
-- STEP 3  逐小时 MAP
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE map_hourly AS
SELECT v.patientunitstayid AS stay_id,
       CAST(floor(v.observationoffset / 60.0) AS INT) AS hr,
       avg(coalesce(v.systemicmean,
                    (v.systemicsystolic + v.systemicdiastolic) / 2.0,
                    (v.pasystolic + v.padiastolic) / 2.0)) AS map_mean,
       min(coalesce(v.systemicmean,
                    (v.systemicsystolic + v.systemicdiastolic) / 2.0,
                    (v.pasystolic + v.padiastolic) / 2.0)) AS map_min
FROM read_csv_auto(base() || '/vitalPeriodic.csv.gz', header = true) v
WHERE v.patientunitstayid IN (SELECT DISTINCT stay_id FROM ne_dose)
  AND v.observationoffset >= 0
GROUP BY 1, 2;

-- ----------------------------------------------------------------------------
-- STEP 4  t0：与 MIMIC 相同的稳定窗
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE hrs2 AS
SELECT h.stay_id, h.hr, h.dose, m.map_mean, m.map_min,
       lead(h.dose)     OVER w AS dose_next,
       lead(m.map_mean) OVER w AS map_next,
       max(h.dose)      OVER (PARTITION BY h.stay_id ORDER BY h.hr
                              ROWS BETWEEN 6 PRECEDING AND 1 PRECEDING) AS prev6_max
FROM ne_hourly h
LEFT JOIN map_hourly m USING (stay_id, hr)
WINDOW w AS (PARTITION BY h.stay_id ORDER BY h.hr);

CREATE OR REPLACE TABLE t0_cases AS
SELECT stay_id, hr AS t0_hr, t0_hr * 60 AS t0_min, dose AS dose_t0, map_mean AS map_t0
FROM hrs2
WHERE hr >= 6
  AND dose BETWEEN 0.0001 AND 0.10
  AND map_mean >= 65 AND map_next >= 65
  AND dose_next <= 0.10
  AND (prev6_max IS NULL OR dose >= 0.9 * prev6_max);

CREATE OR REPLACE TABLE ne_start AS
SELECT stay_id, min(t_off) AS ne_start FROM ne_dose GROUP BY stay_id;

CREATE OR REPLACE TABLE cohort_base AS
SELECT c.*, p.pt_id, p.age, p.gender, p.weight, p.unittype, p.hospitalid,
       p.unit_out, p.hosp_out, p.hosp_status,
       ns.ne_start
FROM t0_cases c
JOIN pat p USING (stay_id)
JOIN ne_start ns USING (stay_id)
JOIN sepsis_stay ss USING (stay_id)
QUALIFY row_number() OVER (PARTITION BY c.stay_id ORDER BY c.t0_hr) = 1;

-- 首次 stay + NE 已用 ≥6h
CREATE OR REPLACE TABLE cohort_v2 AS
WITH a AS (
  SELECT cb.*, row_number() OVER (PARTITION BY cb.pt_id ORDER BY cb.t0_min) AS r_subj,
         (cb.t0_min - cb.ne_start) / 60.0 AS hours_ne_before_t0
  FROM cohort_base cb
)
SELECT * FROM a WHERE r_subj = 1 AND hours_ne_before_t0 >= 6;

-- ----------------------------------------------------------------------------
-- STEP 5  暴露：6 小时内停药或减量 ≥25%
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE exposure_v2 AS
WITH lastseg AS (
  SELECT c.stay_id, max(d.t_off) AS last_obs_all
  FROM cohort_v2 c JOIN ne_dose d USING (stay_id)
  GROUP BY c.stay_id
),
win AS (
  SELECT c.stay_id, c.t0_min, c.dose_t0,
         min(CASE WHEN h.hr > c.t0_hr AND h.hr <= c.t0_hr + 6 THEN h.dose END) AS m6,
         min(CASE WHEN h.hr > c.t0_hr AND h.hr <= c.t0_hr + 2 THEN h.dose END) AS m2,
         min(CASE WHEN h.hr > c.t0_hr AND h.hr <= c.t0_hr + 1 THEN h.dose END) AS m1
  FROM cohort_v2 c
  LEFT JOIN hrs2 h USING (stay_id)
  GROUP BY c.stay_id, c.t0_min, c.dose_t0
)
SELECT w.*, l.last_obs_all,
       CASE WHEN l.last_obs_all <= w.t0_min + 360 THEN 1 ELSE 0 END AS stop6h,
       CASE WHEN l.last_obs_all <= w.t0_min + 360
                 OR (w.m6 IS NOT NULL AND w.m6 < 0.75 * w.dose_t0) THEN 1 ELSE 0 END AS early_wean_6h,
       CASE WHEN l.last_obs_all <= w.t0_min + 360
                 OR (w.m2 IS NOT NULL AND w.m2 < 0.75 * w.dose_t0) THEN 1 ELSE 0 END AS early_wean_2h
FROM win w LEFT JOIN lastseg l USING (stay_id);

-- ----------------------------------------------------------------------------
-- STEP 6  landmark = t0 + 6h
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE landmark AS
SELECT c.*, e.m6, e.stop6h, e.early_wean_6h, e.early_wean_2h,
       c.t0_min + 360 AS lm_min
FROM cohort_v2 c JOIN exposure_v2 e USING (stay_id)
WHERE c.hosp_out > c.t0_min + 360
  AND (c.unit_out IS NULL OR c.unit_out > c.t0_min + 360);

-- ----------------------------------------------------------------------------
-- STEP 7  协变量
-- ----------------------------------------------------------------------------
-- 7.1 检验（t0 前 24h 内最后一次）
CREATE OR REPLACE TABLE labs_v2 AS
WITH raw AS (
  SELECT l.patientunitstayid AS stay_id, lower(l.labname) AS labname,
         TRY_CAST(l.labresult AS DOUBLE) AS v, l.labresultoffset AS t_off
  FROM read_csv_auto(base() || '/lab.csv.gz', header = true) l
  JOIN landmark k ON k.stay_id = l.patientunitstayid
  WHERE l.labresultoffset BETWEEN k.t0_min - 1440 AND k.t0_min
    AND lower(l.labname) IN ('creatinine','wbc x 1000','hgb','platelets x 1000',
                             'bun','total bilirubin','lactate')
),
ranked AS (
  SELECT *, row_number() OVER (PARTITION BY stay_id, labname ORDER BY t_off DESC) AS rn
  FROM raw
)
SELECT stay_id,
       max(CASE WHEN labname = 'lactate' AND rn = 1 THEN v END) AS lactate,
       max(CASE WHEN labname = 'creatinine' AND rn = 1 THEN v END) AS creatinine,
       max(CASE WHEN labname = 'wbc x 1000' AND rn = 1 THEN v END) AS wbc,
       max(CASE WHEN labname = 'hgb' AND rn = 1 THEN v END) AS hemoglobin,
       max(CASE WHEN labname = 'platelets x 1000' AND rn = 1 THEN v END) AS platelet,
       max(CASE WHEN labname = 'bun' AND rn = 1 THEN v END) AS bun,
       max(CASE WHEN labname = 'total bilirubin' AND rn = 1 THEN v END) AS bilirubin
FROM ranked GROUP BY stay_id;

-- 7.2 血管加压素（landmark 时仍在用）
CREATE OR REPLACE TABLE vp_use AS
SELECT DISTINCT k.stay_id
FROM landmark k
JOIN read_csv_auto(base() || '/infusionDrug.csv.gz', header = true) i
  ON i.patientunitstayid = k.stay_id
WHERE lower(i.drugname) LIKE '%vasopressin%'
  AND CAST(i.infusionoffset AS INT) <= k.t0_min + 360;

-- 7.3 机械通气（landmark 时）
CREATE OR REPLACE TABLE vent_use AS
SELECT DISTINCT k.stay_id
FROM landmark k
JOIN read_csv_auto(base() || '/respiratoryCare.csv.gz', header = true) r
  ON r.patientunitstayid = k.stay_id
WHERE TRY_CAST(r.ventstartoffset AS INT) > 0
  AND TRY_CAST(r.ventstartoffset AS INT) <= k.t0_min + 360
  AND (TRY_CAST(r.ventendoffset AS INT) IS NULL
       OR TRY_CAST(r.ventendoffset AS INT) = 0
       OR TRY_CAST(r.ventendoffset AS INT) >= k.t0_min + 360);

-- 7.4 肾脏替代（landmark 前 24h 内）
CREATE OR REPLACE TABLE rrt_use AS
SELECT DISTINCT k.stay_id
FROM landmark k
JOIN read_csv_auto(base() || '/treatment.csv.gz', header = true) t
  ON t.patientunitstayid = k.stay_id
WHERE lower(t.treatmentstring) LIKE '%dialysis%'
  AND TRY_CAST(t.treatmentoffset AS INT) BETWEEN k.t0_min - 1440 AND k.t0_min + 360;

-- 7.5 尿量（t0 前 6h）
CREATE OR REPLACE TABLE urine_6h AS
SELECT k.stay_id, sum(TRY_CAST(io.cellvaluenumeric AS DOUBLE)) AS urine_6h, count(*) AS n_rows
FROM landmark k
JOIN read_csv_auto(base() || '/intakeOutput.csv.gz', header = true) io
  ON io.patientunitstayid = k.stay_id
WHERE lower(io.celllabel) LIKE '%urine%'
  AND TRY_CAST(io.intakeoutputoffset AS INT) BETWEEN k.t0_min - 360 AND k.t0_min
GROUP BY k.stay_id;

-- 7.6 前 6h MAP 负荷
CREATE OR REPLACE TABLE burden AS
SELECT k.stay_id,
       avg(h.map_mean) AS map_mean_6h,
       min(h.map_min)  AS map_min_6h,
       avg(CASE WHEN h.map_mean < 65 THEN 1.0 ELSE 0.0 END) AS frac_map_lt65
FROM landmark k JOIN hrs2 h USING (stay_id)
WHERE h.hr BETWEEN k.t0_hr - 6 AND k.t0_hr - 1
GROUP BY k.stay_id;

-- ----------------------------------------------------------------------------
-- STEP 8  合并导出
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE analytic AS
SELECT k.stay_id, k.pt_id, k.hospitalid, k.unittype,
       k.age, k.gender, k.weight,
       k.early_wean_6h, k.early_wean_2h, k.stop6h,
       k.dose_t0, k.map_t0, k.hours_ne_before_t0,
       bd.map_mean_6h, bd.map_min_6h, bd.frac_map_lt65,
       lb.lactate, lb.creatinine, lb.wbc, lb.hemoglobin, lb.platelet, lb.bun, lb.bilirubin,
       CASE WHEN vp.stay_id IS NOT NULL THEN 1 ELSE 0 END AS on_vasopressin,
       CASE WHEN ve.stay_id IS NOT NULL THEN 1 ELSE 0 END AS on_vent,
       CASE WHEN rr.stay_id IS NOT NULL THEN 1 ELSE 0 END AS on_rrt,
       u.urine_6h,
       -- 结局：院内死亡；时间自 landmark 至出院/死亡
       CASE WHEN k.hosp_status = 'Expired' THEN 1 ELSE 0 END AS event_hosp,
       (k.hosp_out - k.lm_min) / 1440.0 AS time_days,
       (k.hosp_out - k.lm_min) / 1440.0 AS hosp_days_after_lm
FROM landmark k
LEFT JOIN labs_v2   lb USING (stay_id)
LEFT JOIN burden    bd USING (stay_id)
LEFT JOIN urine_6h  u  USING (stay_id)
LEFT JOIN vp_use    vp USING (stay_id)
LEFT JOIN vent_use  ve USING (stay_id)
LEFT JOIN rrt_use   rr USING (stay_id);

COPY analytic TO 'C:/Users/caofu/tte/eicu_analytic.csv' (HEADER, DELIMITER ',');

-- ----------------------------------------------------------------------------
-- STEP 9  核对
-- ----------------------------------------------------------------------------
SELECT '--- 流程 ---' AS section;
SELECT 'NE 记录（可用剂量）' AS step, count(*) AS n FROM ne_dose
UNION ALL SELECT '有 t0 的 stay', count(*) FROM t0_cases
UNION ALL SELECT '首次 stay + NE≥6h', count(*) FROM cohort_v2
UNION ALL SELECT 'landmark 后', count(*) FROM landmark
UNION ALL SELECT '最终入分析', count(*) FROM analytic;

SELECT '--- 暴露与结局 ---' AS section;
SELECT count(*) n, round(100.0*avg(early_wean_6h),1) pct_6h,
       round(100.0*avg(early_wean_2h),1) pct_2h,
       round(100.0*avg(stop6h),1) pct_stop,
       sum(event_hosp) deaths, round(100.0*avg(event_hosp),1) pct_death,
       round(median(time_days),1) med_days
FROM analytic;

SELECT '--- 缺失率 ---' AS section;
SELECT 'lactate' v, round(100.0*avg(CASE WHEN lactate IS NULL THEN 1 ELSE 0 END),1) pct FROM analytic
UNION ALL SELECT 'creatinine', round(100.0*avg(CASE WHEN creatinine IS NULL THEN 1 ELSE 0 END),1) FROM analytic
UNION ALL SELECT 'urine_6h', round(100.0*avg(CASE WHEN urine_6h IS NULL THEN 1 ELSE 0 END),1) FROM analytic
UNION ALL SELECT 'weight', round(100.0*avg(CASE WHEN weight IS NULL THEN 1 ELSE 0 END),1) FROM analytic;

SELECT '--- 粗率 ---' AS section;
SELECT early_wean_6h, count(*) n, round(100.0*avg(event_hosp),1) pct_death,
       round(avg(dose_t0),4) mean_dose
FROM analytic GROUP BY 1 ORDER BY 1;
