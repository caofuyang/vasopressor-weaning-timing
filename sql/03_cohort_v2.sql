-- ============================================================================
-- MIMIC-IV 撤药时机 TTE —— 第 3 步（v2 修正版）
-- 相对 02_cohort_landmark.sql 的六处改动：
--   1. 每位患者只保留首次 ICU stay（原 969 例重复 stay 违反独立性）
--   2. 入组要求 t0 前 NE 已用 ≥6 小时（原 35% 在 2 小时内即入组）
--   3. 暴露 = t0 后 6 小时内「停药」或「减量 ≥25%」（原只看第 1 小时，
--      导致 937 例真正停药者被判为未撤）；1h/2h 窗保留作敏感性分析
--   4. landmark 移到 t0+6h，剔除 6h 内死亡/转出者（原从 t0 起算，存在
--      不朽时间偏倚）
--   5. 尿量与入量不再 coalesce 成 0，缺失保留为 NULL 并附缺失指示变量
--   6. 检验值改为 t0 前最后一次，且显式导出 time_to_event / event_28d
-- 产出：analytic_v2.csv（可直接导入 SPSS 做倾向评分与加权 Cox）
-- 运行：& E:\数据库\duckdb\duckdb.exe "C:/Users/caofu/tte/tte.duckdb" -c ".read 'C:/Users/caofu/tte/03_cohort_v2.sql'"
-- 前置：tte.duckdb 中已有 02 建的 ne_dose / hrs / cohort_base / burden /
--       support / stay_win / outcome 等表（本脚本只重算队列与部分协变量）
-- ============================================================================

CREATE OR REPLACE MACRO base() AS 'E:/数据库/MIMIC-IV/data/mimic-iv-3.1';

SET threads = 6;
SET memory_limit = '4GB';
.timer on

-- ----------------------------------------------------------------------------
-- STEP 0  脓毒症定义（ICD 诊断）
--   ICD-9 : 038.x、995.91/995.92、785.52
--   ICD-10: A40.x/A41.x、R65.2x
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE sepsis_hadm AS
SELECT DISTINCT hadm_id
FROM read_csv_auto(base() || '/hosp/diagnoses_icd.csv.gz', header = true)
WHERE (TRY_CAST(icd_version AS INT) = 9
       AND (icd_code LIKE '038%' OR icd_code IN ('99591','99592','78552')))
   OR (TRY_CAST(icd_version AS INT) = 10
       AND (icd_code LIKE 'A40%' OR icd_code LIKE 'A41%' OR icd_code LIKE 'R652%'));

-- ----------------------------------------------------------------------------
-- STEP 1  队列：首次 stay + NE ≥6h
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE cohort_v2 AS
WITH a AS (
  SELECT cb.*, ns.ne_start,
         row_number() OVER (PARTITION BY cb.subject_id ORDER BY cb.intime) AS r_subj,
         date_diff('minute', ns.ne_start, cb.t0) / 60.0 AS hours_ne_before_t0
  FROM cohort_base cb
  JOIN sepsis_hadm sp USING (hadm_id)
  JOIN (SELECT stay_id, min(starttime) AS ne_start FROM ne_dose GROUP BY stay_id) ns
    USING (stay_id)
)
SELECT * FROM a WHERE r_subj = 1 AND hours_ne_before_t0 >= 6;

-- ----------------------------------------------------------------------------
-- STEP 2  暴露：6 小时窗内停药或减量 ≥25%
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE exposure_v2 AS
WITH lastseg AS (
  SELECT c.stay_id, max(n.endtime) AS end_by_window
  FROM cohort_v2 c
  JOIN ne_dose n USING (stay_id)
  WHERE n.starttime <= c.t0 + INTERVAL 6 HOUR
  GROUP BY c.stay_id
),
win AS (
  SELECT c.stay_id, c.dose_t0, c.t0, c.t0_hr,
         min(CASE WHEN h.hr > c.t0_hr AND h.hr <= c.t0_hr + 1 THEN h.ne_dose END) AS m1,
         min(CASE WHEN h.hr > c.t0_hr AND h.hr <= c.t0_hr + 2 THEN h.ne_dose END) AS m2,
         min(CASE WHEN h.hr > c.t0_hr AND h.hr <= c.t0_hr + 6 THEN h.ne_dose END) AS m6,
         max(CASE WHEN h.hr > c.t0_hr AND h.hr <= c.t0_hr + 24 THEN h.ne_dose END) AS mx24,
         count(CASE WHEN h.hr > c.t0_hr AND h.hr <= c.t0_hr + 6 AND h.ne_dose IS NOT NULL THEN 1 END) AS hours_obs_6h
  FROM cohort_v2 c
  LEFT JOIN hrs h ON h.stay_id = c.stay_id AND h.hr > c.t0_hr AND h.hr <= c.t0_hr + 24
  GROUP BY c.stay_id, c.dose_t0, c.t0, c.t0_hr
)
SELECT w.*, l.end_by_window,
       CASE WHEN l.end_by_window <= w.t0 + INTERVAL 6 HOUR THEN 1 ELSE 0 END AS stop6h,
       CASE WHEN l.end_by_window <= w.t0 + INTERVAL 2 HOUR THEN 1 ELSE 0 END AS stop2h,
       CASE WHEN l.end_by_window <= w.t0 + INTERVAL 1 HOUR THEN 1 ELSE 0 END AS stop1h,
       CASE WHEN l.end_by_window <= w.t0 + INTERVAL 6 HOUR
                 OR (w.m6 IS NOT NULL AND w.m6 < 0.75 * w.dose_t0) THEN 1 ELSE 0 END AS early_wean_6h,
       CASE WHEN l.end_by_window <= w.t0 + INTERVAL 2 HOUR
                 OR (w.m2 IS NOT NULL AND w.m2 < 0.75 * w.dose_t0) THEN 1 ELSE 0 END AS early_wean_2h,
       CASE WHEN l.end_by_window <= w.t0 + INTERVAL 1 HOUR
                 OR (w.m1 IS NOT NULL AND w.m1 < 0.75 * w.dose_t0) THEN 1 ELSE 0 END AS early_wean_1h,
       CASE WHEN w.mx24 > 1.2 * w.dose_t0 THEN 1 ELSE 0 END AS restart_24h
FROM win w
LEFT JOIN lastseg l USING (stay_id);

-- ----------------------------------------------------------------------------
-- STEP 3  landmark = t0 + 6h，剔除窗口内死亡或转出
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE landmark AS
SELECT c.stay_id, c.subject_id, c.hadm_id, c.intime, c.t0, c.t0_hr,
       c.first_careunit, c.hours_ne_before_t0, c.map_t0,
       c.t0 + INTERVAL 6 HOUR AS t_landmark,
       e.* EXCLUDE (stay_id, t0, t0_hr),
       o.death_time, o.dischtime, o.outtime, o.admission_type, o.insurance
FROM cohort_v2 c
JOIN exposure_v2 e USING (stay_id)
JOIN outcome o USING (stay_id)
WHERE (o.death_time IS NULL OR o.death_time > c.t0 + INTERVAL 6 HOUR)
  AND (o.outtime   IS NULL OR o.outtime   > c.t0 + INTERVAL 6 HOUR);

-- ----------------------------------------------------------------------------
-- STEP 4  检验值：t0 前 24 小时内「最后一次」（原用 max，方向不一致）
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE labs_v2 AS
WITH raw AS (
  SELECT i.stay_id, CAST(le.itemid AS INT) AS itemid,
         TRY_CAST(le.valuenum AS DOUBLE) AS v, le.charttime
  FROM read_csv_auto(base() || '/hosp/labevents.csv.gz', header = true) le
  JOIN cohort_v2 i ON i.hadm_id = le.hadm_id
  WHERE TRY_CAST(le.itemid AS INT) IN (50813, 50912, 50885, 51265, 51301, 51222, 51006)
    AND le.charttime BETWEEN i.t0 - INTERVAL 24 HOUR AND i.t0
),
ranked AS (
  SELECT *, row_number() OVER (PARTITION BY stay_id, itemid ORDER BY charttime DESC) AS rn
  FROM raw
)
SELECT stay_id,
       max(CASE WHEN itemid = 50813 AND rn = 1 THEN v END) AS lactate,
       max(CASE WHEN itemid = 50912 AND rn = 1 THEN v END) AS creatinine,
       max(CASE WHEN itemid = 50885 AND rn = 1 THEN v END) AS bilirubin,
       max(CASE WHEN itemid = 51265 AND rn = 1 THEN v END) AS platelet,
       max(CASE WHEN itemid = 51301 AND rn = 1 THEN v END) AS wbc,
       max(CASE WHEN itemid = 51222 AND rn = 1 THEN v END) AS hemoglobin,
       max(CASE WHEN itemid = 51006 AND rn = 1 THEN v END) AS bun
FROM ranked
GROUP BY stay_id;

-- ----------------------------------------------------------------------------
-- STEP 5  液体：缺失保留 NULL（原 coalesce 把 21% 变成「无尿」）
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE fluids_v2 AS
WITH urine_items AS (
  SELECT itemid FROM read_csv_auto(base() || '/icu/d_items.csv.gz', header = true)
  WHERE label ILIKE '%urine%' OR label ILIKE '%foley%'
),
u AS (
  SELECT o.stay_id, sum(TRY_CAST(o.value AS DOUBLE)) AS urine_6h, count(*) AS n_urine_rows
  FROM read_csv_auto(base() || '/icu/outputevents.csv.gz', header = true) o
  JOIN cohort_v2 c USING (stay_id)
  WHERE o.itemid IN (SELECT itemid FROM urine_items)
    AND o.charttime BETWEEN c.t0 - INTERVAL 6 HOUR AND c.t0
  GROUP BY o.stay_id
),
i AS (
  SELECT x.stay_id, sum(TRY_CAST(x.amount AS DOUBLE)) AS input_6h, count(*) AS n_input_rows
  FROM read_csv_auto(base() || '/icu/inputevents.csv.gz', header = true) x
  JOIN cohort_v2 c USING (stay_id)
  WHERE lower(x.amountuom) LIKE '%ml%'
    AND x.starttime BETWEEN c.t0 - INTERVAL 6 HOUR AND c.t0
  GROUP BY x.stay_id
)
SELECT c.stay_id,
       u.urine_6h, CASE WHEN u.urine_6h IS NULL THEN 1 ELSE 0 END AS urine_missing,
       i.input_6h, CASE WHEN i.input_6h IS NULL THEN 1 ELSE 0 END AS input_missing,
       i.input_6h - u.urine_6h AS fluid_balance_6h
FROM cohort_v2 c
LEFT JOIN u USING (stay_id)
LEFT JOIN i USING (stay_id);

-- ----------------------------------------------------------------------------
-- STEP 6  撤药时所在 ICU（原用 first_careunit，中途转科者会错）
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE unit_at_t0 AS
SELECT stay_id, careunit FROM (
  SELECT c.stay_id, t.careunit,
         row_number() OVER (PARTITION BY c.stay_id ORDER BY t.intime DESC) AS rn
  FROM cohort_v2 c
  JOIN read_csv_auto(base() || '/hosp/transfers.csv*', header = true) t
    ON CAST(t.hadm_id AS BIGINT) = c.hadm_id
  WHERE t.intime <= c.t0
    AND t.careunit IS NOT NULL AND t.careunit <> ''
) WHERE rn = 1;

-- ----------------------------------------------------------------------------
-- STEP 7  合并导出
-- ----------------------------------------------------------------------------
CREATE OR REPLACE TABLE analytic_v2 AS
SELECT l.stay_id, l.subject_id, l.hadm_id, l.intime, l.t0, l.t_landmark,
       p.gender,
       p.anchor_age + (year(l.t0) - p.anchor_year) AS age,
       l.admission_type, l.insurance,
       coalesce(un.careunit, l.first_careunit) AS careunit_t0,
       l.early_wean_6h, l.early_wean_2h, l.early_wean_1h, l.stop6h, l.restart_24h,
       l.dose_t0, l.map_t0, l.hours_obs_6h,
       bd.map_mean_6h, bd.map_min_6h, bd.frac_map_lt65,
       l.hours_ne_before_t0, bd.hours_icu_before_t0,
       lb.lactate, lb.creatinine, lb.bilirubin, lb.platelet, lb.wbc, lb.hemoglobin, lb.bun,
       s.gcs_min, s.on_vent, s.on_rrt, s.on_vasopressin,
       f.urine_6h, f.urine_missing, f.input_6h, f.input_missing, f.fluid_balance_6h,
       -- 结局：landmark 起算的 28 天全因死亡
       CASE WHEN l.death_time IS NOT NULL
                 AND l.death_time <= l.t_landmark + INTERVAL 28 DAY THEN 1 ELSE 0 END AS event_28d,
       CASE WHEN l.death_time IS NOT NULL
                 AND l.death_time <= l.t_landmark + INTERVAL 28 DAY
            THEN date_diff('second', l.t_landmark, l.death_time) / 86400.0
            ELSE 28.0 END AS time_days,
       CASE WHEN l.death_time IS NOT NULL
            THEN date_diff('second', l.t_landmark, l.death_time) / 86400.0 END AS days_to_death,
       l.death_time, l.dischtime, l.outtime
FROM landmark l
LEFT JOIN (
  SELECT subject_id, gender, anchor_age, anchor_year
  FROM read_csv_auto(base() || '/hosp/patients.csv.gz', header = true)
) p ON p.subject_id = l.subject_id
LEFT JOIN labs_v2 lb USING (stay_id)
LEFT JOIN fluids_v2 f USING (stay_id)
LEFT JOIN unit_at_t0 un USING (stay_id)
LEFT JOIN support s USING (stay_id)
LEFT JOIN burden bd USING (stay_id);

COPY analytic_v2 TO 'C:/Users/caofu/tte/analytic_v2.csv' (HEADER, DELIMITER ',');

-- ----------------------------------------------------------------------------
-- STEP 8  核对输出
-- ----------------------------------------------------------------------------
SELECT '--- 流程（对照 flow chart）---' AS section;

SELECT '02 的 t0 候选（每 stay 首次）' AS step, count(*) AS n FROM cohort_base
UNION ALL SELECT '首次 stay 的 t0', count(*) FROM cohort_v2
UNION ALL SELECT 'landmark 内死亡/转出剔除后', count(*) FROM landmark
UNION ALL SELECT '最终入分析', count(*) FROM analytic_v2;

SELECT '--- 暴露分布（三套窗口）---' AS section;
SELECT round(100.0*avg(early_wean_6h),1) AS pct_6h,
       round(100.0*avg(early_wean_2h),1) AS pct_2h,
       round(100.0*avg(early_wean_1h),1) AS pct_1h,
       round(100.0*avg(stop6h),1) AS pct_stop6h,
       count(*) AS n, sum(event_28d) AS deaths_28d,
       round(100.0*avg(event_28d),1) AS pct_death
FROM analytic_v2;

SELECT '--- 粗率（方向核对，非结论）---' AS section;
SELECT early_wean_6h, count(*) AS n,
       round(100.0*avg(event_28d),1) AS pct_death,
       round(100.0*avg(restart_24h),1) AS pct_restart,
       round(avg(dose_t0),3) AS mean_dose_t0,
       round(avg(time_days),1) AS mean_time_days
FROM analytic_v2 GROUP BY early_wean_6h ORDER BY early_wean_6h;

SELECT '--- 缺失率 ---' AS section;
SELECT 'lactate' AS var, round(100.0*avg(CASE WHEN lactate IS NULL THEN 1 ELSE 0 END),1) AS missing_pct FROM analytic_v2
UNION ALL SELECT 'bilirubin', round(100.0*avg(CASE WHEN bilirubin IS NULL THEN 1 ELSE 0 END),1) FROM analytic_v2
UNION ALL SELECT 'gcs_min',   round(100.0*avg(CASE WHEN gcs_min   IS NULL THEN 1 ELSE 0 END),1) FROM analytic_v2
UNION ALL SELECT 'urine_6h',  round(100.0*avg(urine_missing),1) FROM analytic_v2
UNION ALL SELECT 'input_6h',  round(100.0*avg(input_missing),1) FROM analytic_v2
UNION ALL SELECT 'creatinine',round(100.0*avg(CASE WHEN creatinine IS NULL THEN 1 ELSE 0 END),1) FROM analytic_v2;

SELECT '--- 各 ICU 撤药实践（工具变量候选）---' AS section;
SELECT careunit_t0, count(*) AS n, round(100.0*avg(early_wean_6h),1) AS pct_early,
       round(100.0*avg(event_28d),1) AS pct_death
FROM analytic_v2 GROUP BY careunit_t0 HAVING count(*) >= 50 ORDER BY pct_early DESC;

SELECT '--- 敏感性：1 小时窗下的阳性性 ---' AS section;
SELECT count(*) AS n, sum(early_wean_1h) AS n_early,
       round(100.0*avg(early_wean_1h),1) AS pct_early
FROM analytic_v2;
