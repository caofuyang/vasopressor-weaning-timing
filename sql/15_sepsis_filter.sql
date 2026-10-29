-- ============================================================================
-- 第 0 步：脓毒症筛选的影响评估（在重跑队列之前先看会掉多少样本）
-- 运行：& E:\数据库\duckdb\duckdb.exe "C:/Users/caofu/tte/tte.duckdb" -c ".read 'C:/Users/caofu/tte/15_sepsis_filter.sql'"
-- 说明：本脚本只做计数，不改动已有表。确认影响后再改 03 / 07 并重跑。
-- ============================================================================

CREATE OR REPLACE MACRO base() AS 'E:/数据库/MIMIC-IV/data/mimic-iv-3.1';
.timer on

-- ---------------------------------------------------------------------------
-- MIMIC：ICD 诊断中的脓毒症
--   ICD-9 : 038.x（败血症）、995.91/995.92（脓毒症/严重脓毒症）、785.52（脓毒性休克）
--   ICD-10: A40.x/A41.x（脓毒症）、R65.2x（严重脓毒症/脓毒性休克）
-- ---------------------------------------------------------------------------
CREATE OR REPLACE TABLE sepsis_hadm AS
SELECT DISTINCT hadm_id
FROM read_csv_auto(base() || '/hosp/diagnoses_icd.csv.gz', header = true)
WHERE (TRY_CAST(icd_version AS INT) = 9
       AND (icd_code LIKE '038%' OR icd_code IN ('99591','99592','78552')))
   OR (TRY_CAST(icd_version AS INT) = 10
       AND (icd_code LIKE 'A40%' OR icd_code LIKE 'A41%' OR icd_code LIKE 'R652%'));

SELECT '--- 脓毒症筛选对 MIMIC 队列的影响 ---' AS section;

SELECT '当前分析队列' AS step, count(*) AS n FROM analytic_v2
UNION ALL
SELECT '其中有脓毒症 ICD 诊断', count(*) FROM analytic_v2 a
  JOIN sepsis_hadm s USING (hadm_id)
UNION ALL
SELECT '当前队列中无脓毒症诊断', count(*) FROM analytic_v2 a
  ANTI JOIN sepsis_hadm s USING (hadm_id);

SELECT '--- 按暴露分组（筛选后）---' AS section;
SELECT a.early_wean_6h, count(*) AS n,
       round(100.0*avg(a.event_28d),1) AS pct_death_28d
FROM analytic_v2 a JOIN sepsis_hadm s USING (hadm_id)
GROUP BY 1 ORDER BY 1;

-- ---------------------------------------------------------------------------
-- eICU：diagnosis 表
-- ---------------------------------------------------------------------------
CREATE OR REPLACE MACRO ebase() AS 'E:/数据库/eICU-CRD/data/eicu-collaborative-research-database-2.0';

CREATE OR REPLACE TABLE sepsis_stay AS
SELECT DISTINCT patientunitstayid AS stay_id
FROM read_csv_auto(ebase() || '/diagnosis.csv.gz', header = true)
WHERE icd9code LIKE '038%'
   OR lower(diagnosisstring) LIKE '%sepsis%'
   OR lower(diagnosisstring) LIKE '%septic%';

SELECT '--- 脓毒症筛选对 eICU 队列的影响 ---' AS section;
SELECT '当前 eICU 队列' AS step, count(*) AS n FROM analytic
UNION ALL
SELECT '其中有脓毒症诊断', count(*) FROM analytic a JOIN sepsis_stay s USING (stay_id)
UNION ALL
SELECT '其中无脓毒症诊断', count(*) FROM analytic a ANTI JOIN sepsis_stay s USING (stay_id);

-- ---------------------------------------------------------------------------
-- 若影响可接受（预期保留 60–75%），下一步：
--   1. 在 03_cohort_v2.sql 的 cohort_v2 定义处加：
--        JOIN sepsis_hadm USING (hadm_id)
--   2. 在 07_eicu_cohort.sql 的 cohort_base 定义处加：
--        JOIN sepsis_stay USING (stay_id)
--   3. 依次重跑 03 → 04 → 05 → 12 → 13 → 10 → 11 → build_docx_iem.js
-- ---------------------------------------------------------------------------
