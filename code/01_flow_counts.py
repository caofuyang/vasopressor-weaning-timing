#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regenerate the cohort flow counts used by Figure 1 and key_numbers.json.

Requires the DuckDB databases produced by sql/03_cohort_v2.sql (MIMIC) and
sql/07_eicu_cohort.sql (eICU). Writes results/flow_mimic.csv and
results/flow_eicu.csv.

Usage:
    python code/01_flow_counts.py --mimic tte.duckdb --eicu eicu.duckdb
"""
import argparse, os
import duckdb, pandas as pd

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
os.makedirs(OUT, exist_ok=True)


def mimic_flow(db):
    con = duckdb.connect(db, read_only=True)
    df = con.sql("""
        select
          (select count(*) from stay_win)                                        as ne_stays,
          (select count(*) from cohort_base)                                     as eligible_hour,
          (select count(*) from cohort_base cb join sepsis_hadm s using (hadm_id)) as sepsis,
          (select count(*) from cohort_v2)                                       as after_filters,
          (select count(*) from landmark)                                        as landmark,
          (select count(*) from analytic_v2 where early_wean_6h = 1)             as early,
          (select count(*) from analytic_v2 where early_wean_6h = 0)             as delayed
    """).df()
    df.to_csv(os.path.join(OUT, "flow_mimic.csv"), index=False)
    return df


def eicu_flow(db):
    con = duckdb.connect(db, read_only=True)
    df = con.sql("""
        select
          (select count(distinct stay_id) from ne_raw)      as ne_stays,
          (select count(distinct stay_id) from ne_dose)     as convertible,
          (select count(distinct stay_id) from t0_cases)    as eligible,
          (select count(*) from cohort_base)                as sepsis_first,
          (select count(*) from cohort_v2)                  as after_ne6h,
          (select count(*) from landmark)                   as landmark,
          (select count(*) from analytic where early_wean_6h = 1) as early,
          (select count(*) from analytic where early_wean_6h = 0) as delayed
    """).df()
    df.to_csv(os.path.join(OUT, "flow_eicu.csv"), index=False)
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mimic", required=True)
    ap.add_argument("--eicu", required=True)
    a = ap.parse_args()
    print("MIMIC:\n", mimic_flow(a.mimic).to_string(index=False))
    print("\neICU:\n", eicu_flow(a.eicu).to_string(index=False))
