#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""构建时校验：稿件里的关键数字必须与 results/key_numbers.json 一致。
   任何一处对不上就退出码 1，防止把不一致的稿子投出去。
"""
import os, re, json, subprocess, sys

BASE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(BASE, "results")
K = json.load(open(os.path.join(R, "key_numbers.json"), encoding="utf-8"))
docx = os.path.join(BASE, "IEM_submission.docx")

txt = subprocess.run(["pandoc", "-t", "plain", docx], capture_output=True, text=True).stdout
flat = re.sub(r"\s+", " ", txt)
errors, checked = [], 0

import math
def r2(x):
    return math.floor(float(x) * 100 + 0.5) / 100

def need(s, label):
    global checked
    checked += 1
    if s not in flat:
        errors.append(f"[缺失] {label}: 找不到「{s}」")

M, F, E = K["main"], K["flow"], K.get("eicu", {})
MH = K.get("mimic_hosp", {})

# 1. 关键结论必须出现
need(f"Of {M['n']:,} patients with sepsis", "样本量")
need(f"{M['early']:,} ({M['pct_early']}%) attempted weaning", "暴露人数")
need(f"{M['deaths']:,} ({M['pct_death']}%) died within 28 days", "死亡人数")
need(f"hazard ratio {M['hr']:.2f}, 95% confidence interval {M['lo']:.2f}-{M['hi']:.2f}", "主 HR（摘要）")
need(f"{M['hr']:.2f} (95% CI {M['lo']:.2f}–{M['hi']:.2f}", "主 HR（结果）")
need(f"E-value was {M['evalue']}", "E-value（摘要）")
need(f"E-value for the primary hazard ratio was {M['evalue']}", "E-value（结果）")
need(f"effective sample size {M['ess']:,} of {M['n']:,}", "有效样本量")
need(f"c-statistic {M['cstat']:.2f}", "c 统计量")
need(f"p={M['ph_p']:.2f} for the exposure", "PH 检验")
need(f"{M['restart_early']}% of the early group and {M['restart_delayed']}%", "重启结局")
need(f"{M['surv_early']}% versus {M['surv_delayed']}%", "28 天生存率")

# 2. 流程
need(f"{F['mimic_ne']:,} ICU stays receiving norepinephrine", "MIMIC 起始")
need(f"{F['mimic_eligible']:,} had at least one eligible hour and {F['mimic_sepsis']:,}", "MIMIC 流程")
need(f"left {F['mimic_after_filters']:,}", "MIMIC 过滤后")
need(f"yielded {F['mimic_analytic']:,} patients", "MIMIC 最终")

# 3. 外部验证
if E:
    need(f"({E['n']} patients, {E['events']} in-hospital deaths)", "eICU 样本")
    need(f"hazard ratio {E['hr']:.2f}, 95% confidence interval {E['lo']:.2f}-{E['hi']:.2f}", "eICU HR")
    need(f"p={E['p']:.2f}", "eICU p 值")
if MH:
    need(f"{r2(MH['hr']):.2f} ({r2(MH['lo']):.2f}–{r2(MH['hi']):.2f}", "MIMIC 院内死亡 HR")

# 4. 表 1 表头
T1 = K["table1"]
need(f"(n={T1['n_early']})", "表 1 早撤组 n")
need(f"(n={T1['n_delayed']})", "表 1 延迟组 n")

# 5. 所有 HR 型数字必须在允许集合内
allowed = set()
def add(hr, lo=None, hi=None):
    allowed.add(r2(hr))
for src in (M, MH, E):
    if src and "hr" in src:
        add(src["hr"])
add(M["unadj_hr"]); add(M["dr_hr"])
for r in K["sensitivity"]:
    add(r["HR"])
for r in K["subgroups"]:
    add(r["HR"])

bad = set()
for m in re.finditer(r"(0\.\d{2}) \(0\.\d{2}–0\.\d{2}", flat):
    ctx = flat[max(0, m.start() - 45):m.start()].lower()
    if not re.search(r"(?<![a-z])(hr|hazard)", ctx):
        continue          # 例如剂量 0.06 (0.04–0.08)，不是风险比
    v = float(m.group(1))
    if v not in allowed:
        bad.add(v)
if bad:
    errors.append(f"[不一致] 正文出现未在数据源中的 HR：{sorted(bad)}")


if E and "aps_med" in E:
    need(f"Median APACHE IV acute physiology score was {int(E['aps_med'])} (IQR {int(E['aps_q1'])}–{int(E['aps_q3'])})", "APACHE 分数")
    need(f"predicted hospital mortality {E['phm_med']}% (IQR {E['phm_q1']}–{E['phm_q3']})", "APACHE 预测死亡率")
    need(f"observed in-hospital mortality of {E['pct_death']}%", "观察死亡率")
need(f"was {M['max_smd']:.2f}", "加权后最大 SMD")
need(f"(1.{int(round(M['evalue_ci']*100))%100:02d} for the confidence limit)", "E-value 置信限")


# 6. 结构断言：小节必须待在正确的章节里
def pos(x):
    return flat.find(x)
im, ir, isec, idis = pos("Methods"), pos("Results"), pos("Secondary outcome"), pos("Discussion")
if not (im < ir < isec < idis):
    errors.append(f"[结构] 章节顺序异常: Methods@{im} Results@{ir} Secondary outcome@{isec} Discussion@{idis}")
checked += 1

print(f"检查项 {checked} 项，允许的 HR 值 {sorted(allowed)}")
if errors:
    print("\n发现 %d 处问题：" % len(errors))
    for e in errors:
        print(" -", e)
    sys.exit(1)
print("\n全部一致：稿件数字与 key_numbers.json 吻合。")
