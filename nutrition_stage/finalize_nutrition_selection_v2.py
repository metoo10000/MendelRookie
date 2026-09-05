from __future__ import annotations

import finalize_nutrition_selection as base

base.CANDIDATES["C-N16"]={
 "status":base.ROOT/"output_c_n16"/"status.json",
 "title_cn":"质子泵抑制剂使用者中膳食镁摄入不足与低镁血症：NHANES 2005–2018 的时间验证研究",
 "title_en":"Inadequate dietary magnesium intake and hypomagnesemia among proton-pump inhibitor users: temporal validation in NHANES 2005–2018",
 "query":"(NHANES OR National Health and Nutrition Examination Survey) AND (dietary magnesium OR magnesium intake) AND (proton pump inhibitor OR PPI) AND (hypomagnesemia OR serum magnesium)",
 "tokens":[["nhanes","national health and nutrition examination survey"],["magnesium"],["proton pump inhibitor"," ppi "],["hypomagnesemia","serum magnesium"]],
}

_old_metric=base.metric_line
def metric_line(code,s):
 if code!="C-N16":return _old_metric(code,s)
 p=s.get("prep") or {};q=s.get("primary") or {};a=q.get("adjusted_difference") or {};b=q.get("adjusted_pr") or {}
 return [
  f"实际纳入：{p.get('actual_n')}；完全调整：{p.get('complete_case_n')}（保留率{p.get('retention')}）",
  f"低镁血症事件：{p.get('events')}；关键组：{p.get('key_group_n')}，关键组事件：{p.get('key_group_events')}",
  f"调整后血清镁差：{a.get('estimate')} mg/dL（95%CI {a.get('lower')}至{a.get('upper')}）",
  f"调整后患病比：{b.get('estimate')}（95%CI {b.get('lower')}至{b.get('upper')}）",
  f"敏感性同向：{s.get('sensitivity_concordant')}/{s.get('sensitivity_total')}",
 ]
base.metric_line=metric_line

if __name__=="__main__":base.main()
