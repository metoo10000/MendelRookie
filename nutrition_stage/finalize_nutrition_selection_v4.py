from __future__ import annotations

import finalize_nutrition_selection_v3 as base

base.base.base.CANDIDATES["C-N18"]={
 "status":base.base.base.ROOT/"output_c_n18"/"status.json",
 "title_cn":"当前吸烟者中非缺乏范围内低叶酸状态与高同型半胱氨酸血症：NHANES 2001–2006的时间复制",
 "title_en":"Low-normal folate status and hyperhomocysteinemia among current smokers: temporal replication in NHANES 2001–2006",
 "query":"(NHANES OR National Health and Nutrition Examination Survey) AND (serum folate OR folate status) AND (current smoker OR smoking) AND (homocysteine OR hyperhomocysteinemia)",
 "tokens":[["nhanes","national health and nutrition examination survey"],["folate"],["smok","tobacco"],["homocysteine","hyperhomocysteinemia"]],
}

_old_metric=base.base.base.metric_line
def metric_line(code,s):
 if code!="C-N18":return _old_metric(code,s)
 p=s.get("prep") or {};q=s.get("primary") or {};a=q.get("adjusted_difference") or {};b=q.get("adjusted_pr") or {}
 return [
  f"实际纳入：{p.get('actual_n')}；完全调整：{p.get('complete_case_n')}（保留率{p.get('retention')}）",
  f"高同型半胱氨酸事件：{p.get('events')}；低正常叶酸组：{p.get('key_group_n')}，组内事件：{p.get('key_group_events')}",
  f"调整后同型半胱氨酸差：{a.get('estimate')} µmol/L（95%CI {a.get('lower')}至{a.get('upper')}）",
  f"调整后患病比：{b.get('estimate')}（95%CI {b.get('lower')}至{b.get('upper')}）",
  f"敏感性同向：{s.get('sensitivity_concordant')}/{s.get('sensitivity_total')}",
 ]
base.base.base.metric_line=metric_line

if __name__=="__main__":base.base.base.main()
