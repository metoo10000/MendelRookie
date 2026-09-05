from __future__ import annotations

import finalize_nutrition_selection_v4 as base

root=base.base.base.base
root.CANDIDATES["C-N19"]={
 "status":root.ROOT/"output_c_n19"/"status.json",
 "title_cn":"已限制钠摄入的治疗中高血压成人中钾摄入不足与血压控制失败：NHANES 2005–2018的时间验证",
 "title_en":"Potassium inadequacy and uncontrolled blood pressure among treated hypertensive adults already meeting sodium-intake targets: temporal validation in NHANES 2005–2018",
 "query":"(NHANES OR National Health and Nutrition Examination Survey) AND (potassium intake OR dietary potassium) AND (sodium intake OR sodium restriction) AND (treated hypertension OR antihypertensive) AND (blood pressure control OR uncontrolled hypertension)",
 "tokens":[["nhanes","national health and nutrition examination survey"],["potassium"],["sodium"],["hypertension","blood pressure"],["treat","antihypertensive"]],
}
_old=root.metric_line
def metric_line(code,s):
 if code!="C-N19":return _old(code,s)
 p=s.get("prep") or {};q=s.get("primary") or {};a=q.get("adjusted_sbp_difference") or {};b=q.get("adjusted_pr") or {}
 return [
  f"实际纳入：{p.get('actual_n')}；完全调整：{p.get('complete_case_n')}（保留率{p.get('retention')}）",
  f"血压未控制事件：{p.get('events')}；钾摄入不足组：{p.get('key_group_n')}，组内事件：{p.get('key_group_events')}；参照组：{p.get('reference_group_n')}，事件：{p.get('reference_group_events')}",
  f"调整后收缩压差：{a.get('estimate')} mmHg（95%CI {a.get('lower')}至{a.get('upper')}）",
  f"调整后未控制患病比：{b.get('estimate')}（95%CI {b.get('lower')}至{b.get('upper')}）",
  f"敏感性同向：{s.get('sensitivity_concordant')}/{s.get('sensitivity_total')}",
 ]
root.metric_line=metric_line
if __name__=="__main__":root.main()
