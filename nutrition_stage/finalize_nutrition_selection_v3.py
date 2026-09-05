from __future__ import annotations

import finalize_nutrition_selection_v2 as base

base.base.CANDIDATES["C-N17"]={
 "status":base.base.ROOT/"output_c_n17"/"status.json",
 "title_cn":"维生素D不足成人中低膳食钙摄入与继发性甲状旁腺激素升高：NHANES 2003–2006的时间复制",
 "title_en":"Low dietary calcium intake and secondary parathyroid hormone elevation among adults with vitamin D insufficiency: temporal replication in NHANES 2003–2006",
 "query":"(NHANES OR National Health and Nutrition Examination Survey) AND (dietary calcium OR calcium intake) AND (vitamin D insufficiency OR vitamin D deficiency) AND (parathyroid hormone OR PTH OR secondary hyperparathyroidism)",
 "tokens":[["nhanes","national health and nutrition examination survey"],["calcium"],["vitamin d"],["parathyroid hormone"," pth ","secondary hyperparathyroidism"]],
}

_old_metric=base.base.metric_line
def metric_line(code,s):
 if code!="C-N17":return _old_metric(code,s)
 p=s.get("prep") or {};q=s.get("primary") or {};a=q.get("adjusted_difference") or {};b=q.get("adjusted_pr") or {}
 return [
  f"实际纳入：{p.get('actual_n')}；完全调整：{p.get('complete_case_n')}（保留率{p.get('retention')}）",
  f"PTH升高事件：{p.get('events')}；低钙摄入组：{p.get('key_group_n')}，组内事件：{p.get('key_group_events')}",
  f"调整后PTH差：{a.get('estimate')} pg/mL（95%CI {a.get('lower')}至{a.get('upper')}）",
  f"调整后PTH升高患病比：{b.get('estimate')}（95%CI {b.get('lower')}至{b.get('upper')}）",
  f"敏感性同向：{s.get('sensitivity_concordant')}/{s.get('sensitivity_total')}",
 ]
base.base.metric_line=metric_line

if __name__=="__main__":base.base.main()
