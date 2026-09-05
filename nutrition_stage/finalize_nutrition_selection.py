from __future__ import annotations

import json
import re
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any
from urllib.parse import quote

import requests

ROOT=Path("nutrition_stage")
OUT=ROOT/"final_selection";OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"nutrition-topic-audit/1.0 (research use)"}

CANDIDATES={
 "C-N11":{
  "status":ROOT/"output_batch_n11"/"batch_status.json",
  "title_cn":"使用排钾利尿剂的成人中低膳食钾摄入与低钾血症：NHANES 2005–2018 的时间验证研究",
  "title_en":"Low dietary potassium intake and hypokalemia among adults using potassium-wasting diuretics: temporal validation in NHANES 2005–2018",
  "query":"(NHANES OR National Health and Nutrition Examination Survey) AND (dietary potassium OR potassium intake) AND (diuretic OR thiazide OR loop diuretic) AND (hypokalemia OR serum potassium)",
  "tokens":[["nhanes","national health and nutrition examination survey"],["potassium"],["diuretic","thiazide","furosemide","loop"],["hypokalemia","serum potassium"]]},
 "C-N12":{
  "status":ROOT/"output_batch_n11"/"batch_status.json",
  "title_cn":"噻嗪类利尿剂使用者中低膳食钠摄入与低钠血症：NHANES 2005–2018 的时间验证研究",
  "title_en":"Low dietary sodium intake and hyponatremia among thiazide users: temporal validation in NHANES 2005–2018",
  "query":"(NHANES OR National Health and Nutrition Examination Survey) AND (dietary sodium OR sodium intake) AND (thiazide) AND (hyponatremia OR serum sodium)",
  "tokens":[["nhanes","national health and nutrition examination survey"],["sodium"],["thiazide"],["hyponatremia","serum sodium"]]},
 "C-N13":{
  "status":ROOT/"output_batch_n11"/"batch_status.json",
  "title_cn":"ACEI或ARB使用者中高膳食钾摄入与高钾血症：NHANES 2005–2018 的时间验证研究",
  "title_en":"High dietary potassium intake and hyperkalemia among ACE inhibitor or ARB users: temporal validation in NHANES 2005–2018",
  "query":"(NHANES OR National Health and Nutrition Examination Survey) AND (dietary potassium OR potassium intake) AND (ACE inhibitor OR angiotensin receptor blocker OR ARB) AND (hyperkalemia OR serum potassium)",
  "tokens":[["nhanes","national health and nutrition examination survey"],["potassium"],["ace inhibitor","angiotensin receptor blocker"," arb "],["hyperkalemia","serum potassium"]]},
 "C-N14":{
  "status":ROOT/"output_c_n14"/"status.json",
  "title_cn":"美国成人血清锌缺乏与听力学测定的高频听力损失：NHANES 2011–2012与2015–2016的时间复制",
  "title_en":"Serum zinc deficiency and audiometrically measured high-frequency hearing loss in US adults: temporal replication across NHANES 2011–2012 and 2015–2016",
  "query":"(NHANES OR National Health and Nutrition Examination Survey) AND (serum zinc OR zinc deficiency) AND (hearing loss OR audiometry OR pure tone)",
  "tokens":[["nhanes","national health and nutrition examination survey"],["zinc"],["hearing loss","audiometry","pure tone","hearing threshold"]]},
 "C-N15":{
  "status":ROOT/"output_c_n15"/"status.json",
  "title_cn":"美国成人血清维生素C缺乏与功能牙列丧失：跨NHANES调查周期的时间复制",
  "title_en":"Serum vitamin C deficiency and loss of functional dentition in US adults: temporal replication across NHANES survey cycles",
  "query":"(NHANES OR National Health and Nutrition Examination Survey) AND (vitamin C OR ascorbic acid) AND (tooth loss OR dentition OR edentulism OR periodontal)",
  "tokens":[["nhanes","national health and nutrition examination survey"],["vitamin c","ascorbic acid","ascorbate"],["tooth loss","dentition","edentulism","periodontal"]]},
}


def load_candidate_status(code:str)->dict[str,Any]|None:
 p=CANDIDATES[code]["status"]
 if not p.exists():return None
 raw=json.loads(p.read_text())
 if code in {"C-N11","C-N12","C-N13"}:
  return (raw.get("candidates") or {}).get(code)
 return raw

def epmc_search(query:str)->list[dict[str,Any]]:
 url="https://www.ebi.ac.uk/europepmc/webservices/rest/search"
 params={"query":query,"format":"json","pageSize":100,"resultType":"core"}
 r=requests.get(url,params=params,headers=UA,timeout=60);r.raise_for_status();return r.json().get("resultList",{}).get("result",[])
def pmc_fulltext(pmcid:str)->tuple[str,str]:
 try:
  r=requests.get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML",headers=UA,timeout=60)
  if r.status_code!=200:return "",f"HTTP{r.status_code}"
  root=ET.fromstring(r.content);text=" ".join(" ".join(root.itertext()).split());return text,"PMC full text"
 except Exception as e:return "",type(e).__name__
def doi_fulltext(doi:str)->tuple[str,str]:
 if not doi:return "","no DOI"
 try:
  r=requests.get(f"https://doi.org/{doi}",headers=UA,timeout=60,allow_redirects=True)
  if r.status_code!=200:return "",f"HTTP{r.status_code}"
  text=re.sub(r"<script.*?</script>|<style.*?</style>"," ",r.text,flags=re.I|re.S)
  text=re.sub(r"<[^>]+>"," ",text);text=" ".join(text.split())
  return (text if len(text)>5000 else ""),("publisher full text" if len(text)>5000 else "publisher page insufficient")
 except Exception as e:return "",type(e).__name__
def token_match(text:str,groups:list[list[str]])->bool:
 t=" "+text.lower()+" ";return all(any(tok.lower() in t for tok in group) for group in groups)
def audit_candidate(code:str)->dict[str,Any]:
 cfg=CANDIDATES[code];hits=epmc_search(cfg["query"]);rows=[];exact=[];unresolved=[]
 for h in hits:
  title=h.get("title") or "";abstract=h.get("abstractText") or "";year=int(str(h.get("pubYear") or "0")[:4] or 0);pmcid=h.get("pmcid") or "";doi=h.get("doi") or ""
  preliminary=token_match(title+" "+abstract,cfg["tokens"])
  full="";source="not required"
  if preliminary:
   if pmcid:full,source=pmc_fulltext(pmcid)
   if not full and doi:full,source=doi_fulltext(doi)
  exact_overlap=bool(full and token_match(full,cfg["tokens"]))
  recent_exact=exact_overlap and year>=2024
  if recent_exact:exact.append(title)
  if preliminary and not full:unresolved.append(title)
  rows.append({"title":title,"year":year,"journal":h.get("journalTitle"),"pmid":h.get("pmid"),"pmcid":pmcid,"doi":doi,"preliminary_relevance":preliminary,"full_text_source":source,"full_text_characters":len(full),"exact_core_overlap":exact_overlap,"recent_exact_overlap":recent_exact})
  time.sleep(.05)
 # Conservative gate: any recent exact overlap fails; any unresolved likely exact paper also blocks passage.
 pass_gate=(len(exact)==0 and len(unresolved)==0)
 return {"candidate_code":code,"query":cfg["query"],"records_screened":len(rows),"full_text_relevant_records":sum(bool(r["full_text_characters"]) for r in rows),"recent_exact_overlap_count":len(exact),"recent_exact_titles":exact,"unresolved_high_relevance_count":len(unresolved),"unresolved_titles":unresolved,"dedup_gate_pass":pass_gate,"records":rows}

def metric_line(code:str,s:dict[str,Any])->list[str]:
 if code in {"C-N11","C-N12","C-N13"}:
  return [f"实际完全调整样本：{s.get('actual_n')}",f"主要事件：{s.get('events')}；关键组：{s.get('key_group_n')}，关键组事件：{s.get('key_group_events')}",f"调整后连续结局差：{(s.get('adjusted_difference') or {}).get('estimate')}（95%CI {(s.get('adjusted_difference') or {}).get('lower')}至{(s.get('adjusted_difference') or {}).get('upper')}）",f"调整后患病比：{(s.get('adjusted_pr') or {}).get('estimate')}（95%CI {(s.get('adjusted_pr') or {}).get('lower')}至{(s.get('adjusted_pr') or {}).get('upper')}）",f"敏感性同向：{s.get('sensitivity_concordant')}/{s.get('sensitivity_total')}"]
 if code=="C-N14":
  p=s.get("prep") or {};q=s.get("primary") or {};a=q.get("adjusted_pta_difference") or {};b=q.get("adjusted_pr") or {}
  return [f"实际纳入：{p.get('actual_n')}；完全调整：{p.get('complete_case_n')}（保留率{p.get('retention')}）",f"听力损失事件：{p.get('events')}；关键组：{p.get('key_group_n')}，关键组事件：{p.get('key_group_events')}",f"调整后高频听阈差：{a.get('estimate')} dB（95%CI {a.get('lower')}至{a.get('upper')}）",f"调整后患病比：{b.get('estimate')}（95%CI {b.get('lower')}至{b.get('upper')}）",f"敏感性同向：{s.get('sensitivity_concordant')}/{s.get('sensitivity_total')}"]
 p=s.get("prep") or {};q=s.get("primary") or {};a=q.get("adjusted_teeth_difference") or {};b=q.get("adjusted_pr") or {}
 return [f"实际纳入：{p.get('actual_n')}；完全调整：{p.get('complete_case_n')}（保留率{p.get('retention')}）",f"功能牙列丧失事件：{p.get('events')}；关键组：{p.get('key_group_n')}，关键组事件：{p.get('key_group_events')}",f"调整后牙数差：{a.get('estimate')}颗（95%CI {a.get('lower')}至{a.get('upper')}）",f"调整后患病比：{b.get('estimate')}（95%CI {b.get('lower')}至{b.get('upper')}）",f"敏感性同向：{s.get('sensitivity_concordant')}/{s.get('sensitivity_total')}"]

def main():
 data_pass=[];audits={};statuses={}
 for code in CANDIDATES:
  s=load_candidate_status(code);statuses[code]=s
  if s and s.get("data_gate_pass") is True:data_pass.append(code)
 selected=None
 for code in data_pass:
  audit=audit_candidate(code);audits[code]=audit
  (OUT/f"{code}_dedup.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2))
  if audit["dedup_gate_pass"]:selected=code;break
 final={"status_files_found":[k for k,v in statuses.items() if v is not None],"data_pass_candidates":data_pass,"dedup_audits":audits,"selected_candidate":selected,"final_state":"QUALIFIED_TOPIC_FOUND" if selected else "NO_QUALIFIED_TOPIC_YET"}
 (OUT/"final_status.json").write_text(json.dumps(final,ensure_ascii=False,indent=2))
 lines=["# 营养选题最终验收",""]
 if selected:
  cfg=CANDIDATES[selected];s=statuses[selected];a=audits[selected]
  lines += ["## 最终判定：GO","",f"### {cfg['title_cn']}","",f"**English title:** {cfg['title_en']}","",f"候选编号：{selected}",""]
  lines += [f"- {x}" for x in metric_line(selected,s)]
  lines += [f"- 数据门禁：全部通过。",f"- 全文排重：检索记录{a['records_screened']}篇；近期同核心全文重复{a['recent_exact_overlap_count']}篇；未解决高相关全文{a['unresolved_high_relevance_count']}篇；排重门禁通过。",""]
  lines += ["### 主张边界","","这是观察性、横断面或调查周期复制结果，只支持关联与风险识别；不解释为营养暴露的因果治疗效应。第二数据库未被强行加入。"]
 else:
  lines += ["## 尚无候选同时通过数据与全文排重门禁","",f"数据门禁通过候选：{', '.join(data_pass) if data_pass else '0'}。未释放任何题名。"]
 (OUT/"final_report.md").write_text("\n".join(lines),encoding="utf-8")
 print(json.dumps(final,ensure_ascii=False,indent=2))
if __name__=="__main__":main()
