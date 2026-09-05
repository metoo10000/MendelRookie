from __future__ import annotations

import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import requests

BASE="https://raw.githubusercontent.com/protobi/nhanes-continuous/main"
TREE_URL="https://api.github.com/repos/protobi/nhanes-continuous/git/trees/main?recursive=1"
OUT=Path("nutrition_stage/output_c_n18");OUT.mkdir(parents=True,exist_ok=True)
CYCLES={"B":("2001-2002",1,"discovery"),"C":("2003-2004",2,"discovery"),"D":("2005-2006",3,"validation")}

def tree_paths():
 r=requests.get(TREE_URL,timeout=120);r.raise_for_status();return [x["path"] for x in r.json()["tree"] if x.get("type")=="blob"]
def read(path):
 r=requests.get(f"{BASE}/{path}",timeout=120);r.raise_for_status();return pd.read_csv(io.BytesIO(r.content),low_memory=False)
def exact(tree,paths,optional=False):
 for p in paths:
  if p not in tree:continue
  try:return read(p),p
  except Exception:pass
 if optional:return pd.DataFrame(columns=["SEQN"]),""
 raise RuntimeError(str(paths))
def find_field(tree,cy,target):
 specs={
  "folate":({"LBXFOL","LBDFOLSI","LBXFOLSI"},re.compile(r"FOL|L06|LAB06",re.I)),
  "hcy":({"LBXHCY","LBDHCY"},re.compile(r"HCY|HOMO|L06|LAB06",re.I)),
  "b12":({"LBXB12","LBDVITB12","LBDVB12"},re.compile(r"B12|COBAL|L06|LAB06",re.I)),
 }
 wanted,rx=specs[target];candidates=[p for p in tree if p.startswith("Laboratory/") and p.endswith(f"_{cy}.csv") and rx.search(p)]
 exacts=[f"Laboratory/L06_{cy}.csv",f"Laboratory/LAB06_{cy}.csv",f"Laboratory/FOLATE_{cy}.csv",f"Laboratory/HCY_{cy}.csv",f"Laboratory/B12_{cy}.csv"]
 for p in dict.fromkeys(exacts+candidates):
  if p not in tree:continue
  try:x=read(p)
  except Exception:continue
  cs=[c for c in x.columns if c.upper() in wanted]
  if not cs:
   patt={"folate":r"FOL","hcy":r"HCY|HOMO","b12":r"B12|COBAL"}[target]
   cs=[c for c in x.columns if re.search(patt,c.upper()) and c.upper().startswith(("LBX","LBD"))]
  for c in cs:
   q=pd.to_numeric(x[c],errors="coerce");m=q.median()
   if target=="folate" and pd.notna(m) and 1<m<100:return x,p,c
   if target=="hcy" and pd.notna(m) and 2<m<50:return x,p,c
   if target=="b12" and pd.notna(m) and 50<m<2000:return x,p,c
 raise RuntimeError(f"{target} unavailable {cy}")
def keep(x,cols):
 u=[c for c in cols if c in x]
 if "SEQN" in x and "SEQN" not in u:u.insert(0,"SEQN")
 return x[u].copy()
def num(x,names):
 o=pd.Series(np.nan,index=x.index,dtype=float)
 for n in names:
  if n in x:o=o.fillna(pd.to_numeric(x[n],errors="coerce"))
 return o
def egfr(scr,age,sex):
 f=sex.eq(2);k=np.where(f,.7,.9);a=np.where(f,-.241,-.302);r=scr/k
 return 142*np.minimum(r,1)**a*np.maximum(r,1)**-1.2*.9938**age*np.where(f,1.012,1)
def merge(tree,cy):
 years,idx,period=CYCLES[cy];audit=[]
 def load(paths,opt=False):
  x,p=exact(tree,paths,opt);audit.append({"cycle":cy,"years":years,"used":p,"optional":opt,"rows":len(x),"columns":x.shape[1]});return x
 demo=load([f"Demographics/DEMO_{cy}.csv"]);bmx=load([f"Examination/BMX_{cy}.csv"]);bio=load([f"Laboratory/BIOPRO_{cy}.csv",f"Laboratory/L40_{cy}.csv"]);smq=load([f"Questionnaire/SMQ_{cy}.csv"]);diq=load([f"Questionnaire/DIQ_{cy}.csv"],True);alq=load([f"Questionnaire/ALQ_{cy}.csv"],True)
 fol,fp,fc=find_field(tree,cy,"folate");hcy,hp,hc=find_field(tree,cy,"hcy");b12,bp,bc=find_field(tree,cy,"b12");audit += [{"cycle":cy,"used":fp,"concept":"folate","field":fc,"rows":len(fol)},{"cycle":cy,"used":hp,"concept":"homocysteine","field":hc,"rows":len(hcy)},{"cycle":cy,"used":bp,"concept":"B12","field":bc,"rows":len(b12)}]
 m=keep(demo,["SEQN","RIAGENDR","RIDAGEYR","RIDRETH1","RIDRETH3","WTMEC2YR","SDMVPSU","SDMVSTRA","INDFMPIR"])
 for p in [keep(bmx,["SEQN","BMXBMI"]),keep(bio,["SEQN","LBXSCR"]),keep(smq,["SEQN","SMQ020","SMQ040"]),keep(diq,["SEQN","DIQ010"]),keep(alq,["SEQN"]+[c for c in alq.columns if c.startswith("ALQ")]),keep(fol,["SEQN",fc]),keep(hcy,["SEQN",hc]),keep(b12,["SEQN",bc])]:m=m.merge(p.drop_duplicates("SEQN"),on="SEQN",how="left")
 m["cycle"]=cy;m["period"]=period;m["age"]=num(m,["RIDAGEYR"]);m["sex"]=num(m,["RIAGENDR"]);m["race"]=num(m,["RIDRETH1","RIDRETH3"]);m["bmi"]=num(m,["BMXBMI"]);m["pir"]=num(m,["INDFMPIR"]);m["weight"]=num(m,["WTMEC2YR"])/3;m["psu_u"]=idx*100+num(m,["SDMVPSU"]);m["strata_u"]=idx*1000+num(m,["SDMVSTRA"]);m["scr"]=num(m,["LBXSCR"]);m["egfr"]=egfr(m.scr,m.age,m.sex)
 f=num(m,[fc]);fm=f.median();m["folate_ng_ml"]=f/2.266 if pd.notna(fm) and fm>30 else f;m["hcy_umol_l"]=num(m,[hc]);m["b12_pg_ml"]=num(m,[bc]);m["smoking"]=np.select([num(m,["SMQ020"]).eq(2),num(m,["SMQ020"]).eq(1)&num(m,["SMQ040"]).eq(3),num(m,["SMQ020"]).eq(1)&num(m,["SMQ040"]).isin([1,2])],["never","former","current"],default=np.nan);m["diabetes"]=np.where(num(m,["DIQ010"]).eq(1),1,np.where(num(m,["DIQ010"]).isin([2,3]),0,np.nan));alc_cols=[c for c in m if c.startswith("ALQ")];m["alcohol_marker"]=num(m,alc_cols[:4]) if alc_cols else np.nan
 return m,audit
def main():
 tree=tree_paths();fs=[];aud=[]
 for cy in CYCLES:
  x,a=merge(tree,cy);fs.append(x);aud.extend(a)
 raw=pd.concat(fs,ignore_index=True,sort=False);flow=[]
 def rec(s,x):flow.append({"step":s,"n":len(x)})
 rec("all_linked",raw);c=raw[raw.age.between(20,79)&raw.smoking.eq("current")].copy();rec("current_smokers_age_20_79",c);c=c[c.folate_ng_ml.between(3,50)&c.hcy_umol_l.between(3,100)&c.b12_pg_ml.between(80,3000)&c.egfr.ge(60)].copy();rec("nondeficient_folate_valid_hcy_b12_no_ckd",c);c=c[c.weight.gt(0)&c.psu_u.notna()&c.strata_u.notna()].copy();rec("valid_survey_design",c);c["low_normal_folate"]=(c.folate_ng_ml<6).astype(int);c["hyperhcy"]=(c.hcy_umol_l>15).astype(int);cov=["age","sex","race","bmi","egfr","b12_pg_ml","diabetes"];c["complete_case"]=c[cov].notna().all(axis=1);rec("complete_case_main",c[c.complete_case]);cols=["SEQN","cycle","period","age","sex","race","weight","psu_u","strata_u","pir","bmi","egfr","b12_pg_ml","diabetes","alcohol_marker","folate_ng_ml","low_normal_folate","hcy_umol_l","hyperhcy","complete_case"];c[cols].to_csv(OUT/"core.csv",index=False);pd.DataFrame(flow).to_csv(OUT/"flow.csv",index=False);pd.DataFrame(aud).to_csv(OUT/"source_audit.csv",index=False)
 miss=[]
 for scope,x in [("overall",c),("discovery",c[c.period=="discovery"]),("validation",c[c.period=="validation"])]:
  for v in cov+["pir","alcohol_marker"]:miss.append({"scope":scope,"variable":v,"n":len(x),"missing_n":int(x[v].isna().sum()),"missing_pct":float(x[v].isna().mean()*100)})
 pd.DataFrame(miss).to_csv(OUT/"missingness.csv",index=False);gr=[]
 for scope,x in [("overall",c),("discovery",c[c.period=="discovery"]),("validation",c[c.period=="validation"])]:
  for z,g in x.groupby("low_normal_folate"):gr.append({"scope":scope,"low_normal":int(z),"n":len(g),"events":int(g.hyperhcy.sum()),"event_pct":float(g.hyperhcy.mean()*100),"mean_hcy":float(g.hcy_umol_l.mean())})
 pd.DataFrame(gr).to_csv(OUT/"group_counts.csv",index=False);pd.DataFrame([{"concept":"low_normal_folate","field":"cycle serum folate","unit":"ng/mL","semantics":"measured serum folate; frank deficiency excluded","rule":"3 to <6 vs >=6"},{"concept":"homocysteine","field":"LBXHCY-equivalent","unit":"umol/L","semantics":"measured total homocysteine","rule":"hyperhomocysteinemia >15; plausible 3-100"},{"concept":"current_smoking","field":"SMQ020/SMQ040","unit":"questionnaire","semantics":">=100 lifetime and now every/some days","rule":"main population"},{"concept":"survey_design","field":"WTMEC2YR/PSU/strata","unit":"survey","semantics":"MEC design","rule":"weight/3; cycle unique IDs"}]).to_csv(OUT/"semantic_audit.csv",index=False)
 prep={"candidate_code":"C-N18","actual_n":len(c),"complete_case_n":int(c.complete_case.sum()),"retention":float(c.complete_case.mean()),"events":int(c.hyperhcy.sum()),"key_group_n":int(c.low_normal_folate.sum()),"key_group_events":int(c.loc[c.low_normal_folate.eq(1),"hyperhcy"].sum()),"validation_key_group_n":int(c.loc[(c.period=="validation")&c.low_normal_folate.eq(1)].shape[0]),"validation_key_group_events":int(c.loc[(c.period=="validation")&c.low_normal_folate.eq(1),"hyperhcy"].sum())};(OUT/"prep_status.json").write_text(json.dumps(prep,indent=2));print(json.dumps(prep,indent=2))
if __name__=="__main__":main()
