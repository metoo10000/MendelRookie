from __future__ import annotations
import io,json
from pathlib import Path
import numpy as np,pandas as pd,requests
BASE="https://raw.githubusercontent.com/protobi/nhanes-continuous/main";OUT=Path("nutrition_stage/output_c_n19");OUT.mkdir(parents=True,exist_ok=True)
CYCLES={"D":1,"E":2,"F":3,"G":4,"H":5,"I":6,"J":7}
def fetch(path,required=True):
 r=requests.get(f"{BASE}/{path}",timeout=120)
 if r.status_code!=200:
  if required:raise RuntimeError(f"{path}:HTTP{r.status_code}")
  return pd.DataFrame(columns=["SEQN"])
 return pd.read_csv(io.BytesIO(r.content),low_memory=False)
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
def bpmean(m,prefix):
 cols=[c for c in [f"{prefix}{i}" for i in [2,3,4,1]] if c in m]
 vals=m[cols].apply(pd.to_numeric,errors="coerce")
 # Mean of all available valid readings, requiring at least two.
 vals=vals.where(vals.between(30 if "DI" in prefix else 60,150 if "DI" in prefix else 260))
 return vals.mean(axis=1).where(vals.notna().sum(axis=1)>=2)
def cycle(cy):
 idx=CYCLES[cy];audit=[]
 def load(path,opt=False):
  x=fetch(path,not opt);audit.append({"cycle":cy,"path":path,"optional":opt,"rows":len(x),"columns":x.shape[1]});return x
 demo=load(f"Demographics/DEMO_{cy}.csv");bmx=load(f"Examination/BMX_{cy}.csv");bpx=load(f"Examination/BPX_{cy}.csv");bpq=load(f"Questionnaire/BPQ_{cy}.csv");d1=load(f"Dietary/DR1TOT_{cy}.csv");d2=load(f"Dietary/DR2TOT_{cy}.csv",True);bio=load(f"Laboratory/BIOPRO_{cy}.csv");smq=load(f"Questionnaire/SMQ_{cy}.csv",True);diq=load(f"Questionnaire/DIQ_{cy}.csv",True);mcq=load(f"Questionnaire/MCQ_{cy}.csv",True)
 m=keep(demo,["SEQN","RIAGENDR","RIDAGEYR","RIDRETH1","RIDRETH3","SDMVPSU","SDMVSTRA","INDFMPIR"])
 for p in [keep(bmx,["SEQN","BMXBMI"]),bpx,keep(bpq,["SEQN","BPQ020","BPQ040A","BPQ050A"]),keep(d1,["SEQN","WTDRD1","DR1TKCAL","DR1TPOTA","DR1TSODI","DR1TMAGN","DR1TFIBE"]),keep(d2,["SEQN","DR2TKCAL","DR2TPOTA","DR2TSODI"]),keep(bio,["SEQN","LBXSCR"]),keep(smq,["SEQN","SMQ020","SMQ040"]),keep(diq,["SEQN","DIQ010"]),keep(mcq,["SEQN","MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F"])]:m=m.merge(p.drop_duplicates("SEQN"),on="SEQN",how="left")
 m["cycle"]=cy;m["period"]="discovery" if cy in "DEFG" else "validation";m["age"]=num(m,["RIDAGEYR"]);m["sex"]=num(m,["RIAGENDR"]);m["race"]=num(m,["RIDRETH1","RIDRETH3"]);m["bmi"]=num(m,["BMXBMI"]);m["pir"]=num(m,["INDFMPIR"]);m["weight"]=num(m,["WTDRD1"])/7;m["psu_u"]=idx*100+num(m,["SDMVPSU"]);m["strata_u"]=idx*1000+num(m,["SDMVSTRA"]);m["sbp"]=bpmean(m,"BPXSY");m["dbp"]=bpmean(m,"BPXDI");m["treated_htn"]=(num(m,["BPQ020"]).eq(1)&num(m,["BPQ050A","BPQ040A"]).eq(1));m["energy1"]=num(m,["DR1TKCAL"]);m["potassium1"]=num(m,["DR1TPOTA"]);m["sodium1"]=num(m,["DR1TSODI"]);m["magnesium1"]=num(m,["DR1TMAGN"]);m["fiber1"]=num(m,["DR1TFIBE"]);m["potassium2"]=num(m,["DR2TPOTA"]);m["sodium2"]=num(m,["DR2TSODI"]);m["potassium_mean2"]=pd.concat([m.potassium1,m.potassium2],axis=1).mean(axis=1,skipna=False);m["sodium_mean2"]=pd.concat([m.sodium1,m.sodium2],axis=1).mean(axis=1,skipna=False);m["scr"]=num(m,["LBXSCR"]);m["egfr"]=egfr(m.scr,m.age,m.sex);m["smoking"]=np.select([num(m,["SMQ020"]).eq(2),num(m,["SMQ020"]).eq(1)&num(m,["SMQ040"]).eq(3),num(m,["SMQ020"]).eq(1)&num(m,["SMQ040"]).isin([1,2])],["never","former","current"],default=np.nan);m["diabetes"]=np.where(num(m,["DIQ010"]).eq(1),1,np.where(num(m,["DIQ010"]).isin([2,3]),0,np.nan));cv=pd.concat([num(m,[v]).eq(1) for v in ["MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F"]],axis=1);m["cvd"]=cv.any(axis=1).astype(int)
 return m,audit
def main():
 fs=[];aud=[]
 for cy in CYCLES:
  x,a=cycle(cy);fs.append(x);aud.extend(a)
 raw=pd.concat(fs,ignore_index=True,sort=False);flow=[]
 def rec(s,x):flow.append({"step":s,"n":len(x)})
 rec("all",raw);c=raw[raw.age.ge(30)&raw.treated_htn].copy();rec("treated_hypertension_age_30_plus",c);c=c[c.energy1.between(500,5000)&c.potassium1.between(100,10000)&c.sodium1.between(100,15000)&c.sodium1.lt(2300)].copy();rec("valid_day1_diet_and_sodium_below_2300",c);c=c[c.sbp.notna()&c.dbp.notna()&c.weight.gt(0)&c.psu_u.notna()&c.strata_u.notna()].copy();rec("valid_objective_bp_and_survey_design",c);c["potassium_inadequate"]=np.where(c.sex.eq(2),c.potassium1<2600,c.potassium1<3400).astype(int);c["uncontrolled_140_90"]=((c.sbp>=140)|(c.dbp>=90)).astype(int);c["uncontrolled_130_80"]=((c.sbp>=130)|(c.dbp>=80)).astype(int);c["k_density"]=c.potassium1/c.energy1*1000;cov=["age","sex","race","bmi","egfr","smoking","diabetes","cvd","energy1","sodium1"];c["complete_case"]=c[cov].notna().all(axis=1);rec("complete_case_main",c[c.complete_case]);cols=["SEQN","cycle","period","age","sex","race","weight","psu_u","strata_u","pir","bmi","egfr","smoking","diabetes","cvd","energy1","potassium1","potassium_mean2","sodium1","sodium_mean2","magnesium1","fiber1","k_density","sbp","dbp","potassium_inadequate","uncontrolled_140_90","uncontrolled_130_80","complete_case"];c[cols].to_csv(OUT/"core.csv",index=False);pd.DataFrame(flow).to_csv(OUT/"flow.csv",index=False);pd.DataFrame(aud).to_csv(OUT/"source_audit.csv",index=False)
 miss=[]
 for scope,x in [("overall",c),("discovery",c[c.period=="discovery"]),("validation",c[c.period=="validation"])]:
  for v in cov+["pir","potassium_mean2","sodium_mean2","magnesium1","fiber1"]:miss.append({"scope":scope,"variable":v,"n":len(x),"missing_n":int(x[v].isna().sum()),"missing_pct":float(x[v].isna().mean()*100)})
 pd.DataFrame(miss).to_csv(OUT/"missingness.csv",index=False);gr=[]
 for scope,x in [("overall",c),("discovery",c[c.period=="discovery"]),("validation",c[c.period=="validation"])]:
  for z,g in x.groupby("potassium_inadequate"):gr.append({"scope":scope,"inadequate":int(z),"n":len(g),"events":int(g.uncontrolled_140_90.sum()),"event_pct":float(g.uncontrolled_140_90.mean()*100),"mean_sbp":float(g.sbp.mean())})
 pd.DataFrame(gr).to_csv(OUT/"group_counts.csv",index=False);pd.DataFrame([{"concept":"treated_hypertension","field":"BPQ020 plus BPQ050A/BPQ040A","unit":"questionnaire","semantics":"self-reported diagnosed hypertension and current prescribed treatment","rule":"main population"},{"concept":"guideline_concordant_sodium","field":"DR1TSODI","unit":"mg/day","semantics":"day-1 24-hour recall","rule":"<2300 mg/day"},{"concept":"potassium_inadequacy","field":"DR1TPOTA","unit":"mg/day","semantics":"day-1 24-hour recall","rule":"<2600 women; <3400 men"},{"concept":"objective_BP","field":"BPXSY/BPXDI repeated readings","unit":"mmHg","semantics":"mean of at least two valid examination readings","rule":"uncontrolled >=140 systolic or >=90 diastolic"},{"concept":"survey_design","field":"WTDRD1/PSU/strata","unit":"survey","semantics":"day-1 dietary design","rule":"weight/7; cycle unique IDs"}]).to_csv(OUT/"semantic_audit.csv",index=False)
 ref=int((c.potassium_inadequate==0).sum());refe=int(c.loc[c.potassium_inadequate.eq(0),"uncontrolled_140_90"].sum());prep={"candidate_code":"C-N19","actual_n":len(c),"complete_case_n":int(c.complete_case.sum()),"retention":float(c.complete_case.mean()),"events":int(c.uncontrolled_140_90.sum()),"key_group_n":int(c.potassium_inadequate.sum()),"key_group_events":int(c.loc[c.potassium_inadequate.eq(1),"uncontrolled_140_90"].sum()),"reference_group_n":ref,"reference_group_events":refe,"validation_reference_n":int(c.loc[(c.period=="validation")&c.potassium_inadequate.eq(0)].shape[0]),"validation_reference_events":int(c.loc[(c.period=="validation")&c.potassium_inadequate.eq(0),"uncontrolled_140_90"].sum())};(OUT/"prep_status.json").write_text(json.dumps(prep,indent=2));print(json.dumps(prep,indent=2))
if __name__=="__main__":main()
