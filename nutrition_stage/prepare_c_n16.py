from __future__ import annotations
import io,json,re
from pathlib import Path
import numpy as np,pandas as pd,requests
BASE="https://raw.githubusercontent.com/protobi/nhanes-continuous/main";OUT=Path("nutrition_stage/output_c_n16");OUT.mkdir(parents=True,exist_ok=True)
CYCLES={"D":1,"E":2,"F":3,"G":4,"H":5,"I":6,"J":7}
PPI=re.compile(r"omeprazole|pantoprazole|esomeprazole|lansoprazole|rabeprazole|dexlansoprazole",re.I);DIURETIC=re.compile(r"hydrochlorothiazide|chlorthalidone|furosemide|bumetanide|torsemide|indapamide|metolazone",re.I);MAGRX=re.compile(r"magnesium oxide|magnesium chloride|magnesium gluconate|magnesium lactate",re.I)
def fetch(path,required=True):
 r=requests.get(f"{BASE}/{path}",timeout=120)
 if r.status_code!=200:
  if required:raise RuntimeError(f"{path}:HTTP{r.status_code}")
  return pd.DataFrame(columns=["SEQN"])
 return pd.read_csv(io.BytesIO(r.content),low_memory=False)
def keep(x,cols):
 use=[c for c in cols if c in x]
 if "SEQN" in x and "SEQN" not in use:use.insert(0,"SEQN")
 return x[use].copy()
def num(x,names):
 o=pd.Series(np.nan,index=x.index,dtype=float)
 for n in names:
  if n in x:o=o.fillna(pd.to_numeric(x[n],errors="coerce"))
 return o
def egfr(scr,age,sex):
 f=sex.eq(2);k=np.where(f,.7,.9);a=np.where(f,-.241,-.302);r=scr/k
 return 142*np.minimum(r,1)**a*np.maximum(r,1)**-1.2*.9938**age*np.where(f,1.012,1)
def meds(rx):
 if rx.empty:return pd.DataFrame(columns=["SEQN","ppi","diuretic","magnesium_rx","ppi_days"])
 cs=[c for c in ["RXDDRUG","RXDINGFL","RXDDRGID"] if c in rx];t=rx[cs].fillna("").astype(str).agg(" ".join,axis=1) if cs else pd.Series("",index=rx.index);z=pd.DataFrame({"SEQN":rx.SEQN,"ppi":t.str.contains(PPI),"diuretic":t.str.contains(DIURETIC),"magnesium_rx":t.str.contains(MAGRX),"ppi_days":np.where(t.str.contains(PPI),num(rx,["RXDDAYS"]),np.nan)});return z.groupby("SEQN",as_index=False).agg(ppi=("ppi","max"),diuretic=("diuretic","max"),magnesium_rx=("magnesium_rx","max"),ppi_days=("ppi_days","max"))
def cycle(cy):
 idx=CYCLES[cy];audit=[]
 def load(path,req=True):
  x=fetch(path,req);audit.append({"cycle":cy,"path":path,"required":req,"rows":len(x),"columns":x.shape[1]});return x
 demo=load(f"Demographics/DEMO_{cy}.csv");bmx=load(f"Examination/BMX_{cy}.csv");d1=load(f"Dietary/DR1TOT_{cy}.csv");d2=load(f"Dietary/DR2TOT_{cy}.csv",False);bio=load(f"Laboratory/BIOPRO_{cy}.csv");rx=load(f"Questionnaire/RXQ_RX_{cy}.csv");smq=load(f"Questionnaire/SMQ_{cy}.csv",False);diq=load(f"Questionnaire/DIQ_{cy}.csv",False);mcq=load(f"Questionnaire/MCQ_{cy}.csv",False)
 m=keep(demo,["SEQN","RIAGENDR","RIDAGEYR","RIDRETH1","RIDRETH3","WTMEC2YR","SDMVPSU","SDMVSTRA","INDFMPIR"])
 for p in [keep(bmx,["SEQN","BMXBMI"]),keep(d1,["SEQN","WTDRD1","DR1TKCAL","DR1TMAGN","DR1TCALC","DR1TPOTA"]),keep(d2,["SEQN","DR2TKCAL","DR2TMAGN"]),keep(bio,["SEQN","LBXSCR","LBXSMGSI","LBXSMSI","LBXSMG"]),meds(rx),keep(smq,["SEQN","SMQ020","SMQ040"]),keep(diq,["SEQN","DIQ010"]),keep(mcq,["SEQN","MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F"])]:m=m.merge(p.drop_duplicates("SEQN"),on="SEQN",how="left")
 m["cycle"]=cy;m["period"]="discovery" if cy in "DEFG" else "validation";m["age"]=num(m,["RIDAGEYR"]);m["sex"]=num(m,["RIAGENDR"]);m["race"]=num(m,["RIDRETH1","RIDRETH3"]);m["bmi"]=num(m,["BMXBMI"]);m["pir"]=num(m,["INDFMPIR"]);m["weight"]=num(m,["WTDRD1"])/7;m["psu_u"]=idx*100+num(m,["SDMVPSU"]);m["strata_u"]=idx*1000+num(m,["SDMVSTRA"]);m["scr"]=num(m,["LBXSCR"]);m["egfr"]=egfr(m.scr,m.age,m.sex)
 raw=num(m,["LBXSMG","LBXSMGSI","LBXSMSI"]);med=raw.median();m["serum_mg_mg_dl"]=raw*2.4305 if pd.notna(med) and med<1.5 else raw;m["energy1"]=num(m,["DR1TKCAL"]);m["magnesium1"]=num(m,["DR1TMAGN"]);m["calcium1"]=num(m,["DR1TCALC"]);m["potassium1"]=num(m,["DR1TPOTA"]);m["magnesium_mean2"]=pd.concat([m.magnesium1,num(m,["DR2TMAGN"])],axis=1).mean(axis=1,skipna=False)
 m["smoking"]=np.select([num(m,["SMQ020"]).eq(2),num(m,["SMQ020"]).eq(1)&num(m,["SMQ040"]).eq(3),num(m,["SMQ020"]).eq(1)&num(m,["SMQ040"]).isin([1,2])],["never","former","current"],default=np.nan);m["diabetes"]=np.where(num(m,["DIQ010"]).eq(1),1,np.where(num(m,["DIQ010"]).isin([2,3]),0,np.nan));cv=pd.concat([num(m,[v]).eq(1) for v in ["MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F"]],axis=1);m["cvd"]=cv.any(axis=1).astype(int)
 for v in ["ppi","diuretic","magnesium_rx"]:m[v]=m.get(v,False);m[v]=m[v].fillna(False).astype(bool)
 return m,audit
def main():
 fs=[];aud=[]
 for cy in CYCLES:
  x,a=cycle(cy);fs.append(x);aud.extend(a)
 raw=pd.concat(fs,ignore_index=True,sort=False);flow=[]
 def rec(s,x):flow.append({"step":s,"n":len(x)})
 rec("all",raw);c=raw[raw.age.ge(40)&raw.ppi].copy();rec("age_40_plus_current_ppi",c);c=c[c.energy1.between(500,5000)&c.magnesium1.between(20,2000)&c.serum_mg_mg_dl.between(1.0,4.0)].copy();rec("valid_diet_and_serum_magnesium",c);c=c[c.weight.gt(0)&c.psu_u.notna()&c.strata_u.notna()].copy();rec("valid_survey_design",c)
 c["low_mg_intake"]=np.where(c.sex.eq(2),c.magnesium1<265,c.magnesium1<350).astype(int);c["hypomagnesemia"]=(c.serum_mg_mg_dl<1.7).astype(int);cov=["age","sex","race","bmi","egfr","smoking","diabetes","cvd","energy1","calcium1","potassium1"];c["complete_case"]=c[cov].notna().all(axis=1);rec("complete_case_main",c[c.complete_case]);cols=["SEQN","cycle","period","age","sex","race","weight","psu_u","strata_u","pir","bmi","egfr","smoking","diabetes","cvd","energy1","magnesium1","magnesium_mean2","calcium1","potassium1","serum_mg_mg_dl","low_mg_intake","hypomagnesemia","diuretic","magnesium_rx","ppi_days","complete_case"];c[cols].to_csv(OUT/"core.csv",index=False);pd.DataFrame(flow).to_csv(OUT/"flow.csv",index=False);pd.DataFrame(aud).to_csv(OUT/"source_audit.csv",index=False)
 miss=[]
 for scope,x in [("overall",c),("discovery",c[c.period=="discovery"]),("validation",c[c.period=="validation"])]:
  for v in cov+["pir","magnesium_mean2","ppi_days"]:miss.append({"scope":scope,"variable":v,"n":len(x),"missing_n":int(x[v].isna().sum()),"missing_pct":float(x[v].isna().mean()*100)})
 pd.DataFrame(miss).to_csv(OUT/"missingness.csv",index=False);gr=[]
 for scope,x in [("overall",c),("discovery",c[c.period=="discovery"]),("validation",c[c.period=="validation"])]:
  for z,g in x.groupby("low_mg_intake"):gr.append({"scope":scope,"low_intake":int(z),"n":len(g),"events":int(g.hypomagnesemia.sum()),"event_pct":float(g.hypomagnesemia.mean()*100),"mean_serum_mg":float(g.serum_mg_mg_dl.mean())})
 pd.DataFrame(gr).to_csv(OUT/"group_counts.csv",index=False);pd.DataFrame([{"concept":"dietary_magnesium","field":"DR1TMAGN","unit":"mg/day","semantics":"day-1 24-hour dietary recall","rule":"below sex-specific EAR: <265 women, <350 men"},{"concept":"current_PPI","field":"RXDDRUG/ingredient text","unit":"prescription","semantics":"current named PPI ingredient","rule":"frozen six-drug regex"},{"concept":"serum_magnesium","field":"cycle biochemistry magnesium","unit":"mg/dL","semantics":"measured serum magnesium; SI converted if required","rule":"hypomagnesemia <1.7; plausible 1-4"},{"concept":"survey_design","field":"WTDRD1/PSU/strata","unit":"survey","semantics":"day-1 dietary design","rule":"weight/7; cycle-unique IDs"}]).to_csv(OUT/"semantic_audit.csv",index=False)
 prep={"candidate_code":"C-N16","actual_n":len(c),"complete_case_n":int(c.complete_case.sum()),"retention":float(c.complete_case.mean()),"events":int(c.hypomagnesemia.sum()),"key_group_n":int(c.low_mg_intake.sum()),"key_group_events":int(c.loc[c.low_mg_intake.eq(1),"hypomagnesemia"].sum()),"validation_key_group_n":int(c.loc[(c.period=="validation")&c.low_mg_intake.eq(1)].shape[0]),"validation_key_group_events":int(c.loc[(c.period=="validation")&c.low_mg_intake.eq(1),"hypomagnesemia"].sum())};(OUT/"prep_status.json").write_text(json.dumps(prep,indent=2));print(json.dumps(prep,indent=2))
if __name__=="__main__":main()
