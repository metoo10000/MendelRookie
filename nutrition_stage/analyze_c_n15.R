options(stringsAsFactors = FALSE)
options(survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(jsonlite)
})

out_dir <- "nutrition_stage/output_c_n15"
dat <- read.csv(file.path(out_dir, "c_n15_analysis_core.csv"), na.strings = c("", "NA"))
num_vars <- c("age","sex","race","education","pir","weight","psu_u","strata_u","mean_kcal","mean_protein","mean_potassium","mean_sodium","mean_magnesium","day1_sodium","day2_sodium","sodium_density","uacr","albuminuria","severe_albuminuria","log_uacr","log_sodium","creatinine","egfr","albumin","bmi","body_weight","height","diabetes","hypertension","cvd","heart_failure","sbp","hba1c")
for (v in intersect(num_vars, names(dat))) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
for (v in c("thiazide","loop","any_diuretic","raas","ppi","complete_case")) dat[[v]] <- as.logical(dat[[v]])
dat$any_diuretic <- as.numeric(dat$any_diuretic)
dat$cycle <- factor(dat$cycle, levels = c("E","F","G","H","I","J"))
dat$period <- factor(dat$period, levels = c("discovery","validation"))
dat$sex <- factor(dat$sex, levels = c(1,2), labels = c("male","female"))
dat$race <- factor(dat$race)
dat$smoking <- factor(dat$smoking, levels = c("never","former","current"))
dat$exposure <- factor(dat$exposure, levels = c("lower","high"))

make_design <- function(d) svydesign(ids=~psu_u, strata=~strata_u, weights=~weight, nest=TRUE, data=d)

safe_effect <- function(model, term, exponentiate=FALSE) {
  cf <- coef(model)
  if (!(term %in% names(cf))) return(rep(NA_real_,4))
  se <- sqrt(diag(vcov(model)))[term]
  est <- unname(cf[term])
  p <- tryCatch(unname(summary(model)$coefficients[term,ncol(summary(model)$coefficients)]), error=function(e) NA_real_)
  z <- c(est, est-1.96*se, est+1.96*se, p)
  if (exponentiate) z[1:3] <- exp(z[1:3])
  z
}

clean_vars <- function(d, vars) {
  out <- vars
  for (v in intersect(c("sex","race","smoking","cycle"), out)) {
    if (length(unique(d[[v]][!is.na(d[[v]])])) < 2) out <- setdiff(out,v)
  }
  out
}

model_rows <- list()
append_model <- function(scope,name,model,term="exposurehigh",exponentiate=FALSE) {
  z <- safe_effect(model,term,exponentiate)
  model_rows[[length(model_rows)+1]] <<- data.frame(candidate_code="C-N15",scope=scope,model=name,term=term,estimate=z[1],lower=z[2],upper=z[3],p=z[4],n=nrow(model$model),stringsAsFactors=FALSE)
  z
}

fit_scope <- function(d,scope) {
  d <- droplevels(d[d$complete_case & !is.na(d$exposure),,drop=FALSE])
  if (nrow(d)<250 || sum(d$exposure=="high")<50 || sum(d$albuminuria==1)<30 || nlevels(d$exposure)<2) return(NULL)
  des <- make_design(d)
  crude_vars <- clean_vars(d,c("exposure","cycle"))
  full_vars <- clean_vars(d,c("exposure","age","sex","race","bmi","egfr","diabetes","hypertension","cvd","sbp","hba1c","mean_kcal","mean_protein","mean_potassium","any_diuretic","smoking","pir","cycle"))
  f <- function(y,x) as.formula(paste(y,"~",paste(x,collapse=" + ")))
  crude_pr <- svyglm(f("albuminuria",crude_vars),design=des,family=quasipoisson(link="log"))
  crude_rd <- svyglm(f("albuminuria",crude_vars),design=des,family=gaussian())
  crude_log <- svyglm(f("log_uacr",crude_vars),design=des,family=gaussian())
  adj_pr <- svyglm(f("albuminuria",full_vars),design=des,family=quasipoisson(link="log"))
  adj_rd <- svyglm(f("albuminuria",full_vars),design=des,family=gaussian())
  adj_log <- svyglm(f("log_uacr",full_vars),design=des,family=gaussian())
  severe_pr <- svyglm(f("severe_albuminuria",full_vars),design=des,family=quasipoisson(link="log"))
  append_model(scope,"crude_albuminuria_pr",crude_pr,exponentiate=TRUE)
  append_model(scope,"crude_albuminuria_rd",crude_rd)
  append_model(scope,"crude_uacr_geometric_ratio",crude_log,exponentiate=TRUE)
  pr <- append_model(scope,"adjusted_albuminuria_pr",adj_pr,exponentiate=TRUE)
  rd <- append_model(scope,"adjusted_albuminuria_rd",adj_rd)
  gmr <- append_model(scope,"adjusted_uacr_geometric_ratio",adj_log,exponentiate=TRUE)
  spr <- append_model(scope,"adjusted_severe_albuminuria_pr",severe_pr,exponentiate=TRUE)
  list(data=d,pr=pr,rd=rd,gmr=gmr,severe_pr=spr)
}

primary <- fit_scope(dat,"overall")
discovery <- fit_scope(dat[dat$period=="discovery",],"discovery")
validation <- fit_scope(dat[dat$period=="validation",],"validation")
if (is.null(primary)) stop("Primary C-N15 model could not be estimated")

# Missingness mechanism.
miss <- droplevels(dat[!is.na(dat$exposure)&!is.na(dat$age)&!is.na(dat$sex)&!is.na(dat$race),,drop=FALSE])
mv <- clean_vars(miss,c("exposure","age","sex","race","cycle"))
miss_mod <- svyglm(as.formula(paste("I(complete_case) ~",paste(mv,collapse=" + "))),design=make_design(miss),family=quasibinomial())
mz <- safe_effect(miss_mod,"exposurehigh",TRUE)
missingness_df <- data.frame(term="exposurehigh",odds_ratio=mz[1],lower=mz[2],upper=mz[3],p=mz[4])

# Restricted cubic spline.
weighted_quantile <- function(x,w,p){ok<-is.finite(x)&is.finite(w)&w>0;x<-x[ok];w<-w[ok];o<-order(x);x<-x[o];w<-w[o];cw<-cumsum(w)/sum(w);sapply(p,function(q)x[which(cw>=q)[1]])}
rcs_basis <- function(x,k){p3<-function(z)pmax(z,0)^3;den<-(k[length(k)]-k[1])^2;h<-sapply(seq_len(length(k)-2),function(j)(p3(x-k[j])-p3(x-k[length(k)-1])*(k[length(k)]-k[j])/(k[length(k)]-k[length(k)-1])+p3(x-k[length(k)])*(k[length(k)-1]-k[j])/(k[length(k)]-k[length(k)-1]))/den);colnames(h)<-paste0("rcs",seq_len(ncol(h)));h}
nl <- primary$data
knots <- weighted_quantile(nl$log_sodium,nl$weight,c(.05,.35,.65,.95))
b <- rcs_basis(nl$log_sodium,knots);nl$xlin<-nl$log_sodium;nl$rcs1<-b[,1];nl$rcs2<-b[,2]
nv <- clean_vars(nl,c("xlin","rcs1","rcs2","age","sex","race","bmi","egfr","diabetes","hypertension","cvd","sbp","hba1c","mean_kcal","mean_protein","mean_potassium","any_diuretic","smoking","pir","cycle"))
fn <- function(y) as.formula(paste(y,"~",paste(nv,collapse=" + ")))
dn <- make_design(nl)
nl_pr <- svyglm(fn("albuminuria"),design=dn,family=quasipoisson(link="log"))
nl_log <- svyglm(fn("log_uacr"),design=dn,family=gaussian())
nonlinear_df <- rbind(
  data.frame(outcome="albuminuria",n=nrow(nl),p_nonlinear=unname(regTermTest(nl_pr,~rcs1+rcs2)$p),knots_mg_day=paste(round(exp(knots)),collapse="|")),
  data.frame(outcome="log_uacr",n=nrow(nl),p_nonlinear=unname(regTermTest(nl_log,~rcs1+rcs2)$p),knots_mg_day=paste(round(exp(knots)),collapse="|"))
)

# Frozen sensitivity analyses.
sens_rows <- list()
run_sens <- function(label,d,low_cut=2300,high_cut=3000,outcome="albuminuria",day1=FALSE,extra=NULL) {
  d$intake <- if(day1) d$day1_sodium else d$mean_sodium
  d$exp2 <- ifelse(d$intake<low_cut,"lower",ifelse(d$intake>=high_cut,"high",NA_character_))
  d <- droplevels(d[!is.na(d$exp2),,drop=FALSE]);d$exp2<-factor(d$exp2,levels=c("lower","high"))
  vars <- clean_vars(d,c("age","sex","race","bmi","egfr","diabetes","hypertension","cvd","sbp","hba1c","mean_kcal","mean_protein","mean_potassium","any_diuretic","smoking","pir","cycle"))
  if(!is.null(extra)) vars<-c(vars,extra)
  req <- setdiff(vars,"cycle")
  d <- droplevels(d[complete.cases(d[,req,drop=FALSE]),,drop=FALSE]);vars<-clean_vars(d,vars)
  if(nrow(d)<200||sum(d$exp2=="high")<40||sum(d[[outcome]]==1)<20||nlevels(d$exp2)<2)return(NULL)
  des<-make_design(d);f<-as.formula(paste(outcome,"~ exp2 +",paste(vars,collapse=" + ")))
  mpr<-svyglm(f,design=des,family=quasipoisson(link="log"));mrd<-svyglm(f,design=des,family=gaussian())
  mlog<-svyglm(as.formula(paste("log_uacr ~ exp2 +",paste(vars,collapse=" + "))),design=des,family=gaussian())
  zp<-safe_effect(mpr,"exp2high",TRUE);zr<-safe_effect(mrd,"exp2high");zg<-safe_effect(mlog,"exp2high",TRUE)
  sens_rows[[length(sens_rows)+1]] <<- data.frame(label=label,outcome=outcome,n=nrow(d),key_group_n=sum(d$exp2=="high"),key_group_events=sum(d[[outcome]][d$exp2=="high"]),pr=zp[1],pr_lower=zp[2],pr_upper=zp[3],rd=zr[1],rd_lower=zr[2],rd_upper=zr[3],uacr_gmr=zg[1],gmr_lower=zg[2],gmr_upper=zg[3],adverse_direction=is.finite(zp[1])&&is.finite(zr[1])&&is.finite(zg[1])&&zp[1]>1&&zr[1]>0&&zg[1]>1,stringsAsFactors=FALSE)
}
run_sens("primary",dat)
run_sens("high_ge_3500",dat,high_cut=3500)
run_sens("lower_lt_2000",dat,low_cut=2000)
run_sens("day1_sodium",dat,day1=TRUE)
run_sens("severe_albuminuria",dat,outcome="severe_albuminuria")
run_sens("egfr_ge_60",dat[dat$egfr>=60,])
run_sens("diabetes_only",dat[dat$diabetes==1,])
run_sens("non_diabetes",dat[dat$diabetes==0,])
run_sens("exclude_diuretics",dat[dat$any_diuretic==0,])
run_sens("sbp_lt_140",dat[dat$sbp<140,])
run_sens("hypertension_only",dat[dat$hypertension==1,])
run_sens("exclude_current_smokers",dat[dat$smoking!="current",])
run_sens("add_albumin",dat,extra="albumin")
for(cy in levels(dat$cycle)) run_sens(paste0("leave_out_",cy),dat[dat$cycle!=cy,])

model_df<-do.call(rbind,model_rows);sens_df<-do.call(rbind,sens_rows)
write.csv(model_df,file.path(out_dir,"c_n15_models.csv"),row.names=FALSE)
write.csv(sens_df,file.path(out_dir,"c_n15_sensitivities.csv"),row.names=FALSE)
write.csv(nonlinear_df,file.path(out_dir,"c_n15_nonlinearity.csv"),row.names=FALSE)
write.csv(missingness_df,file.path(out_dir,"c_n15_missingness_mechanism.csv"),row.names=FALSE)

prep<-fromJSON(file.path(out_dir,"c_n15_prep_status.json"));pr<-primary$pr;rd<-primary$rd;gmr<-primary$gmr;spr<-primary$severe_pr
dpr<-if(is.null(discovery))rep(NA_real_,4)else discovery$pr;vpr<-if(is.null(validation))rep(NA_real_,4)else validation$pr
drd<-if(is.null(discovery))rep(NA_real_,4)else discovery$rd;vrd<-if(is.null(validation))rep(NA_real_,4)else validation$rd
dg<-if(is.null(discovery))rep(NA_real_,4)else discovery$gmr;vg<-if(is.null(validation))rep(NA_real_,4)else validation$gmr
sample_gate<-prep$actual_n>=1500
retention_gate<-prep$complete_case_retention>=.80
key_gate<-prep$key_group_n>=500&&prep$key_group_events>=100&&prep$validation_key_group_n>=200&&prep$validation_key_group_events>=40
clinical_gate<-(pr[1]>=1.30&&pr[2]>=1.08&&rd[1]>=.04&&rd[2]>=.01)||(gmr[1]>=1.25&&gmr[2]>=1.08)
temporal_gate<-all(is.finite(c(dpr[1],vpr[1],drd[1],vrd[1],dg[1],vg[1])))&&dpr[1]>1&&vpr[1]>1&&drd[1]>0&&vrd[1]>0&&dg[1]>1&&vg[1]>1
sens_gate<-nrow(sens_df)>=10&&sum(sens_df$adverse_direction,na.rm=TRUE)>=ceiling(.70*nrow(sens_df))
nl_gate<-nrow(nonlinear_df)==2&&all(is.finite(nonlinear_df$p_nonlinear));miss_gate<-is.finite(missingness_df$odds_ratio[1])
pass<-all(c(sample_gate,retention_gate,key_gate,clinical_gate,temporal_gate,sens_gate,nl_gate,miss_gate))
status<-list(candidate_code="C-N15",data_analysis_completed=TRUE,prep=prep,primary=list(adjusted_albuminuria_pr=list(estimate=pr[1],lower=pr[2],upper=pr[3],p=pr[4]),adjusted_albuminuria_rd=list(estimate=rd[1],lower=rd[2],upper=rd[3],p=rd[4]),adjusted_uacr_gmr=list(estimate=gmr[1],lower=gmr[2],upper=gmr[3],p=gmr[4]),severe_albuminuria_pr=list(estimate=spr[1],lower=spr[2],upper=spr[3],p=spr[4])),discovery=list(pr=dpr[1],rd=drd[1],gmr=dg[1]),validation=list(pr=vpr[1],rd=vrd[1],gmr=vg[1]),nonlinearity=nonlinear_df,sensitivity_concordant=sum(sens_df$adverse_direction,na.rm=TRUE),sensitivity_total=nrow(sens_df),gates=list(sample_gate=sample_gate,complete_case_retention_gate=retention_gate,key_group_event_gate=key_gate,clinical_magnitude_precision_gate=clinical_gate,temporal_replication_gate=temporal_gate,sensitivity_direction_gate=sens_gate,nonlinearity_completed=nl_gate,missingness_mechanism_completed=miss_gate),data_gate_pass=pass,final_state=if(pass)"DATA_PASS_PENDING_FULL_TEXT_DEDUP"else"NO_GO_DATA_GATE")
write_json(status,file.path(out_dir,"c_n15_status.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
writeLines(c("# C-N15 blinded report",paste0("Actual N: ",prep$actual_n),paste0("Primary events: ",prep$primary_events),paste0("Key group/events: ",prep$key_group_n,"/",prep$key_group_events),paste0("Adjusted PR: ",round(pr[1],3)," [",round(pr[2],3),", ",round(pr[3],3),"]"),paste0("Adjusted RD: ",round(100*rd[1],2)," pp [",round(100*rd[2],2),", ",round(100*rd[3],2),"]"),paste0("UACR GMR: ",round(gmr[1],3)," [",round(gmr[2],3),", ",round(gmr[3],3),"]"),paste0("Discovery/validation PR: ",round(dpr[1],3),"/",round(vpr[1],3)),paste0("Sensitivities: ",sum(sens_df$adverse_direction,na.rm=TRUE),"/",nrow(sens_df)),paste0("State: ",status$final_state)),file.path(out_dir,"c_n15_report.md"))
cat(toJSON(status,pretty=TRUE,auto_unbox=TRUE),"\n")
