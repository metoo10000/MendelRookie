options(stringsAsFactors = FALSE)
options(survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(jsonlite)
})

out_dir <- "nutrition_stage/output_c_n13"
dat <- read.csv(file.path(out_dir, "c_n13_analysis_core.csv"), na.strings = c("", "NA"))
num_vars <- c("age","sex","race","education","pir","weight","psu_u","strata_u","mean_kcal","mean_protein","mean_potassium","mean_sodium","mean_magnesium","day1_potassium","day2_potassium","potassium_density","sodium_potassium_ratio","serum_potassium","serum_sodium","serum_chloride","serum_bicarbonate","creatinine","egfr","albumin","bmi","body_weight","height","diabetes","hypertension","cvd","heart_failure","wasting_days","low_k_3_6","hypokalemia_3_5","low_normal_k_4_0","log_potassium_intake")
for (v in intersect(num_vars, names(dat))) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
for (v in c("thiazide","loop","chronic_diuretic","raas","ppi","complete_case")) dat[[v]] <- as.logical(dat[[v]])
dat$cycle <- factor(dat$cycle, levels=c("E","F","G","H","I","J"))
dat$period <- factor(dat$period, levels=c("discovery","validation"))
dat$sex <- factor(dat$sex, levels=c(1,2), labels=c("male","female"))
dat$race <- factor(dat$race)
dat$diuretic_class <- factor(dat$diuretic_class, levels=c("thiazide","loop","both"))
dat$exposure <- factor(dat$exposure, levels=c("higher","low"))

make_design <- function(d) svydesign(ids=~psu_u, strata=~strata_u, weights=~weight, nest=TRUE, data=d)

safe_effect <- function(model, term, expit=FALSE) {
  cf <- coef(model)
  if (!(term %in% names(cf))) return(c(NA_real_,NA_real_,NA_real_,NA_real_))
  se <- sqrt(diag(vcov(model)))[term]
  est <- unname(cf[term])
  p <- tryCatch(unname(summary(model)$coefficients[term,ncol(summary(model)$coefficients)]), error=function(e) NA_real_)
  z <- c(est,est-1.96*se,est+1.96*se,p)
  if (expit) z[1:3] <- exp(z[1:3])
  z
}

model_rows <- list()
append_model <- function(scope,name,model,term="exposurelow",expit=FALSE) {
  z <- safe_effect(model,term,expit)
  model_rows[[length(model_rows)+1]] <<- data.frame(candidate_code="C-N13",scope=scope,model=name,term=term,estimate=z[1],lower=z[2],upper=z[3],p=z[4],n=nrow(model$model),stringsAsFactors=FALSE)
  z
}

clean_rhs <- function(d, vars) {
  out <- vars
  for (v in intersect(c("sex","race","diuretic_class","cycle"),out)) {
    if (length(unique(d[[v]][!is.na(d[[v]])])) < 2) out <- setdiff(out,v)
  }
  out
}

fit_scope <- function(d,scope) {
  d <- droplevels(d[d$complete_case & !is.na(d$exposure),,drop=FALSE])
  if (nrow(d)<250 || sum(d$exposure=="low")<50 || sum(d$low_k_3_6==1)<20 || nlevels(d$exposure)<2) return(NULL)
  des <- make_design(d)
  crude_vars <- clean_rhs(d,c("exposure","cycle"))
  full_vars <- clean_rhs(d,c("exposure","age","sex","race","bmi","egfr","diabetes","hypertension","cvd","mean_kcal","mean_sodium","mean_magnesium","diuretic_class","raas","pir","cycle"))
  f <- function(y,vars) as.formula(paste(y,"~",paste(vars,collapse=" + ")))
  crude_pr <- svyglm(f("low_k_3_6",crude_vars),design=des,family=quasipoisson(link="log"))
  crude_rd <- svyglm(f("low_k_3_6",crude_vars),design=des,family=gaussian())
  crude_k <- svyglm(f("serum_potassium",crude_vars),design=des,family=gaussian())
  adj_pr <- svyglm(f("low_k_3_6",full_vars),design=des,family=quasipoisson(link="log"))
  adj_rd <- svyglm(f("low_k_3_6",full_vars),design=des,family=gaussian())
  adj_k <- svyglm(f("serum_potassium",full_vars),design=des,family=gaussian())
  strict_pr <- svyglm(f("hypokalemia_3_5",full_vars),design=des,family=quasipoisson(link="log"))
  strict_rd <- svyglm(f("hypokalemia_3_5",full_vars),design=des,family=gaussian())
  append_model(scope,"crude_low_k_pr",crude_pr,expit=TRUE)
  append_model(scope,"crude_low_k_rd",crude_rd)
  append_model(scope,"crude_serum_k_difference",crude_k)
  pr <- append_model(scope,"adjusted_low_k_pr",adj_pr,expit=TRUE)
  rd <- append_model(scope,"adjusted_low_k_rd",adj_rd)
  kd <- append_model(scope,"adjusted_serum_k_difference",adj_k)
  spr <- append_model(scope,"adjusted_strict_hypokalemia_pr",strict_pr,expit=TRUE)
  srd <- append_model(scope,"adjusted_strict_hypokalemia_rd",strict_rd)
  list(data=d,design=des,pr=pr,rd=rd,kdiff=kd,strict_pr=spr,strict_rd=srd)
}

primary <- fit_scope(dat,"overall")
discovery <- fit_scope(dat[dat$period=="discovery",],"discovery")
validation <- fit_scope(dat[dat$period=="validation",],"validation")
if (is.null(primary)) stop("Primary C-N13 model could not be estimated")

# Missingness mechanism.
miss <- droplevels(dat[!is.na(dat$exposure)&!is.na(dat$age)&!is.na(dat$sex)&!is.na(dat$race),,drop=FALSE])
miss_vars <- clean_rhs(miss,c("exposure","age","sex","race","cycle"))
miss_mod <- svyglm(as.formula(paste("I(complete_case) ~",paste(miss_vars,collapse=" + "))),design=make_design(miss),family=quasibinomial())
miss_z <- safe_effect(miss_mod,"exposurelow",TRUE)
missingness_df <- data.frame(term="exposurelow",odds_ratio=miss_z[1],lower=miss_z[2],upper=miss_z[3],p=miss_z[4])

# Restricted cubic spline using four weighted knots.
weighted_quantile <- function(x,w,probs) {
  ok <- is.finite(x)&is.finite(w)&w>0; x<-x[ok];w<-w[ok];ord<-order(x);x<-x[ord];w<-w[ord];cw<-cumsum(w)/sum(w)
  sapply(probs,function(p)x[which(cw>=p)[1]])
}
rcs_basis <- function(x,k) {
  p3<-function(z)pmax(z,0)^3;den<-(k[length(k)]-k[1])^2
  h<-sapply(seq_len(length(k)-2),function(j)(p3(x-k[j])-p3(x-k[length(k)-1])*(k[length(k)]-k[j])/(k[length(k)]-k[length(k)-1])+p3(x-k[length(k)])*(k[length(k)-1]-k[j])/(k[length(k)]-k[length(k)-1]))/den)
  colnames(h)<-paste0("rcs",seq_len(ncol(h)));h
}
nl <- primary$data
knots <- weighted_quantile(nl$log_potassium_intake,nl$weight,c(.05,.35,.65,.95))
b<-rcs_basis(nl$log_potassium_intake,knots);nl$xlin<-nl$log_potassium_intake;nl$rcs1<-b[,1];nl$rcs2<-b[,2]
nl_vars <- clean_rhs(nl,c("xlin","rcs1","rcs2","age","sex","race","bmi","egfr","diabetes","hypertension","cvd","mean_kcal","mean_sodium","mean_magnesium","diuretic_class","raas","pir","cycle"))
f_nl <- function(y) as.formula(paste(y,"~",paste(nl_vars,collapse=" + ")))
des_nl<-make_design(nl)
nl_pr<-svyglm(f_nl("low_k_3_6"),design=des_nl,family=quasipoisson(link="log"))
nl_k<-svyglm(f_nl("serum_potassium"),design=des_nl,family=gaussian())
nonlinear_df<-rbind(
  data.frame(outcome="low_k_3_6",n=nrow(nl),p_nonlinear=unname(regTermTest(nl_pr,~rcs1+rcs2)$p),knots_mg_day=paste(round(exp(knots)),collapse="|")),
  data.frame(outcome="serum_potassium",n=nrow(nl),p_nonlinear=unname(regTermTest(nl_k,~rcs1+rcs2)$p),knots_mg_day=paste(round(exp(knots)),collapse="|"))
)

sens_rows <- list()
run_sens <- function(label,d,low_cut=2000,high_cut=2500,outcome="low_k_3_6",extra=NULL,use_day1=FALSE) {
  if (use_day1) d$intake <- d$day1_potassium else d$intake <- d$mean_potassium
  d$exp2 <- ifelse(d$intake<low_cut,"low",ifelse(d$intake>=high_cut,"higher",NA_character_))
  d <- droplevels(d[!is.na(d$exp2),,drop=FALSE]);d$exp2<-factor(d$exp2,levels=c("higher","low"))
  vars<-c("age","sex","race","bmi","egfr","diabetes","hypertension","cvd","mean_kcal","mean_sodium","mean_magnesium","diuretic_class","raas","pir","cycle")
  if (!is.null(extra)) vars<-c(vars,extra)
  vars<-clean_rhs(d,vars); req<-setdiff(vars,"cycle")
  d<-droplevels(d[complete.cases(d[,req,drop=FALSE]),,drop=FALSE]);vars<-clean_rhs(d,vars)
  if(nrow(d)<200||sum(d$exp2=="low")<40||sum(d[[outcome]]==1)<15||nlevels(d$exp2)<2)return(NULL)
  f<-as.formula(paste(outcome,"~ exp2 +",paste(vars,collapse=" + ")))
  des<-make_design(d);mpr<-svyglm(f,design=des,family=quasipoisson(link="log"));mrd<-svyglm(f,design=des,family=gaussian())
  fk<-as.formula(paste("serum_potassium ~ exp2 +",paste(vars,collapse=" + ")));mk<-svyglm(fk,design=des,family=gaussian())
  zpr<-safe_effect(mpr,"exp2low",TRUE);zrd<-safe_effect(mrd,"exp2low");zk<-safe_effect(mk,"exp2low")
  sens_rows[[length(sens_rows)+1]]<<-data.frame(label=label,outcome=outcome,n=nrow(d),key_group_n=sum(d$exp2=="low"),key_group_events=sum(d[[outcome]][d$exp2=="low"]),pr=zpr[1],pr_lower=zpr[2],pr_upper=zpr[3],rd=zrd[1],rd_lower=zrd[2],rd_upper=zrd[3],serum_k_difference=zk[1],k_lower=zk[2],k_upper=zk[3],adverse_direction=is.finite(zpr[1])&&is.finite(zrd[1])&&is.finite(zk[1])&&zpr[1]>1&&zrd[1]>0&&zk[1]<0,stringsAsFactors=FALSE)
}
run_sens("primary_reestimated",dat)
run_sens("strict_hypokalemia_3_5",dat,outcome="hypokalemia_3_5")
run_sens("low_1800_vs_2500",dat,low_cut=1800)
run_sens("low_2000_vs_3000",dat,high_cut=3000)
run_sens("day1_intake",dat,use_day1=TRUE)
run_sens("egfr_ge_60",dat[dat$egfr>=60,])
run_sens("exclude_raas",dat[!dat$raas,])
run_sens("chronic_diuretic_30d",dat[dat$chronic_diuretic,])
run_sens("exclude_heart_failure",dat[dat$heart_failure==0,])
run_sens("thiazide_only",dat[dat$diuretic_class=="thiazide",])
run_sens("loop_only",dat[dat$diuretic_class=="loop",])
run_sens("women_only",dat[dat$sex=="female",])
run_sens("men_only",dat[dat$sex=="male",])
run_sens("add_bicarbonate",dat,extra="serum_bicarbonate")
for(cy in levels(dat$cycle))run_sens(paste0("leave_out_",cy),dat[dat$cycle!=cy,])

model_df<-do.call(rbind,model_rows);sens_df<-do.call(rbind,sens_rows)
write.csv(model_df,file.path(out_dir,"c_n13_models.csv"),row.names=FALSE)
write.csv(sens_df,file.path(out_dir,"c_n13_sensitivities.csv"),row.names=FALSE)
write.csv(nonlinear_df,file.path(out_dir,"c_n13_nonlinearity.csv"),row.names=FALSE)
write.csv(missingness_df,file.path(out_dir,"c_n13_missingness_mechanism.csv"),row.names=FALSE)

prep<-fromJSON(file.path(out_dir,"c_n13_prep_status.json"));pr<-primary$pr;rd<-primary$rd;kd<-primary$kdiff;spr<-primary$strict_pr
dpr<-if(is.null(discovery))rep(NA_real_,4)else discovery$pr;vpr<-if(is.null(validation))rep(NA_real_,4)else validation$pr
drd<-if(is.null(discovery))rep(NA_real_,4)else discovery$rd;vrd<-if(is.null(validation))rep(NA_real_,4)else validation$rd
dkd<-if(is.null(discovery))rep(NA_real_,4)else discovery$kdiff;vkd<-if(is.null(validation))rep(NA_real_,4)else validation$kdiff
sample_gate<-prep$actual_n>=2500
retention_gate<-prep$complete_case_retention>=.80
key_gate<-prep$key_group_n>=800&&prep$key_group_events>=60&&prep$validation_key_group_n>=300&&prep$validation_key_group_events>=20
clinical_gate<-(is.finite(pr[1])&&pr[1]>=1.50&&pr[2]>=1.10&&is.finite(rd[1])&&rd[1]>=.02&&rd[2]>.002)||(is.finite(kd[1])&&kd[1]<=-.10&&kd[3]<=-.03)
temporal_gate<-all(is.finite(c(dpr[1],vpr[1],drd[1],vrd[1],dkd[1],vkd[1])))&&dpr[1]>1&&vpr[1]>1&&drd[1]>0&&vrd[1]>0&&dkd[1]<0&&vkd[1]<0
sens_gate<-nrow(sens_df)>=10&&sum(sens_df$adverse_direction,na.rm=TRUE)>=ceiling(.70*nrow(sens_df))
nonlinear_complete<-nrow(nonlinear_df)==2&&all(is.finite(nonlinear_df$p_nonlinear));missing_complete<-is.finite(missingness_df$odds_ratio[1])
pass<-all(c(sample_gate,retention_gate,key_gate,clinical_gate,temporal_gate,sens_gate,nonlinear_complete,missing_complete))
status<-list(candidate_code="C-N13",data_analysis_completed=TRUE,prep=prep,primary=list(adjusted_low_k_pr=list(estimate=pr[1],lower=pr[2],upper=pr[3],p=pr[4]),adjusted_low_k_rd=list(estimate=rd[1],lower=rd[2],upper=rd[3],p=rd[4]),adjusted_serum_k_difference=list(estimate=kd[1],lower=kd[2],upper=kd[3],p=kd[4]),strict_hypokalemia_pr=list(estimate=spr[1],lower=spr[2],upper=spr[3],p=spr[4])),discovery=list(pr=dpr[1],rd=drd[1],k_difference=dkd[1]),validation=list(pr=vpr[1],rd=vrd[1],k_difference=vkd[1]),nonlinearity=nonlinear_df,sensitivity_concordant=sum(sens_df$adverse_direction,na.rm=TRUE),sensitivity_total=nrow(sens_df),gates=list(sample_gate=sample_gate,complete_case_retention_gate=retention_gate,key_group_event_gate=key_gate,clinical_magnitude_precision_gate=clinical_gate,temporal_replication_gate=temporal_gate,sensitivity_direction_gate=sens_gate,nonlinearity_completed=nonlinear_complete,missingness_mechanism_completed=missing_complete),data_gate_pass=pass,final_state=if(pass)"DATA_PASS_PENDING_FULL_TEXT_DEDUP"else"NO_GO_DATA_GATE")
write_json(status,file.path(out_dir,"c_n13_status.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
writeLines(c("# C-N13 blinded data-gate report",paste0("Actual N: ",prep$actual_n),paste0("Primary events: ",prep$primary_events),paste0("Key group/events: ",prep$key_group_n,"/",prep$key_group_events),paste0("Adjusted PR: ",round(pr[1],3)," [",round(pr[2],3),", ",round(pr[3],3),"]"),paste0("Adjusted RD: ",round(100*rd[1],2)," pp [",round(100*rd[2],2),", ",round(100*rd[3],2),"]"),paste0("Serum K difference: ",round(kd[1],3)," [",round(kd[2],3),", ",round(kd[3],3),"]"),paste0("Discovery/validation PR: ",round(dpr[1],3),"/",round(vpr[1],3)),paste0("Sensitivities: ",sum(sens_df$adverse_direction,na.rm=TRUE),"/",nrow(sens_df)),paste0("State: ",status$final_state)),file.path(out_dir,"c_n13_report.md"))
cat(toJSON(status,pretty=TRUE,auto_unbox=TRUE),"\n")
