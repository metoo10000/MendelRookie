options(stringsAsFactors = FALSE)
options(survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(jsonlite)
})

out_dir <- "nutrition_stage/output_batch_n11"
d <- read.csv(file.path(out_dir, "batch_n11_core.csv"), na.strings = c("", "NA"))

numeric_vars <- c("age","sex","race","weight","psu_u","strata_u","pir","bmi","diabetes","hypertension","cvd","egfr","albumin","energy1","potassium1","sodium1","potassium_mean2","sodium_mean2","k_serum","na_serum","hypokalemia","hyponatremia","hyperkalemia","low_k_intake","low_na_intake","high_k_intake","k_per_1000kcal","na_per_1000kcal")
for (v in intersect(numeric_vars, names(d))) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
for (v in c("diuretic","thiazide","loop","k_sparing","potassium_rx","ace_arb","complete_case_k","complete_case_na")) d[[v]] <- as.logical(d[[v]])
d$cycle <- factor(d$cycle, levels = c("D","E","F","G","H","I","J"))
d$period <- factor(d$period, levels = c("discovery","validation"))
d$race <- factor(d$race)
d$sex <- factor(d$sex)
d$smoking <- factor(d$smoking, levels = c("never","former","current"))

make_design <- function(x) svydesign(ids = ~psu_u, strata = ~strata_u, weights = ~weight, nest = TRUE, data = x)
ci_term <- function(model, term, exponentiate = FALSE) {
  if (!(term %in% names(coef(model)))) return(c(est=NA, lo=NA, hi=NA, p=NA))
  est <- unname(coef(model)[term]); se <- sqrt(diag(vcov(model)))[term]
  z <- c(est=est, lo=est-1.96*se, hi=est+1.96*se, p=summary(model)$coefficients[term,ncol(summary(model)$coefficients)])
  if (exponentiate) z[1:3] <- exp(z[1:3])
  z
}

weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]; ord <- order(x); x <- x[ord]; w <- w[ord]
  cw <- cumsum(w)/sum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}
rcs_basis <- function(x, knots) {
  k <- knots; pos3 <- function(z) pmax(z,0)^3; den <- (k[length(k)]-k[1])^2
  h <- sapply(seq_len(length(k)-2), function(j) {
    (pos3(x-k[j]) - pos3(x-k[length(k)-1])*(k[length(k)]-k[j])/(k[length(k)]-k[length(k)-1]) + pos3(x-k[length(k)])*(k[length(k)-1]-k[j])/(k[length(k)]-k[length(k)-1]))/den
  })
  colnames(h) <- paste0("rcs", seq_len(ncol(h))); h
}

all_models <- list(); all_sens <- list(); all_nl <- list(); all_miss <- list(); statuses <- list()

candidate_specs <- list(
  list(code="C-N11", subset=quote(diuretic & !k_sparing), cc="complete_case_k", exposure="low_k_intake", continuous="potassium1", ycont="k_serum", ybin="hypokalemia", diff_dir=-1, diff_cut=0.10, diff_ci=0.03, pr_cut=1.40, pr_lo=1.10, sample_n=2500, key_n=500, event_n=100, val_event_n=30, cov_other="sodium1"),
  list(code="C-N12", subset=quote(thiazide & !k_sparing), cc="complete_case_na", exposure="low_na_intake", continuous="sodium1", ycont="na_serum", ybin="hyponatremia", diff_dir=-1, diff_cut=1.00, diff_ci=0.30, pr_cut=1.50, pr_lo=1.15, sample_n=2000, key_n=400, event_n=80, val_event_n=25, cov_other="potassium1"),
  list(code="C-N13", subset=quote(ace_arb & !diuretic), cc="complete_case_k", exposure="high_k_intake", continuous="potassium1", ycont="k_serum", ybin="hyperkalemia", diff_dir=1, diff_cut=0.10, diff_ci=0.03, pr_cut=1.50, pr_lo=1.15, sample_n=2500, key_n=500, event_n=80, val_event_n=25, cov_other="sodium1")
)

fit_one <- function(spec) {
  code <- spec$code
  base <- d[eval(spec$subset, d), , drop=FALSE]
  base <- base[base[[spec$cc]] & is.finite(base[[spec$exposure]]) & is.finite(base[[spec$ycont]]) & is.finite(base[[spec$ybin]]), , drop=FALSE]
  base$expf <- factor(base[[spec$exposure]], levels=c(0,1))
  term <- "expf1"
  covars <- c("age","sex","race","bmi","egfr","diabetes","hypertension","cvd","smoking","energy1",spec$cov_other,"cycle")
  cov_string <- paste(covars, collapse=" + ")

  fit_scope <- function(x, scope) {
    x <- droplevels(x)
    if (nrow(x)<100 || length(unique(x$psu_u))<4 || length(unique(x$expf))<2) return(NULL)
    des <- make_design(x)
    m_crude_cont <- svyglm(as.formula(paste(spec$ycont, "~ expf + cycle")), design=des, family=gaussian())
    m_crude_bin <- svyglm(as.formula(paste(spec$ybin, "~ expf + cycle")), design=des, family=quasipoisson(link="log"))
    m_adj_cont <- svyglm(as.formula(paste(spec$ycont, "~ expf +", cov_string)), design=des, family=gaussian())
    m_adj_bin <- svyglm(as.formula(paste(spec$ybin, "~ expf +", cov_string)), design=des, family=quasipoisson(link="log"))
    vals <- list(
      crude_diff=ci_term(m_crude_cont,term,FALSE), crude_pr=ci_term(m_crude_bin,term,TRUE),
      adjusted_diff=ci_term(m_adj_cont,term,FALSE), adjusted_pr=ci_term(m_adj_bin,term,TRUE)
    )
    for (nm in names(vals)) {
      z <- vals[[nm]]
      all_models[[length(all_models)+1]] <<- data.frame(candidate_code=code,scope=scope,model=nm,estimate=z[1],lower=z[2],upper=z[3],p=z[4],n=nrow(x),events=sum(x[[spec$ybin]]),key_n=sum(x[[spec$exposure]]==1),key_events=sum(x[[spec$ybin]][x[[spec$exposure]]==1]),stringsAsFactors=FALSE)
    }
    vals
  }

  overall <- fit_scope(base,"overall")
  discovery <- fit_scope(base[base$period=="discovery",],"discovery")
  validation <- fit_scope(base[base$period=="validation",],"validation")
  if (is.null(overall)) return(list(code=code,final_state="NO_GO_MODEL_NOT_ESTIMABLE",data_gate_pass=FALSE))

  # Missingness mechanism in the broader medication-defined source population.
  source <- d[eval(spec$subset,d),,drop=FALSE]
  source$expf <- factor(source[[spec$exposure]], levels=c(0,1))
  source$observed <- as.numeric(source[[spec$cc]])
  mm <- source[is.finite(source$observed) & is.finite(source[[spec$exposure]]) & is.finite(source$age) & !is.na(source$race),]
  if (nrow(mm)>100 && length(unique(mm$expf))==2) {
    md <- make_design(mm)
    mf <- svyglm(observed ~ expf + age + race + cycle, design=md, family=quasibinomial())
    z <- ci_term(mf,"expf1",TRUE)
    all_miss[[length(all_miss)+1]] <<- data.frame(candidate_code=code,term="exposure_predicts_complete_case",odds_ratio=z[1],lower=z[2],upper=z[3],p=z[4],n=nrow(mm))
  }

  # Restricted cubic spline for continuous nutrient exposure, with a nonlinear joint test.
  nl <- base[is.finite(base[[spec$continuous]]),]
  knots <- weighted_quantile(nl[[spec$continuous]],nl$weight,c(.05,.35,.65,.95))
  basis <- rcs_basis(nl[[spec$continuous]],knots)
  nl$xlin <- nl[[spec$continuous]]/1000; nl$rcs1 <- basis[,1]/1e9; nl$rcs2 <- basis[,2]/1e9
  nd <- make_design(nl)
  fm_cont <- as.formula(paste(spec$ycont,"~ xlin + rcs1 + rcs2 +",cov_string))
  fm_bin <- as.formula(paste(spec$ybin,"~ xlin + rcs1 + rcs2 +",cov_string))
  mn_cont <- svyglm(fm_cont,design=nd,family=gaussian())
  mn_bin <- svyglm(fm_bin,design=nd,family=quasipoisson(link="log"))
  pn_cont <- tryCatch(unname(regTermTest(mn_cont,~rcs1+rcs2)$p),error=function(e) NA_real_)
  pn_bin <- tryCatch(unname(regTermTest(mn_bin,~rcs1+rcs2)$p),error=function(e) NA_real_)
  all_nl[[length(all_nl)+1]] <<- data.frame(candidate_code=code,outcome=spec$ycont,n=nrow(nl),p_nonlinear=pn_cont,knots=paste(round(knots,1),collapse="|"))
  all_nl[[length(all_nl)+1]] <<- data.frame(candidate_code=code,outcome=spec$ybin,n=nrow(nl),p_nonlinear=pn_bin,knots=paste(round(knots,1),collapse="|"))

  run_sens <- function(label,x,exp_var=spec$exposure,extra_cov=NULL) {
    x <- x[x[[spec$cc]] & is.finite(x[[exp_var]]) & is.finite(x[[spec$ycont]]) & is.finite(x[[spec$ybin]]),,drop=FALSE]
    x$exp2 <- factor(x[[exp_var]],levels=c(0,1))
    if (nrow(x)<100 || length(unique(x$exp2))<2 || sum(x[[spec$ybin]])<10) return()
    des <- make_design(x)
    fs <- paste(c(covars,extra_cov),collapse=" + ")
    mc <- svyglm(as.formula(paste(spec$ycont,"~ exp2 +",fs)),design=des,family=gaussian())
    mb <- svyglm(as.formula(paste(spec$ybin,"~ exp2 +",fs)),design=des,family=quasipoisson(link="log"))
    zd <- ci_term(mc,"exp21",FALSE); zp <- ci_term(mb,"exp21",TRUE)
    concord <- if (spec$diff_dir<0) is.finite(zd[1])&&is.finite(zp[1])&&zd[1]<0&&zp[1]>1 else is.finite(zd[1])&&is.finite(zp[1])&&zd[1]>0&&zp[1]>1
    all_sens[[length(all_sens)+1]] <<- data.frame(candidate_code=code,label=label,n=nrow(x),events=sum(x[[spec$ybin]]),key_n=sum(x[[exp_var]]==1),key_events=sum(x[[spec$ybin]][x[[exp_var]]==1]),difference=zd[1],diff_lower=zd[2],diff_upper=zd[3],pr=zp[1],pr_lower=zp[2],pr_upper=zp[3],direction_concordant=concord,stringsAsFactors=FALSE)
  }

  run_sens("primary_reestimated",base)
  run_sens("exclude_potassium_prescription",base[!base$potassium_rx,])
  run_sens("exclude_egfr_lt_30",base[base$egfr>=30,])
  run_sens("exclude_cvd",base[base$cvd==0,])
  run_sens("add_albumin",base[is.finite(base$albumin),],extra_cov="albumin")
  run_sens("add_pir",base[is.finite(base$pir),],extra_cov="pir")
  for (cy in levels(d$cycle)) run_sens(paste0("leave_out_",cy),base[base$cycle!=cy,])
  if (code=="C-N11") {
    base$alt2500 <- as.numeric(base$potassium1<2500); run_sens("potassium_cut_2500",base,"alt2500")
    base$densitylow <- as.numeric(base$k_per_1000kcal<1000); run_sens("potassium_density_cut",base,"densitylow")
    run_sens("thiazide_only",base[base$thiazide & !base$loop,])
    run_sens("loop_only",base[base$loop,])
  }
  if (code=="C-N12") {
    base$alt2300 <- as.numeric(base$sodium1<2300); run_sens("sodium_cut_2300",base,"alt2300")
    base$densitylow <- as.numeric(base$na_per_1000kcal<1000); run_sens("sodium_density_cut",base,"densitylow")
  }
  if (code=="C-N13") {
    base$alt3000 <- as.numeric(base$potassium1>=3000); run_sens("potassium_cut_3000",base,"alt3000")
    run_sens("exclude_ckd",base[base$egfr>=60,])
  }

  sr <- do.call(rbind,all_sens)
  sr <- sr[sr$candidate_code==code,]
  diff <- overall$adjusted_diff; pr <- overall$adjusted_pr
  n <- nrow(base); events <- sum(base[[spec$ybin]]); key_n <- sum(base[[spec$exposure]]==1); key_events <- sum(base[[spec$ybin]][base[[spec$exposure]]==1])
  val <- base[base$period=="validation",]; val_events <- sum(val[[spec$ybin]]); val_key_events <- sum(val[[spec$ybin]][val[[spec$exposure]]==1])
  source_n <- nrow(d[eval(spec$subset,d),,drop=FALSE]); retention <- n/max(source_n,1)
  sample_gate <- n>=spec$sample_n
  retention_gate <- retention>=0.80
  event_gate <- events>=spec$event_n && key_n>=spec$key_n && key_events>=40 && val_events>=spec$val_event_n && val_key_events>=10
  if (spec$diff_dir<0) {
    clinical_diff <- is.finite(diff[1]) && diff[1]<=-spec$diff_cut && diff[3]<=-spec$diff_ci
    temporal <- !is.null(discovery)&&!is.null(validation)&&discovery$adjusted_diff[1]<0&&validation$adjusted_diff[1]<0&&discovery$adjusted_pr[1]>1&&validation$adjusted_pr[1]>1
  } else {
    clinical_diff <- is.finite(diff[1]) && diff[1]>=spec$diff_cut && diff[2]>=spec$diff_ci
    temporal <- !is.null(discovery)&&!is.null(validation)&&discovery$adjusted_diff[1]>0&&validation$adjusted_diff[1]>0&&discovery$adjusted_pr[1]>1&&validation$adjusted_pr[1]>1
  }
  clinical_pr <- is.finite(pr[1]) && pr[1]>=spec$pr_cut && pr[2]>=spec$pr_lo
  clinical_gate <- clinical_diff || clinical_pr
  sensitivity_gate <- nrow(sr)>=8 && mean(sr$direction_concordant,na.rm=TRUE)>=0.70
  nl_gate <- all(is.finite(c(pn_cont,pn_bin)))
  gates <- list(sample_gate=sample_gate,complete_case_retention_gate=retention_gate,key_group_event_gate=event_gate,clinical_magnitude_precision_gate=clinical_gate,temporal_replication_gate=temporal,sensitivity_direction_gate=sensitivity_gate,nonlinearity_completed=nl_gate)
  pass <- all(unlist(gates))
  list(code=code,actual_n=n,source_eligible_n=source_n,complete_case_retention=retention,events=events,key_group_n=key_n,key_group_events=key_events,validation_events=val_events,validation_key_events=val_key_events,adjusted_difference=list(estimate=diff[1],lower=diff[2],upper=diff[3],p=diff[4]),adjusted_pr=list(estimate=pr[1],lower=pr[2],upper=pr[3],p=pr[4]),discovery=list(difference=if(is.null(discovery)) NA else discovery$adjusted_diff[1],pr=if(is.null(discovery)) NA else discovery$adjusted_pr[1]),validation=list(difference=if(is.null(validation)) NA else validation$adjusted_diff[1],pr=if(is.null(validation)) NA else validation$adjusted_pr[1]),sensitivity_concordant=sum(sr$direction_concordant,na.rm=TRUE),sensitivity_total=nrow(sr),gates=gates,data_gate_pass=pass,final_state=if(pass)"DATA_PASS_PENDING_FULL_TEXT_DEDUP" else "NO_GO_DATA_GATE")
}

for (spec in candidate_specs) statuses[[spec$code]] <- fit_one(spec)
models_df <- if(length(all_models)) do.call(rbind,all_models) else data.frame()
sens_df <- if(length(all_sens)) do.call(rbind,all_sens) else data.frame()
nl_df <- if(length(all_nl)) do.call(rbind,all_nl) else data.frame()
miss_df <- if(length(all_miss)) do.call(rbind,all_miss) else data.frame()
write.csv(models_df,file.path(out_dir,"models.csv"),row.names=FALSE)
write.csv(sens_df,file.path(out_dir,"sensitivities.csv"),row.names=FALSE)
write.csv(nl_df,file.path(out_dir,"nonlinearity.csv"),row.names=FALSE)
write.csv(miss_df,file.path(out_dir,"missingness_mechanism.csv"),row.names=FALSE)

passes <- names(statuses)[vapply(statuses,function(x)isTRUE(x$data_gate_pass),logical(1))]
batch <- list(batch="batch_n11",analysis_completed=TRUE,candidates=statuses,data_pass_candidates=passes,selected_for_dedup=if(length(passes))passes[[1]] else NULL,final_state=if(length(passes))"DATA_PASS_FOUND" else "NO_DATA_PASS")
write_json(batch,file.path(out_dir,"batch_status.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
cat(toJSON(batch,pretty=TRUE,auto_unbox=TRUE,na="null"),"\n")
