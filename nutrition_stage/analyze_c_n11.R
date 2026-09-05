options(stringsAsFactors = FALSE)
options(survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(jsonlite)
})

out_dir <- "nutrition_stage/output_c_n11"
dat <- read.csv(file.path(out_dir, "c_n11_analysis_core.csv"), na.strings = c("", "NA"))
num_vars <- c("age", "sex", "race", "pir", "weight", "psu_u", "strata_u", "vitamin_c_umol_l", "ferritin", "stfr", "hgb", "mcv", "rdw", "mch", "crp_mg_l", "creatinine", "egfr", "albumin", "total_protein", "alt", "ast", "bmi", "body_weight", "height", "diabetes", "known_cancer", "known_liver", "anemia", "high_rdw", "microcytosis", "log_vitamin_c", "log_crp", "log_ferritin")
for (v in intersect(num_vars, names(dat))) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
dat$complete_case <- as.logical(dat$complete_case)
dat$cycle <- factor(dat$cycle, levels = c("C", "D"))
dat$period <- factor(dat$period, levels = c("discovery", "validation"))
dat$sex <- factor(dat$sex, levels = c(1, 2), labels = c("male", "female"))
dat$race <- factor(dat$race)
dat$smoking <- factor(dat$smoking, levels = c("never", "former", "current"))
dat$exposure <- factor(dat$exposure, levels = c("adequate", "depleted"))

make_design <- function(d) svydesign(ids = ~psu_u, strata = ~strata_u, weights = ~weight, nest = TRUE, data = d)

safe_effect <- function(model, term, exponentiate = FALSE) {
  cf <- coef(model)
  if (!(term %in% names(cf))) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  se <- sqrt(diag(vcov(model)))[term]
  est <- unname(cf[term])
  p <- tryCatch(unname(summary(model)$coefficients[term, ncol(summary(model)$coefficients)]), error = function(e) NA_real_)
  z <- c(est, est - 1.96 * se, est + 1.96 * se, p)
  if (exponentiate) z[1:3] <- exp(z[1:3])
  z
}

model_rows <- list()
append_model <- function(scope, name, model, term = "exposuredepleted", exponentiate = FALSE) {
  z <- safe_effect(model, term, exponentiate)
  model_rows[[length(model_rows) + 1]] <<- data.frame(candidate_code = "C-N11", scope = scope, model = name, term = term, estimate = z[1], lower = z[2], upper = z[3], p = z[4], n = nrow(model$model), stringsAsFactors = FALSE)
  z
}

fit_scope <- function(d, scope) {
  d <- d[d$complete_case & !is.na(d$exposure), , drop = FALSE]
  d$exposure <- droplevels(factor(d$exposure, levels = c("adequate", "depleted")))
  if (nrow(d) < 250 || sum(d$exposure == "depleted") < 30 || sum(d$anemia == 1) < 20 || nlevels(d$exposure) < 2) return(NULL)
  des <- make_design(d)
  core_rhs <- "exposure + age + sex + race + bmi + log_crp + log_ferritin + albumin + smoking + diabetes + pir"
  if (length(unique(d$cycle)) > 1) core_rhs <- paste(core_rhs, "+ cycle")
  crude_rhs <- "exposure"
  if (length(unique(d$cycle)) > 1) crude_rhs <- paste(crude_rhs, "+ cycle")
  crude_hgb <- svyglm(as.formula(paste("hgb ~", crude_rhs)), design = des, family = gaussian())
  crude_anemia <- svyglm(as.formula(paste("anemia ~", crude_rhs)), design = des, family = quasipoisson(link = "log"))
  adj_hgb <- svyglm(as.formula(paste("hgb ~", core_rhs)), design = des, family = gaussian())
  adj_anemia <- svyglm(as.formula(paste("anemia ~", core_rhs)), design = des, family = quasipoisson(link = "log"))
  adj_rdw <- svyglm(as.formula(paste("high_rdw ~", core_rhs)), design = des, family = quasipoisson(link = "log"))
  adj_mcv <- svyglm(as.formula(paste("mcv ~", core_rhs)), design = des, family = gaussian())
  append_model(scope, "crude_hgb_difference_g_dl", crude_hgb)
  append_model(scope, "crude_anemia_prevalence_ratio", crude_anemia, exponentiate = TRUE)
  hgb <- append_model(scope, "adjusted_hgb_difference_g_dl", adj_hgb)
  anemia <- append_model(scope, "adjusted_anemia_prevalence_ratio", adj_anemia, exponentiate = TRUE)
  rdw <- append_model(scope, "adjusted_high_rdw_prevalence_ratio", adj_rdw, exponentiate = TRUE)
  mcv <- append_model(scope, "adjusted_mcv_difference_fl", adj_mcv)
  list(data = d, design = des, hgb = hgb, anemia = anemia, rdw = rdw, mcv = mcv)
}

primary <- fit_scope(dat, "overall")
discovery <- fit_scope(dat[dat$period == "discovery", ], "discovery")
validation <- fit_scope(dat[dat$period == "validation", ], "validation")
if (is.null(primary)) stop("Primary C-N11 model could not be estimated")

# Missingness mechanism
miss <- dat[!is.na(dat$exposure) & !is.na(dat$age) & !is.na(dat$sex) & !is.na(dat$race), , drop = FALSE]
miss_des <- make_design(miss)
miss_rhs <- "exposure + age + sex + race"
if (length(unique(miss$cycle)) > 1) miss_rhs <- paste(miss_rhs, "+ cycle")
miss_mod <- svyglm(as.formula(paste("I(complete_case) ~", miss_rhs)), design = miss_des, family = quasibinomial())
miss_z <- safe_effect(miss_mod, "exposuredepleted", exponentiate = TRUE)
missingness_df <- data.frame(term = "exposuredepleted", odds_ratio = miss_z[1], lower = miss_z[2], upper = miss_z[3], p = miss_z[4])

weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]; ord <- order(x); x <- x[ord]; w <- w[ord]
  cw <- cumsum(w) / sum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}
rcs_basis <- function(x, k) {
  p3 <- function(z) pmax(z, 0)^3
  den <- (k[length(k)] - k[1])^2
  h <- sapply(seq_len(length(k) - 2), function(j) (p3(x-k[j]) - p3(x-k[length(k)-1])*(k[length(k)]-k[j])/(k[length(k)]-k[length(k)-1]) + p3(x-k[length(k)])*(k[length(k)-1]-k[j])/(k[length(k)]-k[length(k)-1]))/den)
  colnames(h) <- paste0("rcs", seq_len(ncol(h))); h
}

nl <- primary$data
knots <- weighted_quantile(nl$log_vitamin_c, nl$weight, c(.05, .35, .65, .95))
b <- rcs_basis(nl$log_vitamin_c, knots)
nl$xlin <- nl$log_vitamin_c; nl$rcs1 <- b[,1]; nl$rcs2 <- b[,2]
des_nl <- make_design(nl)
nl_rhs <- "xlin + rcs1 + rcs2 + age + sex + race + bmi + log_crp + log_ferritin + albumin + smoking + diabetes + pir + cycle"
nl_hgb <- svyglm(as.formula(paste("hgb ~", nl_rhs)), design = des_nl, family = gaussian())
nl_anemia <- svyglm(as.formula(paste("anemia ~", nl_rhs)), design = des_nl, family = quasipoisson(link = "log"))
nonlinear_df <- rbind(
  data.frame(outcome = "hemoglobin", n = nrow(nl), p_nonlinear = unname(regTermTest(nl_hgb, ~rcs1 + rcs2)$p), knots_umol_l = paste(round(exp(knots),2), collapse="|")),
  data.frame(outcome = "anemia", n = nrow(nl), p_nonlinear = unname(regTermTest(nl_anemia, ~rcs1 + rcs2)$p), knots_umol_l = paste(round(exp(knots),2), collapse="|"))
)

sens_rows <- list()
run_sens <- function(label, d, low_cut = 23, adequate_cut = 50, ferritin_lo = 30, ferritin_hi = 500, crp_max = 10, egfr_min = 60, extra = NULL) {
  d <- d[d$ferritin >= ferritin_lo & d$ferritin <= ferritin_hi & d$crp_mg_l <= crp_max & d$egfr >= egfr_min, , drop=FALSE]
  d$exp2 <- ifelse(d$vitamin_c_umol_l < low_cut, "depleted", ifelse(d$vitamin_c_umol_l >= adequate_cut, "adequate", NA_character_))
  req <- c("age","sex","race","bmi","log_crp","log_ferritin","albumin","smoking","diabetes","pir")
  d <- d[!is.na(d$exp2) & complete.cases(d[,req,drop=FALSE]),,drop=FALSE]
  d$exp2 <- factor(d$exp2, levels=c("adequate","depleted"))
  if (nrow(d)<200 || sum(d$exp2=="depleted")<25 || sum(d$anemia==1)<15) return(NULL)
  rhs <- c("exp2","age","sex","race","bmi","log_crp","log_ferritin","albumin","smoking","diabetes","pir")
  if (length(unique(d$cycle))>1) rhs <- c(rhs,"cycle")
  if (!is.null(extra)) rhs <- c(rhs, extra)
  des <- make_design(d)
  mh <- svyglm(as.formula(paste("hgb ~",paste(rhs,collapse=" + "))),design=des,family=gaussian())
  ma <- svyglm(as.formula(paste("anemia ~",paste(rhs,collapse=" + "))),design=des,family=quasipoisson(link="log"))
  mr <- svyglm(as.formula(paste("high_rdw ~",paste(rhs,collapse=" + "))),design=des,family=quasipoisson(link="log"))
  zh <- safe_effect(mh,"exp2depleted"); za <- safe_effect(ma,"exp2depleted",TRUE); zr <- safe_effect(mr,"exp2depleted",TRUE)
  sens_rows[[length(sens_rows)+1]] <<- data.frame(label=label,n=nrow(d),key_group_n=sum(d$exp2=="depleted"),key_group_anemia=sum(d$anemia[d$exp2=="depleted"]),hgb_difference=zh[1],hgb_lower=zh[2],hgb_upper=zh[3],anemia_pr=za[1],anemia_lower=za[2],anemia_upper=za[3],high_rdw_pr=zr[1],rdw_lower=zr[2],rdw_upper=zr[3],adverse_direction=is.finite(zh[1])&&is.finite(za[1])&&zh[1]<0&&za[1]>1,stringsAsFactors=FALSE)
}

run_sens("primary_reestimated",dat)
run_sens("severe_depletion_lt_11_4",dat,low_cut=11.4)
run_sens("adequate_ge_45",dat,adequate_cut=45)
run_sens("strict_ferritin_50_300",dat,ferritin_lo=50,ferritin_hi=300)
run_sens("low_inflammation_crp_le_3",dat,crp_max=3)
run_sens("normal_kidney_egfr_ge_90",dat,egfr_min=90)
run_sens("exclude_diabetes",dat[dat$diabetes==0,])
run_sens("exclude_current_smokers",dat[dat$smoking!="current",])
run_sens("non_obese",dat[dat$bmi<30,])
run_sens("women_only",dat[dat$sex=="female",])
run_sens("men_only",dat[dat$sex=="male",])
run_sens("adjust_total_protein",dat,extra="total_protein")

model_df <- do.call(rbind, model_rows)
sens_df <- do.call(rbind, sens_rows)
write.csv(model_df,file.path(out_dir,"c_n11_models.csv"),row.names=FALSE)
write.csv(sens_df,file.path(out_dir,"c_n11_sensitivities.csv"),row.names=FALSE)
write.csv(nonlinear_df,file.path(out_dir,"c_n11_nonlinearity.csv"),row.names=FALSE)
write.csv(missingness_df,file.path(out_dir,"c_n11_missingness_mechanism.csv"),row.names=FALSE)

prep <- fromJSON(file.path(out_dir,"c_n11_prep_status.json"))
ph <- primary$hgb; pa <- primary$anemia; pr <- primary$rdw; pm <- primary$mcv
dh <- if(is.null(discovery)) rep(NA_real_,4) else discovery$hgb
da <- if(is.null(discovery)) rep(NA_real_,4) else discovery$anemia
vh <- if(is.null(validation)) rep(NA_real_,4) else validation$hgb
va <- if(is.null(validation)) rep(NA_real_,4) else validation$anemia
sample_gate <- prep$actual_n >= 2500
retention_gate <- prep$complete_case_retention >= .80
key_gate <- prep$key_group_n >= 250 && prep$key_group_anemia_events >= 30 && prep$validation_key_group_n >= 100 && prep$validation_key_group_anemia_events >= 12
clinical_gate <- (is.finite(ph[1])&&ph[1]<=-.30&&ph[3]<=-.10) || (is.finite(pa[1])&&pa[1]>=1.50&&pa[2]>=1.15) || (is.finite(pr[1])&&pr[1]>=1.35&&pr[2]>=1.05)
temporal_gate <- all(is.finite(c(dh[1],da[1],vh[1],va[1]))) && dh[1]<0 && da[1]>1 && vh[1]<0 && va[1]>1
sens_gate <- nrow(sens_df)>=8 && sum(sens_df$adverse_direction,na.rm=TRUE) >= ceiling(.75*nrow(sens_df))
nonlinear_complete <- nrow(nonlinear_df)==2 && all(is.finite(nonlinear_df$p_nonlinear))
missing_complete <- is.finite(missingness_df$odds_ratio[1])
pass <- all(c(sample_gate,retention_gate,key_gate,clinical_gate,temporal_gate,sens_gate,nonlinear_complete,missing_complete))
status <- list(candidate_code="C-N11",data_analysis_completed=TRUE,prep=prep,primary=list(adjusted_hgb_difference=list(estimate=ph[1],lower=ph[2],upper=ph[3],p=ph[4]),adjusted_anemia_pr=list(estimate=pa[1],lower=pa[2],upper=pa[3],p=pa[4]),adjusted_high_rdw_pr=list(estimate=pr[1],lower=pr[2],upper=pr[3],p=pr[4]),adjusted_mcv_difference=list(estimate=pm[1],lower=pm[2],upper=pm[3],p=pm[4])),discovery=list(hgb_difference=dh[1],anemia_pr=da[1]),validation=list(hgb_difference=vh[1],anemia_pr=va[1]),nonlinearity=nonlinear_df,sensitivity_concordant=sum(sens_df$adverse_direction,na.rm=TRUE),sensitivity_total=nrow(sens_df),gates=list(sample_gate=sample_gate,complete_case_retention_gate=retention_gate,key_group_event_gate=key_gate,clinical_magnitude_precision_gate=clinical_gate,temporal_replication_gate=temporal_gate,sensitivity_direction_gate=sens_gate,nonlinearity_completed=nonlinear_complete,missingness_mechanism_completed=missing_complete),data_gate_pass=pass,final_state=if(pass) "DATA_PASS_PENDING_FULL_TEXT_DEDUP" else "NO_GO_DATA_GATE")
write_json(status,file.path(out_dir,"c_n11_status.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
writeLines(c("# C-N11 blinded data-gate report",paste0("Actual N: ",prep$actual_n),paste0("Key group/events: ",prep$key_group_n,"/",prep$key_group_anemia_events),paste0("Adjusted Hgb difference: ",round(ph[1],3)," [",round(ph[2],3),", ",round(ph[3],3),"]"),paste0("Adjusted anemia PR: ",round(pa[1],3)," [",round(pa[2],3),", ",round(pa[3],3),"]"),paste0("Discovery/validation anemia PR: ",round(da[1],3),"/",round(va[1],3)),paste0("Sensitivity concordance: ",sum(sens_df$adverse_direction,na.rm=TRUE),"/",nrow(sens_df)),paste0("State: ",status$final_state)),file.path(out_dir,"c_n11_report.md"))
cat(toJSON(status,pretty=TRUE,auto_unbox=TRUE),"\n")
