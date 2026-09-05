options(stringsAsFactors = FALSE)
options(survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(jsonlite)
})

out_dir <- "nutrition_stage/output_c_n12"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dat <- read.csv(file.path(out_dir, "c_n12_analysis_core.csv"), na.strings = c("", "NA"))

numeric_vars <- c(
  "age", "sex", "race", "education", "pir", "weight", "psu_u", "strata_u",
  "b12_pg_ml", "mma_umol_l", "log_b12", "log_mma", "folate_nmol_l", "folate_ng_ml",
  "hba1c", "creatinine_mg_dl", "egfr", "albumin_g_dl", "total_protein_g_dl",
  "bmi", "body_weight_kg", "height_cm", "diabetes", "left_insensate", "right_insensate",
  "pn_sites", "pn_any", "pn_two_plus"
)
for (v in intersect(numeric_vars, names(dat))) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
dat$complete_case <- as.logical(dat$complete_case)
dat$cycle <- factor(dat$cycle, levels = c("A", "B", "C"))
dat$period <- factor(dat$period, levels = c("discovery", "validation"))
dat$sex <- factor(dat$sex, levels = c(1, 2), labels = c("male", "female"))
dat$race <- factor(dat$race)
dat$smoking <- factor(dat$smoking, levels = c("never", "former", "current"))
dat$group <- factor(dat$group, levels = c("concordant_adequate", "metabolic_discordance"))

make_design <- function(d) {
  svydesign(ids = ~psu_u, strata = ~strata_u, weights = ~weight, nest = TRUE, data = d)
}

safe_effect <- function(model, term, exponentiate = FALSE) {
  cf <- coef(model)
  if (!(term %in% names(cf))) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  se <- sqrt(diag(vcov(model)))[term]
  est <- unname(cf[term])
  p <- tryCatch(unname(summary(model)$coefficients[term, ncol(summary(model)$coefficients)]), error = function(e) NA_real_)
  out <- c(est, est - 1.96 * se, est + 1.96 * se, p)
  if (exponentiate) out[1:3] <- exp(out[1:3])
  out
}

model_rows <- list()
append_model <- function(scope, model_name, model, term = "groupmetabolic_discordance", exponentiate = FALSE) {
  z <- safe_effect(model, term, exponentiate)
  model_rows[[length(model_rows) + 1]] <<- data.frame(
    candidate_code = "C-N12", scope = scope, model = model_name, term = term,
    estimate = z[1], lower = z[2], upper = z[3], p = z[4], n = nrow(model$model),
    stringsAsFactors = FALSE
  )
  z
}

fit_scope <- function(d, scope) {
  d <- d[d$complete_case & !is.na(d$group), , drop = FALSE]
  d$group <- droplevels(factor(d$group, levels = c("concordant_adequate", "metabolic_discordance")))
  if (nrow(d) < 250 || sum(d$group == "metabolic_discordance") < 30 || sum(d$pn_any == 1) < 20 || nlevels(d$group) < 2) return(NULL)
  des <- make_design(d)
  cycle_term <- if (length(unique(d$cycle)) > 1) " + cycle" else ""
  crude_rhs <- paste0("group", cycle_term)
  full_rhs <- paste0("group + age + sex + race + bmi + smoking + hba1c + egfr + folate_nmol_l + pir + albumin_g_dl", cycle_term)
  crude_pr <- svyglm(as.formula(paste("pn_any ~", crude_rhs)), design = des, family = quasipoisson(link = "log"))
  crude_rd <- svyglm(as.formula(paste("pn_any ~", crude_rhs)), design = des, family = gaussian())
  adjusted_pr <- svyglm(as.formula(paste("pn_any ~", full_rhs)), design = des, family = quasipoisson(link = "log"))
  adjusted_rd <- svyglm(as.formula(paste("pn_any ~", full_rhs)), design = des, family = gaussian())
  strict_pr <- svyglm(as.formula(paste("pn_two_plus ~", full_rhs)), design = des, family = quasipoisson(link = "log"))
  strict_rd <- svyglm(as.formula(paste("pn_two_plus ~", full_rhs)), design = des, family = gaussian())
  append_model(scope, "crude_pn_any_prevalence_ratio", crude_pr, exponentiate = TRUE)
  append_model(scope, "crude_pn_any_prevalence_difference", crude_rd)
  pr <- append_model(scope, "adjusted_pn_any_prevalence_ratio", adjusted_pr, exponentiate = TRUE)
  rd <- append_model(scope, "adjusted_pn_any_prevalence_difference", adjusted_rd)
  spr <- append_model(scope, "adjusted_pn_two_plus_prevalence_ratio", strict_pr, exponentiate = TRUE)
  srd <- append_model(scope, "adjusted_pn_two_plus_prevalence_difference", strict_rd)
  list(data = d, design = des, pr = pr, rd = rd, strict_pr = spr, strict_rd = srd)
}

primary <- fit_scope(dat, "overall")
discovery <- fit_scope(dat[dat$period == "discovery", ], "discovery")
validation <- fit_scope(dat[dat$period == "validation", ], "validation")
if (is.null(primary)) stop("Primary C-N12 model could not be estimated")

# Missingness mechanism.
miss_base <- dat[!is.na(dat$group) & !is.na(dat$age) & !is.na(dat$sex) & !is.na(dat$race), , drop = FALSE]
miss_design <- make_design(miss_base)
miss_model <- svyglm(I(complete_case) ~ group + age + sex + race + cycle, design = miss_design, family = quasibinomial())
miss_z <- safe_effect(miss_model, "groupmetabolic_discordance", exponentiate = TRUE)
missingness_df <- data.frame(
  term = "groupmetabolic_discordance", odds_ratio = miss_z[1], lower = miss_z[2], upper = miss_z[3], p = miss_z[4], stringsAsFactors = FALSE
)

# Restricted cubic spline for MMA among participants with non-low serum B12.
weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  ord <- order(x); x <- x[ord]; w <- w[ord]
  cw <- cumsum(w) / sum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}
rcs_basis <- function(x, knots) {
  pos3 <- function(z) pmax(z, 0)^3
  k <- knots
  denom <- (k[length(k)] - k[1])^2
  h <- sapply(seq_len(length(k) - 2), function(j) {
    (pos3(x - k[j])
     - pos3(x - k[length(k) - 1]) * (k[length(k)] - k[j]) / (k[length(k)] - k[length(k) - 1])
     + pos3(x - k[length(k)]) * (k[length(k) - 1] - k[j]) / (k[length(k)] - k[length(k) - 1])) / denom
  })
  colnames(h) <- paste0("rcs", seq_len(ncol(h)))
  h
}

nl <- primary$data[primary$data$b12_pg_ml >= 200, , drop = FALSE]
knots <- weighted_quantile(nl$log_mma, nl$weight, c(0.05, 0.35, 0.65, 0.95))
basis <- rcs_basis(nl$log_mma, knots)
nl$xlin <- nl$log_mma
nl$rcs1 <- basis[, 1]
nl$rcs2 <- basis[, 2]
des_nl <- make_design(nl)
nl_pr <- svyglm(
  pn_any ~ xlin + rcs1 + rcs2 + log_b12 + age + sex + race + bmi + smoking + hba1c + egfr + folate_nmol_l + pir + albumin_g_dl + cycle,
  design = des_nl, family = quasipoisson(link = "log")
)
nl_strict <- svyglm(
  pn_two_plus ~ xlin + rcs1 + rcs2 + log_b12 + age + sex + race + bmi + smoking + hba1c + egfr + folate_nmol_l + pir + albumin_g_dl + cycle,
  design = des_nl, family = quasipoisson(link = "log")
)
test_pr <- regTermTest(nl_pr, ~rcs1 + rcs2)
test_strict <- regTermTest(nl_strict, ~rcs1 + rcs2)
nonlinear_df <- rbind(
  data.frame(outcome = "pn_any", n = nrow(nl), p_nonlinear = unname(test_pr$p), knots_mma_umol_l = paste(round(exp(knots), 3), collapse = "|")),
  data.frame(outcome = "pn_two_plus", n = nrow(nl), p_nonlinear = unname(test_strict$p), knots_mma_umol_l = paste(round(exp(knots), 3), collapse = "|"))
)

# Frozen sensitivities; thresholds were specified before execution.
sensitivity_rows <- list()
run_sensitivity <- function(label, d, mma_cut = 0.26, ref_b12 = 300, exposed_b12 = 200,
                            egfr_cut = 60, folate_cut = 7, strict_outcome = FALSE,
                            extra_adjust = NULL) {
  d <- d[
    d$egfr >= egfr_cut & d$folate_nmol_l >= folate_cut & d$b12_pg_ml >= exposed_b12,
    , drop = FALSE
  ]
  d$group2 <- ifelse(
    d$mma_umol_l > mma_cut,
    "metabolic_discordance",
    ifelse(d$b12_pg_ml >= ref_b12 & d$mma_umol_l <= mma_cut, "concordant_adequate", NA_character_)
  )
  d <- d[!is.na(d$group2), , drop = FALSE]
  required <- c("age", "sex", "race", "bmi", "smoking", "hba1c", "egfr", "folate_nmol_l", "pir", "albumin_g_dl")
  d <- d[complete.cases(d[, required, drop = FALSE]), , drop = FALSE]
  d$group2 <- factor(d$group2, levels = c("concordant_adequate", "metabolic_discordance"))
  outcome <- if (strict_outcome) "pn_two_plus" else "pn_any"
  if (nrow(d) < 250 || sum(d$group2 == "metabolic_discordance") < 25 || sum(d[[outcome]] == 1) < 20 || nlevels(d$group2) < 2) return(NULL)
  rhs <- c("group2", "age", "sex", "race", "bmi", "smoking", "hba1c", "egfr", "folate_nmol_l", "pir", "albumin_g_dl")
  if (length(unique(d$cycle)) > 1) rhs <- c(rhs, "cycle")
  if (!is.null(extra_adjust)) rhs <- c(rhs, extra_adjust)
  des <- make_design(d)
  formula <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
  model_pr <- svyglm(formula, design = des, family = quasipoisson(link = "log"))
  model_rd <- svyglm(formula, design = des, family = gaussian())
  zpr <- safe_effect(model_pr, "group2metabolic_discordance", exponentiate = TRUE)
  zrd <- safe_effect(model_rd, "group2metabolic_discordance")
  sensitivity_rows[[length(sensitivity_rows) + 1]] <<- data.frame(
    label = label, outcome = outcome, n = nrow(d),
    key_group_n = sum(d$group2 == "metabolic_discordance"),
    key_group_events = sum(d[[outcome]][d$group2 == "metabolic_discordance"]),
    prevalence_ratio = zpr[1], pr_lower = zpr[2], pr_upper = zpr[3], pr_p = zpr[4],
    prevalence_difference = zrd[1], rd_lower = zrd[2], rd_upper = zrd[3], rd_p = zrd[4],
    adverse_direction = is.finite(zpr[1]) && is.finite(zrd[1]) && zpr[1] > 1 && zrd[1] > 0,
    stringsAsFactors = FALSE
  )
}

run_sensitivity("primary_reestimated", dat)
run_sensitivity("strict_neuropathy_two_plus_sites", dat, strict_outcome = TRUE)
run_sensitivity("mma_cut_0_30", dat, mma_cut = 0.30)
run_sensitivity("mma_cut_0_40", dat, mma_cut = 0.40)
run_sensitivity("reference_b12_ge_400", dat, ref_b12 = 400)
run_sensitivity("exposed_b12_ge_250", dat, exposed_b12 = 250)
run_sensitivity("normal_kidney_egfr_ge_90", dat, egfr_cut = 90)
run_sensitivity("folate_ge_10", dat, folate_cut = 10)
run_sensitivity("exclude_current_smokers", dat[dat$smoking != "current", ])
run_sensitivity("normal_glycemia_hba1c_lt_5_7", dat[dat$hba1c < 5.7, ])
run_sensitivity("non_obese_bmi_lt_30", dat[dat$bmi < 30, ])
run_sensitivity("women_only", dat[dat$sex == "female", ])
run_sensitivity("men_only", dat[dat$sex == "male", ])
run_sensitivity("age_40_59", dat[dat$age < 60, ])
run_sensitivity("age_60_plus", dat[dat$age >= 60, ])
run_sensitivity("adjust_total_protein", dat, extra_adjust = "total_protein_g_dl")
for (cy in levels(dat$cycle)) run_sensitivity(paste0("leave_out_cycle_", cy), dat[dat$cycle != cy, ])

model_df <- do.call(rbind, model_rows)
sensitivity_df <- do.call(rbind, sensitivity_rows)
write.csv(model_df, file.path(out_dir, "c_n12_models.csv"), row.names = FALSE)
write.csv(sensitivity_df, file.path(out_dir, "c_n12_sensitivities.csv"), row.names = FALSE)
write.csv(nonlinear_df, file.path(out_dir, "c_n12_nonlinearity.csv"), row.names = FALSE)
write.csv(missingness_df, file.path(out_dir, "c_n12_missingness_mechanism.csv"), row.names = FALSE)

prep <- fromJSON(file.path(out_dir, "c_n12_prep_status.json"))
pr <- primary$pr
rd <- primary$rd
spr <- primary$strict_pr
srd <- primary$strict_rd
dpr <- if (is.null(discovery)) rep(NA_real_, 4) else discovery$pr
vpr <- if (is.null(validation)) rep(NA_real_, 4) else validation$pr
drd <- if (is.null(discovery)) rep(NA_real_, 4) else discovery$rd
vrd <- if (is.null(validation)) rep(NA_real_, 4) else validation$rd

sample_gate <- prep$actual_n >= 3000
retention_gate <- prep$complete_case_retention >= 0.80
key_group_gate <- prep$key_group_n >= 300 && prep$key_group_events >= 50 && prep$validation_key_group_n >= 100 && prep$validation_key_group_events >= 20
clinical_gate <- is.finite(pr[1]) && is.finite(rd[1]) && pr[1] >= 1.35 && pr[2] >= 1.05 && rd[1] >= 0.03 && rd[2] >= 0.005
temporal_gate <- all(is.finite(c(dpr[1], vpr[1], drd[1], vrd[1]))) && dpr[1] > 1 && vpr[1] >= 1.15 && drd[1] > 0 && vrd[1] > 0
sensitivity_gate <- nrow(sensitivity_df) >= 10 && sum(sensitivity_df$adverse_direction, na.rm = TRUE) >= ceiling(0.75 * nrow(sensitivity_df))
nonlinearity_completed <- nrow(nonlinear_df) == 2 && all(is.finite(nonlinear_df$p_nonlinear))
missingness_completed <- is.finite(missingness_df$odds_ratio[1])

status <- list(
  candidate_code = "C-N12",
  data_analysis_completed = TRUE,
  prep = prep,
  primary = list(
    adjusted_pn_any_pr = list(estimate = pr[1], lower = pr[2], upper = pr[3], p = pr[4]),
    adjusted_pn_any_rd = list(estimate = rd[1], lower = rd[2], upper = rd[3], p = rd[4]),
    adjusted_pn_two_plus_pr = list(estimate = spr[1], lower = spr[2], upper = spr[3], p = spr[4]),
    adjusted_pn_two_plus_rd = list(estimate = srd[1], lower = srd[2], upper = srd[3], p = srd[4])
  ),
  discovery = list(pr = dpr[1], rd = drd[1]),
  validation = list(pr = vpr[1], rd = vrd[1]),
  nonlinearity = nonlinear_df,
  sensitivity_concordant = sum(sensitivity_df$adverse_direction, na.rm = TRUE),
  sensitivity_total = nrow(sensitivity_df),
  gates = list(
    sample_gate = sample_gate,
    complete_case_retention_gate = retention_gate,
    key_group_event_gate = key_group_gate,
    clinical_magnitude_precision_gate = clinical_gate,
    temporal_replication_gate = temporal_gate,
    sensitivity_direction_gate = sensitivity_gate,
    nonlinearity_completed = nonlinearity_completed,
    missingness_mechanism_completed = missingness_completed
  ),
  data_gate_pass = all(c(sample_gate, retention_gate, key_group_gate, clinical_gate, temporal_gate, sensitivity_gate, nonlinearity_completed, missingness_completed)),
  final_state = if (all(c(sample_gate, retention_gate, key_group_gate, clinical_gate, temporal_gate, sensitivity_gate, nonlinearity_completed, missingness_completed))) "DATA_PASS_PENDING_FULL_TEXT_DEDUP" else "NO_GO_DATA_GATE"
)
write_json(status, file.path(out_dir, "c_n12_status.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")

report <- c(
  "# Blinded nutrition candidate C-N12: actual data gate report",
  "",
  paste0("- Actual analytic cohort: ", prep$actual_n),
  paste0("- Complete-case full model: ", prep$complete_case_n, " (", round(100 * prep$complete_case_retention, 1), "%)"),
  paste0("- Key group: ", prep$key_group_n, "; primary events: ", prep$key_group_events),
  paste0("- Validation key group: ", prep$validation_key_group_n, "; primary events: ", prep$validation_key_group_events),
  paste0("- Adjusted primary PR: ", round(pr[1], 3), " (95% CI ", round(pr[2], 3), " to ", round(pr[3], 3), ")"),
  paste0("- Adjusted primary prevalence difference: ", round(100 * rd[1], 2), " percentage points (95% CI ", round(100 * rd[2], 2), " to ", round(100 * rd[3], 2), ")"),
  paste0("- Strict-outcome PR: ", round(spr[1], 3), " (95% CI ", round(spr[2], 3), " to ", round(spr[3], 3), ")"),
  paste0("- Discovery/validation PR: ", round(dpr[1], 3), " / ", round(vpr[1], 3)),
  paste0("- Concordant sensitivities: ", sum(sensitivity_df$adverse_direction, na.rm = TRUE), "/", nrow(sensitivity_df)),
  paste0("- Data gate: ", status$final_state),
  "",
  "The candidate title remains blinded unless all data gates and full-text deduplication pass."
)
writeLines(report, file.path(out_dir, "c_n12_report.md"))
cat(toJSON(status, pretty = TRUE, auto_unbox = TRUE), "\n")
