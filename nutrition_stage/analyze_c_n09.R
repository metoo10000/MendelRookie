options(stringsAsFactors = FALSE)
options(survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(jsonlite)
})

out_dir <- "nutrition_stage/output_c_n09"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dat <- read.csv(file.path(out_dir, "c_n09_analysis_core.csv"), na.strings = c("", "NA"))

numeric_vars <- c(
  "age", "sex", "race", "education", "pir", "weight", "psu_u", "strata_u",
  "plp", "four_pa", "ferritin", "hgb", "mcv", "mch", "rdw", "crp_mg_l",
  "creatinine_mg_dl", "egfr", "albumin_g_dl", "total_protein_g_dl", "alt_u_l",
  "ast_u_l", "bmi", "body_weight_kg", "height_cm", "diabetes", "known_cancer",
  "known_liver_disease", "anemia", "nonmacrocytic_anemia", "microcytosis",
  "log_plp", "log_crp", "log_ferritin"
)
for (v in intersect(numeric_vars, names(dat))) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
dat$complete_case <- as.logical(dat$complete_case)
dat$cycle <- factor(dat$cycle, levels = c("D", "E", "F"))
dat$period <- factor(dat$period, levels = c("discovery", "validation"))
dat$sex <- factor(dat$sex, levels = c(1, 2), labels = c("male", "female"))
dat$race <- factor(dat$race)
dat$smoking <- factor(dat$smoking, levels = c("never", "former", "current"))
dat$exposure <- factor(dat$exposure, levels = c("adequate", "deficient"))

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
append_model <- function(scope, model_name, model, term = "exposuredeficient", exponentiate = FALSE) {
  z <- safe_effect(model, term, exponentiate)
  model_rows[[length(model_rows) + 1]] <<- data.frame(
    candidate_code = "C-N09", scope = scope, model = model_name, term = term,
    estimate = z[1], lower = z[2], upper = z[3], p = z[4],
    n = nrow(model$model), stringsAsFactors = FALSE
  )
  z
}

full_hgb_formula <- hgb ~ exposure + age + sex + race + bmi + log_crp + log_ferritin + albumin_g_dl + smoking + diabetes + pir + cycle
full_mcv_formula <- mcv ~ exposure + age + sex + race + bmi + log_crp + log_ferritin + albumin_g_dl + smoking + diabetes + pir + cycle
full_anemia_formula <- anemia ~ exposure + age + sex + race + bmi + log_crp + log_ferritin + albumin_g_dl + smoking + diabetes + pir + cycle
full_nonmacro_formula <- nonmacrocytic_anemia ~ exposure + age + sex + race + bmi + log_crp + log_ferritin + albumin_g_dl + smoking + diabetes + pir + cycle
full_micro_formula <- microcytosis ~ exposure + age + sex + race + bmi + log_crp + log_ferritin + albumin_g_dl + smoking + diabetes + pir + cycle

fit_scope <- function(d, scope) {
  d <- d[d$complete_case & !is.na(d$exposure), , drop = FALSE]
  d$exposure <- droplevels(factor(d$exposure, levels = c("adequate", "deficient")))
  if (nrow(d) < 250 || sum(d$exposure == "deficient") < 30 || sum(d$anemia == 1) < 20) return(NULL)
  des <- make_design(d)
  crude_hgb <- svyglm(hgb ~ exposure + cycle, design = des, family = gaussian())
  crude_anemia <- svyglm(anemia ~ exposure + cycle, design = des, family = quasipoisson(link = "log"))
  adj_hgb <- svyglm(full_hgb_formula, design = des, family = gaussian())
  adj_mcv <- svyglm(full_mcv_formula, design = des, family = gaussian())
  adj_anemia <- svyglm(full_anemia_formula, design = des, family = quasipoisson(link = "log"))
  adj_nonmacro <- svyglm(full_nonmacro_formula, design = des, family = quasipoisson(link = "log"))
  adj_micro <- svyglm(full_micro_formula, design = des, family = quasipoisson(link = "log"))
  append_model(scope, "crude_hgb_difference_g_dl", crude_hgb)
  append_model(scope, "crude_anemia_prevalence_ratio", crude_anemia, exponentiate = TRUE)
  hgb <- append_model(scope, "adjusted_hgb_difference_g_dl", adj_hgb)
  mcv <- append_model(scope, "adjusted_mcv_difference_fl", adj_mcv)
  anemia <- append_model(scope, "adjusted_anemia_prevalence_ratio", adj_anemia, exponentiate = TRUE)
  nonmacro <- append_model(scope, "adjusted_nonmacrocytic_anemia_prevalence_ratio", adj_nonmacro, exponentiate = TRUE)
  micro <- append_model(scope, "adjusted_microcytosis_prevalence_ratio", adj_micro, exponentiate = TRUE)
  list(data = d, design = des, hgb = hgb, mcv = mcv, anemia = anemia, nonmacro = nonmacro, micro = micro)
}

primary <- fit_scope(dat, "overall")
discovery <- fit_scope(dat[dat$period == "discovery", ], "discovery")
validation <- fit_scope(dat[dat$period == "validation", ], "validation")
if (is.null(primary)) stop("Primary C-N09 model could not be estimated")

# Complete-case missingness mechanism.
miss_base <- dat[!is.na(dat$exposure) & !is.na(dat$age) & !is.na(dat$sex) & !is.na(dat$race), , drop = FALSE]
miss_design <- make_design(miss_base)
miss_model <- svyglm(I(complete_case) ~ exposure + age + sex + race + cycle, design = miss_design, family = quasibinomial())
missingness_rows <- list()
for (term in c("exposuredeficient")) {
  z <- safe_effect(miss_model, term, exponentiate = TRUE)
  missingness_rows[[length(missingness_rows) + 1]] <- data.frame(
    term = term, odds_ratio = z[1], lower = z[2], upper = z[3], p = z[4], stringsAsFactors = FALSE
  )
}

# Restricted cubic spline on continuous PLP with four weighted knots.
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

nl <- primary$data
knots <- weighted_quantile(nl$log_plp, nl$weight, c(0.05, 0.35, 0.65, 0.95))
b <- rcs_basis(nl$log_plp, knots)
nl$xlin <- nl$log_plp
nl$rcs1 <- b[, 1]
nl$rcs2 <- b[, 2]
des_nl <- make_design(nl)
nl_hgb <- svyglm(hgb ~ xlin + rcs1 + rcs2 + age + sex + race + bmi + log_crp + log_ferritin + albumin_g_dl + smoking + diabetes + pir + cycle, design = des_nl, family = gaussian())
nl_anemia <- svyglm(anemia ~ xlin + rcs1 + rcs2 + age + sex + race + bmi + log_crp + log_ferritin + albumin_g_dl + smoking + diabetes + pir + cycle, design = des_nl, family = quasipoisson(link = "log"))
test_hgb <- regTermTest(nl_hgb, ~rcs1 + rcs2)
test_anemia <- regTermTest(nl_anemia, ~rcs1 + rcs2)
nonlinear_df <- rbind(
  data.frame(outcome = "hemoglobin", n = nrow(nl), p_nonlinear = unname(test_hgb$p), knots_plp_nmol_l = paste(round(exp(knots), 2), collapse = "|")),
  data.frame(outcome = "anemia", n = nrow(nl), p_nonlinear = unname(test_anemia$p), knots_plp_nmol_l = paste(round(exp(knots), 2), collapse = "|"))
)

# Frozen sensitivity analyses. No data-driven cut-point search.
sensitivity_rows <- list()
run_sensitivity <- function(label, d, deficient_cut = 20, adequate_cut = 30, ferritin_low = 30, ferritin_high = 500,
                            crp_max = 10, egfr_min = 60, include_albumin = TRUE, adjust_four_pa = FALSE) {
  d <- d[
    d$ferritin >= ferritin_low & d$ferritin <= ferritin_high &
      d$crp_mg_l <= crp_max & d$egfr >= egfr_min,
    , drop = FALSE
  ]
  d$exposure2 <- ifelse(d$plp < deficient_cut, "deficient", ifelse(d$plp >= adequate_cut, "adequate", NA_character_))
  d <- d[!is.na(d$exposure2), , drop = FALSE]
  required <- c("age", "sex", "race", "bmi", "log_crp", "log_ferritin", "smoking", "diabetes", "pir")
  if (include_albumin) required <- c(required, "albumin_g_dl")
  if (adjust_four_pa) required <- c(required, "four_pa")
  d <- d[complete.cases(d[, required, drop = FALSE]), , drop = FALSE]
  d$exposure2 <- factor(d$exposure2, levels = c("adequate", "deficient"))
  if (nrow(d) < 250 || sum(d$exposure2 == "deficient") < 25 || sum(d$anemia == 1) < 20) return(NULL)
  des <- make_design(d)
  rhs <- c("exposure2", "age", "sex", "race", "bmi", "log_crp", "log_ferritin", "smoking", "diabetes", "pir", "cycle")
  if (include_albumin) rhs <- c(rhs, "albumin_g_dl")
  if (adjust_four_pa) rhs <- c(rhs, "log1p(four_pa)")
  fh <- as.formula(paste("hgb ~", paste(rhs, collapse = " + ")))
  fa <- as.formula(paste("anemia ~", paste(rhs, collapse = " + ")))
  fn <- as.formula(paste("nonmacrocytic_anemia ~", paste(rhs, collapse = " + ")))
  fm <- as.formula(paste("mcv ~", paste(rhs, collapse = " + ")))
  mh <- svyglm(fh, design = des, family = gaussian())
  ma <- svyglm(fa, design = des, family = quasipoisson(link = "log"))
  mn <- svyglm(fn, design = des, family = quasipoisson(link = "log"))
  mm <- svyglm(fm, design = des, family = gaussian())
  zh <- safe_effect(mh, "exposure2deficient")
  za <- safe_effect(ma, "exposure2deficient", exponentiate = TRUE)
  zn <- safe_effect(mn, "exposure2deficient", exponentiate = TRUE)
  zm <- safe_effect(mm, "exposure2deficient")
  sensitivity_rows[[length(sensitivity_rows) + 1]] <<- data.frame(
    label = label, n = nrow(d), key_group_n = sum(d$exposure2 == "deficient"),
    key_group_anemia = sum(d$anemia[d$exposure2 == "deficient"]),
    hgb_difference = zh[1], hgb_lower = zh[2], hgb_upper = zh[3], hgb_p = zh[4],
    anemia_pr = za[1], anemia_lower = za[2], anemia_upper = za[3], anemia_p = za[4],
    nonmacro_pr = zn[1], nonmacro_lower = zn[2], nonmacro_upper = zn[3], nonmacro_p = zn[4],
    mcv_difference = zm[1], mcv_lower = zm[2], mcv_upper = zm[3], mcv_p = zm[4],
    adverse_direction = is.finite(zh[1]) && is.finite(za[1]) && zh[1] < 0 && za[1] > 1,
    stringsAsFactors = FALSE
  )
}

run_sensitivity("primary_reestimated", dat)
run_sensitivity("strict_ferritin_50_300", dat, ferritin_low = 50, ferritin_high = 300)
run_sensitivity("low_inflammation_crp_le_3", dat, crp_max = 3)
run_sensitivity("normal_kidney_egfr_ge_90", dat, egfr_min = 90)
run_sensitivity("deficient_lt_20_adequate_ge_50", dat, adequate_cut = 50)
run_sensitivity("deficient_lt_15_adequate_ge_30", dat, deficient_cut = 15)
run_sensitivity("exclude_diabetes", dat[dat$diabetes == 0, ])
run_sensitivity("exclude_current_smokers", dat[dat$smoking != "current", ])
run_sensitivity("non_obese_bmi_lt_30", dat[dat$bmi < 30, ])
run_sensitivity("women_only", dat[dat$sex == "female", ])
run_sensitivity("men_only", dat[dat$sex == "male", ])
run_sensitivity("without_albumin_adjustment", dat, include_albumin = FALSE)
run_sensitivity("adjust_for_four_pyridoxic_acid", dat, adjust_four_pa = TRUE)
for (cy in levels(dat$cycle)) run_sensitivity(paste0("leave_out_cycle_", cy), dat[dat$cycle != cy, ])

model_df <- do.call(rbind, model_rows)
sensitivity_df <- do.call(rbind, sensitivity_rows)
missingness_df <- do.call(rbind, missingness_rows)
write.csv(model_df, file.path(out_dir, "c_n09_models.csv"), row.names = FALSE)
write.csv(sensitivity_df, file.path(out_dir, "c_n09_sensitivities.csv"), row.names = FALSE)
write.csv(nonlinear_df, file.path(out_dir, "c_n09_nonlinearity.csv"), row.names = FALSE)
write.csv(missingness_df, file.path(out_dir, "c_n09_missingness_mechanism.csv"), row.names = FALSE)

prep <- fromJSON(file.path(out_dir, "c_n09_prep_status.json"))
ph <- primary$hgb
pa <- primary$anemia
pn <- primary$nonmacro
pm <- primary$mcv
dh <- if (is.null(discovery)) rep(NA_real_, 4) else discovery$hgb
da <- if (is.null(discovery)) rep(NA_real_, 4) else discovery$anemia
vh <- if (is.null(validation)) rep(NA_real_, 4) else validation$hgb
va <- if (is.null(validation)) rep(NA_real_, 4) else validation$anemia

sample_gate <- prep$actual_n >= 3000
retention_gate <- prep$complete_case_retention >= 0.80
key_group_gate <- prep$key_group_n >= 300 && prep$key_group_anemia_events >= 40 && prep$validation_key_group_n >= 100 && prep$validation_key_group_anemia_events >= 15
clinical_gate <- (
  is.finite(ph[1]) && ph[1] <= -0.30 && ph[3] <= -0.10
) || (
  is.finite(pa[1]) && pa[1] >= 1.50 && pa[2] >= 1.15
) || (
  is.finite(pn[1]) && pn[1] >= 1.60 && pn[2] >= 1.15
)
temporal_gate <- all(is.finite(c(dh[1], da[1], vh[1], va[1]))) && dh[1] < 0 && da[1] > 1 && vh[1] < 0 && va[1] > 1
sensitivity_gate <- nrow(sensitivity_df) >= 8 && sum(sensitivity_df$adverse_direction, na.rm = TRUE) >= ceiling(0.75 * nrow(sensitivity_df))
nonlinearity_completed <- nrow(nonlinear_df) == 2 && all(is.finite(nonlinear_df$p_nonlinear))
missingness_completed <- nrow(missingness_df) >= 1 && all(is.finite(missingness_df$odds_ratio))

status <- list(
  candidate_code = "C-N09",
  data_analysis_completed = TRUE,
  prep = prep,
  primary = list(
    adjusted_hgb_difference = list(estimate = ph[1], lower = ph[2], upper = ph[3], p = ph[4]),
    adjusted_mcv_difference = list(estimate = pm[1], lower = pm[2], upper = pm[3], p = pm[4]),
    adjusted_anemia_pr = list(estimate = pa[1], lower = pa[2], upper = pa[3], p = pa[4]),
    adjusted_nonmacrocytic_anemia_pr = list(estimate = pn[1], lower = pn[2], upper = pn[3], p = pn[4])
  ),
  discovery = list(hgb_difference = dh[1], anemia_pr = da[1]),
  validation = list(hgb_difference = vh[1], anemia_pr = va[1]),
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
write_json(status, file.path(out_dir, "c_n09_status.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")

report <- c(
  "# Blinded nutrition candidate C-N09: actual data gate report",
  "",
  paste0("- Actual analytic cohort: ", prep$actual_n),
  paste0("- Complete-case full model: ", prep$complete_case_n, " (", round(100 * prep$complete_case_retention, 1), "%)"),
  paste0("- Key group: ", prep$key_group_n, "; anemia events: ", prep$key_group_anemia_events),
  paste0("- Validation key group: ", prep$validation_key_group_n, "; anemia events: ", prep$validation_key_group_anemia_events),
  paste0("- Adjusted hemoglobin difference: ", round(ph[1], 3), " g/dL (95% CI ", round(ph[2], 3), " to ", round(ph[3], 3), ")"),
  paste0("- Adjusted MCV difference: ", round(pm[1], 3), " fL (95% CI ", round(pm[2], 3), " to ", round(pm[3], 3), ")"),
  paste0("- Adjusted anemia PR: ", round(pa[1], 3), " (95% CI ", round(pa[2], 3), " to ", round(pa[3], 3), ")"),
  paste0("- Adjusted nonmacrocytic anemia PR: ", round(pn[1], 3), " (95% CI ", round(pn[2], 3), " to ", round(pn[3], 3), ")"),
  paste0("- Discovery/validation hemoglobin differences: ", round(dh[1], 3), " / ", round(vh[1], 3)),
  paste0("- Discovery/validation anemia PRs: ", round(da[1], 3), " / ", round(va[1], 3)),
  paste0("- Concordant sensitivities: ", sum(sensitivity_df$adverse_direction, na.rm = TRUE), "/", nrow(sensitivity_df)),
  paste0("- Data gate: ", status$final_state),
  "",
  "The candidate title remains blinded unless all data gates and full-text deduplication pass."
)
writeLines(report, file.path(out_dir, "c_n09_report.md"))
cat(toJSON(status, pretty = TRUE, auto_unbox = TRUE), "\n")
