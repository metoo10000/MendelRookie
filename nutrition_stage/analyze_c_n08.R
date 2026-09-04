options(stringsAsFactors = FALSE)
options(survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(survival)
  library(jsonlite)
})

out_dir <- "nutrition_stage/output_c_n08"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dat <- read.csv(file.path(out_dir, "c_n08_analysis_core.csv"), na.strings = c("", "NA"))

numeric_vars <- c(
  "age", "sex", "race", "education", "pir", "weight_2day", "psu_u", "strata_u",
  "protein1_g", "protein2_g", "protein1_gkg", "protein2_gkg", "mean_protein_gkg",
  "minimum_protein_gkg", "absolute_protein_difference_gkg", "energy1_kcal", "energy2_kcal",
  "mean_energy_kcal", "absolute_energy_difference_kcal", "body_weight_kg", "adjusted_weight_kg",
  "height_cm", "bmi", "protein1_gkg_adjw", "protein2_gkg_adjw", "mean_protein_gkg_adjw",
  "active_recreation", "sedentary_minutes", "diabetes", "hba1c", "creatinine_mg_dl", "egfr",
  "cvd", "cancer", "albumin_g_dl", "special_diet", "death", "follow_up_months", "follow_up_years"
)
for (v in intersect(numeric_vars, names(dat))) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
dat$complete_case <- as.logical(dat$complete_case)
dat$cycle <- factor(dat$cycle, levels = c("E", "F", "G", "H", "I", "J"))
dat$period <- factor(dat$period, levels = c("discovery", "validation"))
dat$sex <- factor(dat$sex)
dat$race <- factor(dat$race)
dat$smoking <- factor(dat$smoking, levels = c("never", "former", "current"))
dat$pattern <- factor(dat$pattern, levels = c("consistent_rda", "episodic_sub_ear"))
dat$death5 <- ifelse(dat$death == 1 & dat$follow_up_years <= 5, 1,
                     ifelse(dat$follow_up_years >= 5, 0, NA))

make_design <- function(d) {
  svydesign(ids = ~psu_u, strata = ~strata_u, weights = ~weight_2day, nest = TRUE, data = d)
}

safe_effect <- function(model, term, exponentiate = TRUE) {
  b <- coef(model)
  if (!(term %in% names(b)) || !is.finite(b[[term]])) {
    return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_))
  }
  se <- sqrt(diag(vcov(model)))[[term]]
  values <- c(unname(b[[term]]), unname(b[[term]]) - 1.96 * se, unname(b[[term]]) + 1.96 * se)
  if (exponentiate) values <- exp(values)
  p <- tryCatch(unname(summary(model)$coefficients[term, ncol(summary(model)$coefficients)]), error = function(e) NA_real_)
  c(estimate = values[1], lower = values[2], upper = values[3], p = p)
}

full_rhs <- function(d, include_mean_protein = TRUE, include_albumin = FALSE) {
  rhs <- c("pattern", "age", "sex", "race", "bmi", "height_cm", "mean_energy_kcal",
           "active_recreation", "sedentary_minutes", "smoking", "diabetes", "egfr",
           "cvd", "cancer", "pir")
  if (include_mean_protein) rhs <- c(rhs, "mean_protein_gkg")
  if (include_albumin) rhs <- c(rhs, "albumin_g_dl")
  if (length(unique(d$cycle[!is.na(d$cycle)])) > 1) rhs <- c(rhs, "cycle")
  rhs
}

cox_formula <- function(d, include_mean_protein = TRUE, include_albumin = FALSE) {
  as.formula(paste("Surv(follow_up_years, death) ~", paste(full_rhs(d, include_mean_protein, include_albumin), collapse = " + ")))
}

logistic_formula <- function(d, include_mean_protein = TRUE, include_albumin = FALSE) {
  as.formula(paste("death5 ~", paste(full_rhs(d, include_mean_protein, include_albumin), collapse = " + ")))
}

predictive_margin_difference <- function(model, d) {
  beta <- coef(model)
  keep <- is.finite(beta)
  beta <- beta[keep]
  V <- vcov(model)[keep, keep, drop = FALSE]
  tt <- delete.response(terms(model))
  evaluate <- function(level) {
    nd <- d
    nd$pattern <- factor(level, levels = levels(d$pattern))
    X <- model.matrix(tt, nd)
    X <- X[, names(beta), drop = FALSE]
    eta <- drop(X %*% beta)
    p <- plogis(eta)
    w <- d$weight_2day
    ok <- is.finite(p) & is.finite(w) & w > 0
    p <- p[ok]; X <- X[ok, , drop = FALSE]; w <- w[ok]
    mean_p <- sum(w * p) / sum(w)
    grad <- colSums((w * p * (1 - p)) * X) / sum(w)
    list(mean = mean_p, grad = grad)
  }
  r <- evaluate("consistent_rda")
  e <- evaluate("episodic_sub_ear")
  diff <- e$mean - r$mean
  grad <- e$grad - r$grad
  se <- sqrt(drop(t(grad) %*% V %*% grad))
  c(reference = r$mean, exposed = e$mean, difference = diff,
    lower = diff - 1.96 * se, upper = diff + 1.96 * se)
}

model_rows <- list()
margin_rows <- list()
missing_rows <- list()
nonlinear_rows <- list()
sensitivity_rows <- list()

append_model <- function(scope, name, model, n, exponentiate = TRUE) {
  z <- safe_effect(model, "patternepisodic_sub_ear", exponentiate)
  model_rows[[length(model_rows) + 1]] <<- data.frame(
    scope = scope, model = name, term = "patternepisodic_sub_ear", n = n,
    estimate = z[["estimate"]], lower = z[["lower"]], upper = z[["upper"]], p = z[["p"]],
    stringsAsFactors = FALSE
  )
  z
}

fit_scope <- function(d, scope, include_mean_protein = TRUE, include_albumin = FALSE) {
  needed <- c("follow_up_years", "death", full_rhs(d, include_mean_protein, include_albumin),
              "weight_2day", "psu_u", "strata_u")
  d <- d[complete.cases(d[, unique(needed)]), , drop = FALSE]
  if (nrow(d) < 300 || sum(d$pattern == "episodic_sub_ear") < 30 || sum(d$death == 1) < 30) return(NULL)
  d$pattern <- droplevels(factor(d$pattern, levels = c("consistent_rda", "episodic_sub_ear")))
  d$sex <- droplevels(factor(d$sex)); d$race <- droplevels(factor(d$race))
  d$smoking <- droplevels(factor(d$smoking, levels = c("never", "former", "current")))
  d$cycle <- droplevels(factor(d$cycle, levels = c("E", "F", "G", "H", "I", "J")))
  design <- make_design(d)
  crude_formula <- if (length(unique(d$cycle)) > 1) {
    Surv(follow_up_years, death) ~ pattern + cycle
  } else {
    Surv(follow_up_years, death) ~ pattern
  }
  crude <- svycoxph(crude_formula, design = design)
  adjusted <- svycoxph(cox_formula(d, include_mean_protein, include_albumin), design = design)
  crude_effect <- append_model(scope, "crude_survey_cox", crude, nrow(d))
  adjusted_effect <- append_model(scope, "fully_adjusted_survey_cox", adjusted, nrow(d))

  d5 <- d[!is.na(d$death5), , drop = FALSE]
  margin <- rep(NA_real_, 5); names(margin) <- c("reference", "exposed", "difference", "lower", "upper")
  if (nrow(d5) >= 300 && sum(d5$death5 == 1) >= 30 && sum(d5$pattern == "episodic_sub_ear") >= 30) {
    design5 <- make_design(d5)
    logistic <- svyglm(logistic_formula(d5, include_mean_protein, include_albumin), design = design5, family = quasibinomial())
    append_model(scope, "fully_adjusted_5y_logistic_or", logistic, nrow(d5))
    margin <- predictive_margin_difference(logistic, d5)
    margin_rows[[length(margin_rows) + 1]] <<- data.frame(
      scope = scope, n = nrow(d5), deaths_5y = sum(d5$death5 == 1),
      reference_risk = margin[["reference"]], exposed_risk = margin[["exposed"]],
      risk_difference = margin[["difference"]], lower = margin[["lower"]], upper = margin[["upper"]],
      stringsAsFactors = FALSE
    )
  }
  list(data = d, adjusted = adjusted_effect, crude = crude_effect, margin5 = margin)
}

primary <- fit_scope(dat, "overall")
discovery <- fit_scope(dat[dat$period == "discovery", , drop = FALSE], "discovery")
validation <- fit_scope(dat[dat$period == "validation", , drop = FALSE], "validation")
if (is.null(primary)) stop("Primary model could not be estimated")

# Missingness mechanism audit.
miss <- dat[complete.cases(dat[, c("pattern", "age", "sex", "race", "cycle", "weight_2day", "psu_u", "strata_u")]), , drop = FALSE]
miss$retained <- as.integer(miss$complete_case)
miss$pattern <- factor(miss$pattern, levels = c("consistent_rda", "episodic_sub_ear"))
miss$sex <- factor(miss$sex); miss$race <- factor(miss$race); miss$cycle <- factor(miss$cycle)
miss_design <- make_design(miss)
miss_model <- svyglm(retained ~ pattern + age + sex + race + cycle, design = miss_design, family = quasibinomial())
missing_effect <- safe_effect(miss_model, "patternepisodic_sub_ear", TRUE)
missing_rows[[1]] <- data.frame(
  term = "patternepisodic_sub_ear", odds_ratio = missing_effect[["estimate"]],
  lower = missing_effect[["lower"]], upper = missing_effect[["upper"]], p = missing_effect[["p"]],
  stringsAsFactors = FALSE
)

# Four-knot restricted cubic spline for absolute day-to-day protein difference.
weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  ord <- order(x); x <- x[ord]; w <- w[ord]
  cw <- cumsum(w) / sum(w)
  sapply(probs, function(prob) x[which(cw >= prob)[1]])
}
rcs_basis <- function(x, knots) {
  pos3 <- function(z) pmax(z, 0)^3
  K <- length(knots); den <- (knots[K] - knots[1])^2
  ans <- sapply(seq_len(K - 2), function(j) {
    (pos3(x - knots[j]) - pos3(x - knots[K - 1]) * (knots[K] - knots[j]) / (knots[K] - knots[K - 1]) +
       pos3(x - knots[K]) * (knots[K - 1] - knots[j]) / (knots[K] - knots[K - 1])) / den
  })
  colnames(ans) <- paste0("rcs", seq_len(ncol(ans)))
  ans
}
nl <- primary$data
knots <- weighted_quantile(nl$absolute_protein_difference_gkg, nl$weight_2day, c(0.05, 0.35, 0.65, 0.95))
basis <- rcs_basis(nl$absolute_protein_difference_gkg, knots)
nl$xlin <- nl$absolute_protein_difference_gkg; nl$rcs1 <- basis[, 1]; nl$rcs2 <- basis[, 2]
nl_design <- make_design(nl)
nl_formula <- Surv(follow_up_years, death) ~ xlin + rcs1 + rcs2 + mean_protein_gkg + age + sex + race + bmi +
  height_cm + mean_energy_kcal + active_recreation + sedentary_minutes + smoking + diabetes + egfr + cvd + cancer + pir + cycle
nl_model <- svycoxph(nl_formula, design = nl_design)
nl_test <- regTermTest(nl_model, ~rcs1 + rcs2)
nonlinear_rows[[1]] <- data.frame(
  outcome = "all_cause_mortality", n = nrow(nl), p_nonlinear = unname(nl_test$p),
  knots_absolute_difference_gkg = paste(round(knots, 4), collapse = "|"), stringsAsFactors = FALSE
)

run_sensitivity <- function(label, d, include_mean_protein = TRUE, include_albumin = FALSE, adjusted_weight = FALSE) {
  if (adjusted_weight) {
    d$mean_protein_gkg <- d$mean_protein_gkg_adjw
    d$min_adj <- pmin(d$protein1_gkg_adjw, d$protein2_gkg_adjw)
    d <- d[d$mean_protein_gkg >= 0.8 & d$mean_protein_gkg <= 3, , drop = FALSE]
    d$pattern <- ifelse(d$protein1_gkg_adjw >= 0.8 & d$protein2_gkg_adjw >= 0.8, "consistent_rda",
                        ifelse(d$min_adj < 0.66, "episodic_sub_ear", NA))
    d <- d[!is.na(d$pattern), , drop = FALSE]
    d$pattern <- factor(d$pattern, levels = c("consistent_rda", "episodic_sub_ear"))
  }
  fit <- fit_scope(d, paste0("sensitivity_", label), include_mean_protein, include_albumin)
  if (is.null(fit)) return(NULL)
  z <- fit$adjusted; m <- fit$margin5
  sensitivity_rows[[length(sensitivity_rows) + 1]] <<- data.frame(
    label = label, n = nrow(fit$data), exposed_n = sum(fit$data$pattern == "episodic_sub_ear"),
    exposed_deaths = sum(fit$data$death[fit$data$pattern == "episodic_sub_ear"] == 1),
    adjusted_hr = z[["estimate"]], lower = z[["lower"]], upper = z[["upper"]], p = z[["p"]],
    risk_difference_5y = m[["difference"]], rd_lower = m[["lower"]], rd_upper = m[["upper"]],
    direction_concordant = is.finite(z[["estimate"]]) && z[["estimate"]] > 1,
    stringsAsFactors = FALSE
  )
}

run_sensitivity("primary_reestimated", dat)
run_sensitivity("exclude_deaths_first_24_months", dat[!(dat$death == 1 & dat$follow_up_months < 24), , drop = FALSE])
run_sensitivity("age_60_plus", dat[dat$age >= 60, , drop = FALSE])
run_sensitivity("egfr_60_plus", dat[dat$egfr >= 60, , drop = FALSE])
run_sensitivity("exclude_baseline_cancer", dat[dat$cancer == 0, , drop = FALSE])
run_sensitivity("exclude_baseline_cvd", dat[dat$cvd == 0, , drop = FALSE])
run_sensitivity("plausible_energy_both_days", dat[dat$energy1_kcal >= 500 & dat$energy1_kcal <= 5000 & dat$energy2_kcal >= 500 & dat$energy2_kcal <= 5000, , drop = FALSE])
run_sensitivity("adjusted_weight_in_obesity", dat, adjusted_weight = TRUE)
run_sensitivity("without_mean_protein_adjustment", dat, include_mean_protein = FALSE)
run_sensitivity("exclude_special_diet", dat[is.na(dat$special_diet) | dat$special_diet != 1, , drop = FALSE])
if (mean(!is.na(dat$albumin_g_dl)) >= 0.80) run_sensitivity("additional_albumin_adjustment", dat, include_albumin = TRUE)

model_df <- do.call(rbind, model_rows)
margin_df <- if (length(margin_rows)) do.call(rbind, margin_rows) else data.frame()
missing_df <- do.call(rbind, missing_rows)
nonlinear_df <- do.call(rbind, nonlinear_rows)
sensitivity_df <- if (length(sensitivity_rows)) do.call(rbind, sensitivity_rows) else data.frame()
write.csv(model_df, file.path(out_dir, "c_n08_models.csv"), row.names = FALSE)
write.csv(margin_df, file.path(out_dir, "c_n08_five_year_margins.csv"), row.names = FALSE)
write.csv(missing_df, file.path(out_dir, "c_n08_missingness_mechanism.csv"), row.names = FALSE)
write.csv(nonlinear_df, file.path(out_dir, "c_n08_nonlinearity.csv"), row.names = FALSE)
write.csv(sensitivity_df, file.path(out_dir, "c_n08_sensitivities.csv"), row.names = FALSE)

prep <- fromJSON(file.path(out_dir, "c_n08_prep_status.json"))
phr <- primary$adjusted; pm <- primary$margin5
dhr <- if (is.null(discovery)) rep(NA_real_, 4) else discovery$adjusted
vhr <- if (is.null(validation)) rep(NA_real_, 4) else validation$adjusted

sample_gate <- prep$actual_n >= 4000
retention_gate <- prep$complete_case_retention >= 0.80
key_group_gate <- prep$key_group_n >= 500 && prep$key_group_deaths >= 100
validation_gate <- prep$validation_key_group_n >= 200 && prep$validation_key_group_deaths >= 30
clinical_gate <- (is.finite(phr[["estimate"]]) && phr[["estimate"]] >= 1.20 && phr[["lower"]] >= 1.05) ||
  (is.finite(pm[["difference"]]) && pm[["difference"]] >= 0.03 && pm[["lower"]] >= 0.005)
temporal_gate <- all(is.finite(c(dhr[["estimate"]], vhr[["estimate"]]))) && dhr[["estimate"]] > 1.10 && vhr[["estimate"]] > 1.10
sensitivity_gate <- nrow(sensitivity_df) >= 8 && sum(sensitivity_df$direction_concordant, na.rm = TRUE) >= ceiling(0.75 * nrow(sensitivity_df))
nonlinearity_completed <- nrow(nonlinear_df) == 1 && is.finite(nonlinear_df$p_nonlinear[1])

status <- list(
  candidate_code = "C-N08", data_analysis_completed = TRUE, prep = prep,
  primary = list(
    adjusted_hazard_ratio = list(estimate = phr[["estimate"]], lower = phr[["lower"]], upper = phr[["upper"]], p = phr[["p"]]),
    adjusted_five_year_risk_difference = list(estimate = pm[["difference"]], lower = pm[["lower"]], upper = pm[["upper"]]),
    five_year_reference_risk = pm[["reference"]], five_year_exposed_risk = pm[["exposed"]]
  ),
  discovery_adjusted_hr = dhr[["estimate"]], validation_adjusted_hr = vhr[["estimate"]],
  sensitivity_concordant = if (nrow(sensitivity_df)) sum(sensitivity_df$direction_concordant, na.rm = TRUE) else 0,
  sensitivity_total = nrow(sensitivity_df), nonlinearity = nonlinear_df,
  gates = list(sample_gate = sample_gate, complete_case_retention_gate = retention_gate,
               key_group_event_gate = key_group_gate, validation_precision_gate = validation_gate,
               clinical_magnitude_precision_gate = clinical_gate, temporal_replication_gate = temporal_gate,
               sensitivity_direction_gate = sensitivity_gate, nonlinearity_completed = nonlinearity_completed),
  data_gate_pass = all(c(sample_gate, retention_gate, key_group_gate, validation_gate, clinical_gate, temporal_gate, sensitivity_gate, nonlinearity_completed)),
  final_state = if (all(c(sample_gate, retention_gate, key_group_gate, validation_gate, clinical_gate, temporal_gate, sensitivity_gate, nonlinearity_completed))) "DATA_PASS_PENDING_FULL_TEXT_DEDUP" else "NO_GO_DATA_GATE"
)
write_json(status, file.path(out_dir, "c_n08_status.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
writeLines(c(
  "# Blind candidate C-N08 actual-data gate report", "",
  paste0("- Actual analytic N: ", prep$actual_n),
  paste0("- Complete-case N: ", prep$complete_case_n, " (", round(100 * prep$complete_case_retention, 1), "%)"),
  paste0("- Deaths: ", prep$deaths),
  paste0("- Key group N/deaths: ", prep$key_group_n, "/", prep$key_group_deaths),
  paste0("- Validation key group N/deaths: ", prep$validation_key_group_n, "/", prep$validation_key_group_deaths),
  paste0("- Adjusted HR: ", round(phr[["estimate"]], 3), " (95% CI ", round(phr[["lower"]], 3), " to ", round(phr[["upper"]], 3), ")"),
  paste0("- Adjusted 5-year risk difference: ", round(100 * pm[["difference"]], 2), " percentage points (95% CI ", round(100 * pm[["lower"]], 2), " to ", round(100 * pm[["upper"]], 2), ")"),
  paste0("- Discovery/validation HR: ", round(dhr[["estimate"]], 3), "/", round(vhr[["estimate"]], 3)),
  paste0("- Concordant sensitivities: ", status$sensitivity_concordant, "/", status$sensitivity_total),
  paste0("- Final data state: ", status$final_state), "",
  "The scientific title remains blinded unless full-text deduplication also passes."
), file.path(out_dir, "c_n08_report.md"))
cat(toJSON(status, pretty = TRUE, auto_unbox = TRUE), "\n")
