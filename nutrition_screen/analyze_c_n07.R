options(stringsAsFactors = FALSE)
options(survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(jsonlite)
})

out_dir <- "nutrition_screen/output_c_n07"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dat <- read.csv(file.path(out_dir, "c_n07_analysis_core.csv"), na.strings = c("", "NA"))

numeric_vars <- c(
  "age", "sex", "race", "education", "pir", "weight_2day", "psu_u", "strata_u",
  "protein1_g", "protein2_g", "protein1_gkg", "protein2_gkg", "mean_protein_gkg",
  "absolute_protein_difference_gkg", "protein_cv", "energy1_kcal", "energy2_kcal",
  "mean_energy_kcal", "absolute_energy_difference_kcal", "body_weight_kg", "adjusted_weight_kg",
  "height_cm", "bmi", "protein1_gkg_adjw", "protein2_gkg_adjw", "mean_protein_gkg_adjw",
  "grip_max_kg", "grip_combined_kg", "weak_fnih", "weak_ewgsop2", "grip_bmi_ratio",
  "weak_fnih_bmi", "active_recreation", "sedentary_minutes", "diabetes", "hba1c",
  "creatinine_mg_dl", "egfr", "albumin_g_dl", "serum_total_protein_g_dl", "special_diet"
)
for (variable in intersect(numeric_vars, names(dat))) {
  dat[[variable]] <- suppressWarnings(as.numeric(dat[[variable]]))
}
dat$complete_case <- as.logical(dat$complete_case)
dat$cycle <- factor(dat$cycle, levels = c("G", "H"))
dat$period <- factor(dat$period, levels = c("discovery", "validation"))
dat$sex <- factor(dat$sex)
dat$race <- factor(dat$race)
dat$smoking <- factor(dat$smoking, levels = c("never", "former", "current"))
dat$pattern <- factor(dat$pattern, levels = c("consistent_adequate", "unstable_adequacy"))
dat$pattern_adjw <- factor(dat$pattern_adjw, levels = c("consistent_adequate", "unstable_adequacy", "below_mean_requirement"))

make_design <- function(frame) {
  svydesign(
    ids = ~psu_u,
    strata = ~strata_u,
    weights = ~weight_2day,
    nest = TRUE,
    data = frame
  )
}

safe_effect <- function(model, term, exponentiate = TRUE) {
  coefficients <- coef(model)
  if (!(term %in% names(coefficients)) || !is.finite(coefficients[[term]])) {
    return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_))
  }
  standard_error <- sqrt(diag(vcov(model)))[[term]]
  beta <- unname(coefficients[[term]])
  values <- c(beta, beta - 1.96 * standard_error, beta + 1.96 * standard_error)
  if (exponentiate) values <- exp(values)
  p_value <- tryCatch(
    unname(summary(model)$coefficients[term, ncol(summary(model)$coefficients)]),
    error = function(e) NA_real_
  )
  c(estimate = values[1], lower = values[2], upper = values[3], p = p_value)
}

full_formula <- function(outcome, frame, exposure = "pattern", include_mean_protein = TRUE) {
  rhs <- c(
    exposure,
    "age",
    "sex",
    "race",
    "bmi",
    "height_cm",
    "mean_energy_kcal",
    "active_recreation",
    "sedentary_minutes",
    "smoking",
    "diabetes",
    "egfr",
    "pir"
  )
  if (include_mean_protein) rhs <- c(rhs, "mean_protein_gkg")
  if (length(unique(frame$cycle[!is.na(frame$cycle)])) > 1) rhs <- c(rhs, "cycle")
  as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
}

crude_formula <- function(outcome, frame, exposure = "pattern") {
  rhs <- exposure
  if (length(unique(frame$cycle[!is.na(frame$cycle)])) > 1) rhs <- c(rhs, "cycle")
  as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
}

predictive_margin_difference <- function(model, frame, exposure_name, reference_level, exposed_level) {
  beta <- coef(model)
  finite <- is.finite(beta)
  beta <- beta[finite]
  covariance <- vcov(model)[finite, finite, drop = FALSE]
  terms_without_outcome <- delete.response(terms(model))

  calculate <- function(level) {
    new_data <- frame
    new_data[[exposure_name]] <- factor(level, levels = levels(frame[[exposure_name]]))
    matrix <- model.matrix(terms_without_outcome, new_data)
    matrix <- matrix[, names(beta), drop = FALSE]
    linear_predictor <- drop(matrix %*% beta)
    probability <- plogis(linear_predictor)
    weights <- frame$weight_2day
    valid <- is.finite(probability) & is.finite(weights) & weights > 0
    probability <- probability[valid]
    matrix <- matrix[valid, , drop = FALSE]
    weights <- weights[valid]
    mean_probability <- sum(weights * probability) / sum(weights)
    gradient <- colSums((weights * probability * (1 - probability)) * matrix) / sum(weights)
    list(mean = mean_probability, gradient = gradient)
  }

  reference <- calculate(reference_level)
  exposed <- calculate(exposed_level)
  difference <- exposed$mean - reference$mean
  gradient <- exposed$gradient - reference$gradient
  standard_error <- sqrt(drop(t(gradient) %*% covariance %*% gradient))
  c(
    reference_prevalence = reference$mean,
    exposed_prevalence = exposed$mean,
    difference = difference,
    lower = difference - 1.96 * standard_error,
    upper = difference + 1.96 * standard_error
  )
}

model_rows <- list()
margin_rows <- list()
missing_rows <- list()
nonlinear_rows <- list()
sensitivity_rows <- list()

append_model <- function(scope, model_name, model, term, exponentiate = TRUE) {
  effect <- safe_effect(model, term, exponentiate = exponentiate)
  model_rows[[length(model_rows) + 1]] <<- data.frame(
    scope = scope,
    model = model_name,
    term = term,
    estimate = effect[["estimate"]],
    lower = effect[["lower"]],
    upper = effect[["upper"]],
    p = effect[["p"]],
    n = nrow(model$model),
    stringsAsFactors = FALSE
  )
  effect
}

fit_scope <- function(
  frame,
  scope,
  outcome = "weak_fnih",
  exposure = "pattern",
  reference_level = "consistent_adequate",
  exposed_level = "unstable_adequacy",
  include_mean_protein = TRUE
) {
  required <- c(
    exposure,
    outcome,
    "age",
    "sex",
    "race",
    "bmi",
    "height_cm",
    "mean_energy_kcal",
    "active_recreation",
    "sedentary_minutes",
    "smoking",
    "diabetes",
    "egfr",
    "pir",
    "weight_2day",
    "psu_u",
    "strata_u"
  )
  if (include_mean_protein) required <- c(required, "mean_protein_gkg")
  frame <- frame[complete.cases(frame[, required]), , drop = FALSE]
  if (nrow(frame) < 200) return(NULL)
  frame[[exposure]] <- droplevels(factor(frame[[exposure]], levels = c(reference_level, exposed_level)))
  frame$sex <- droplevels(factor(frame$sex))
  frame$race <- droplevels(factor(frame$race))
  frame$smoking <- droplevels(factor(frame$smoking, levels = c("never", "former", "current")))
  frame$cycle <- droplevels(factor(frame$cycle, levels = c("G", "H")))
  if (sum(frame[[exposure]] == exposed_level, na.rm = TRUE) < 20 || sum(frame[[outcome]] == 1, na.rm = TRUE) < 20) return(NULL)

  design <- make_design(frame)
  crude_pr_model <- svyglm(
    crude_formula(outcome, frame, exposure),
    design = design,
    family = quasipoisson(link = "log")
  )
  adjusted_pr_model <- svyglm(
    full_formula(outcome, frame, exposure, include_mean_protein),
    design = design,
    family = quasipoisson(link = "log")
  )
  adjusted_logistic_model <- svyglm(
    full_formula(outcome, frame, exposure, include_mean_protein),
    design = design,
    family = quasibinomial(link = "logit")
  )

  term <- paste0(exposure, exposed_level)
  crude_effect <- append_model(scope, paste0("crude_", outcome, "_pr"), crude_pr_model, term)
  adjusted_effect <- append_model(scope, paste0("adjusted_", outcome, "_pr"), adjusted_pr_model, term)
  adjusted_or <- append_model(scope, paste0("adjusted_", outcome, "_or"), adjusted_logistic_model, term)
  margins <- predictive_margin_difference(
    adjusted_logistic_model,
    frame,
    exposure,
    reference_level,
    exposed_level
  )
  margin_rows[[length(margin_rows) + 1]] <<- data.frame(
    scope = scope,
    outcome = outcome,
    exposure = exposure,
    reference_level = reference_level,
    exposed_level = exposed_level,
    n = nrow(frame),
    exposed_n = sum(frame[[exposure]] == exposed_level, na.rm = TRUE),
    exposed_events = sum(frame[[outcome]][frame[[exposure]] == exposed_level] == 1, na.rm = TRUE),
    reference_prevalence = margins[["reference_prevalence"]],
    exposed_prevalence = margins[["exposed_prevalence"]],
    difference = margins[["difference"]],
    lower = margins[["lower"]],
    upper = margins[["upper"]],
    stringsAsFactors = FALSE
  )
  list(
    data = frame,
    design = design,
    crude_pr = crude_effect,
    adjusted_pr = adjusted_effect,
    adjusted_or = adjusted_or,
    margins = margins
  )
}

primary <- fit_scope(dat, "overall")
discovery <- fit_scope(dat[dat$period == "discovery", , drop = FALSE], "discovery")
validation <- fit_scope(dat[dat$period == "validation", , drop = FALSE], "validation")
if (is.null(primary)) stop("Primary model could not be estimated")

# Missingness mechanism: probability of being retained in the complete adjusted model.
missing_base <- dat[complete.cases(dat[, c("pattern", "age", "sex", "race", "cycle", "weight_2day", "psu_u", "strata_u")]), , drop = FALSE]
missing_base$retained <- as.integer(missing_base$complete_case)
missing_base$pattern <- factor(missing_base$pattern, levels = c("consistent_adequate", "unstable_adequacy"))
missing_base$sex <- factor(missing_base$sex)
missing_base$race <- factor(missing_base$race)
missing_base$cycle <- factor(missing_base$cycle)
missing_design <- make_design(missing_base)
missing_model <- svyglm(retained ~ pattern + age + sex + race + cycle, design = missing_design, family = quasibinomial())
missing_effect <- safe_effect(missing_model, "patternunstable_adequacy", exponentiate = TRUE)
missing_rows[[1]] <- data.frame(
  term = "patternunstable_adequacy",
  odds_ratio = missing_effect[["estimate"]],
  lower = missing_effect[["lower"]],
  upper = missing_effect[["upper"]],
  p = missing_effect[["p"]],
  stringsAsFactors = FALSE
)

# Restricted cubic spline for the absolute two-day protein difference.
weighted_quantile <- function(values, weights, probabilities) {
  valid <- is.finite(values) & is.finite(weights) & weights > 0
  values <- values[valid]
  weights <- weights[valid]
  order_index <- order(values)
  values <- values[order_index]
  weights <- weights[order_index]
  cumulative <- cumsum(weights) / sum(weights)
  sapply(probabilities, function(probability) values[which(cumulative >= probability)[1]])
}

rcs_basis <- function(values, knots) {
  positive_cube <- function(value) pmax(value, 0)^3
  number_of_knots <- length(knots)
  denominator <- (knots[number_of_knots] - knots[1])^2
  basis <- sapply(seq_len(number_of_knots - 2), function(index) {
    (
      positive_cube(values - knots[index]) -
      positive_cube(values - knots[number_of_knots - 1]) *
        (knots[number_of_knots] - knots[index]) /
        (knots[number_of_knots] - knots[number_of_knots - 1]) +
      positive_cube(values - knots[number_of_knots]) *
        (knots[number_of_knots - 1] - knots[index]) /
        (knots[number_of_knots] - knots[number_of_knots - 1])
    ) / denominator
  })
  colnames(basis) <- paste0("rcs", seq_len(ncol(basis)))
  basis
}

nonlinear_data <- primary$data
knots <- weighted_quantile(
  nonlinear_data$absolute_protein_difference_gkg,
  nonlinear_data$weight_2day,
  c(0.05, 0.35, 0.65, 0.95)
)
basis <- rcs_basis(nonlinear_data$absolute_protein_difference_gkg, knots)
nonlinear_data$x_linear <- nonlinear_data$absolute_protein_difference_gkg
nonlinear_data$rcs1 <- basis[, 1]
nonlinear_data$rcs2 <- basis[, 2]
nonlinear_data$cycle <- droplevels(nonlinear_data$cycle)
nonlinear_design <- make_design(nonlinear_data)
nonlinear_model <- svyglm(
  weak_fnih ~ x_linear + rcs1 + rcs2 + mean_protein_gkg + age + sex + race + bmi + height_cm +
    mean_energy_kcal + active_recreation + sedentary_minutes + smoking + diabetes + egfr + pir + cycle,
  design = nonlinear_design,
  family = quasibinomial()
)
nonlinear_test <- regTermTest(nonlinear_model, ~rcs1 + rcs2)
nonlinear_rows[[1]] <- data.frame(
  outcome = "weak_fnih",
  n = nrow(nonlinear_data),
  p_nonlinear = unname(nonlinear_test$p),
  knots_absolute_difference_gkg = paste(round(knots, 4), collapse = "|"),
  stringsAsFactors = FALSE
)

# Frozen sensitivity analyses.
run_sensitivity <- function(
  label,
  frame,
  threshold = 0.8,
  outcome = "weak_fnih",
  use_adjusted_weight = FALSE,
  include_mean_protein = TRUE
) {
  frame <- frame[!is.na(frame$protein1_g) & !is.na(frame$protein2_g), , drop = FALSE]
  if (use_adjusted_weight) {
    frame$p1 <- frame$protein1_gkg_adjw
    frame$p2 <- frame$protein2_gkg_adjw
    frame$mean_protein_gkg <- frame$mean_protein_gkg_adjw
  } else {
    frame$p1 <- frame$protein1_gkg
    frame$p2 <- frame$protein2_gkg
  }
  frame <- frame[is.finite(frame$p1) & is.finite(frame$p2) & frame$mean_protein_gkg >= threshold & frame$mean_protein_gkg <= 3, , drop = FALSE]
  frame$pattern <- ifelse(
    frame$p1 >= threshold & frame$p2 >= threshold,
    "consistent_adequate",
    "unstable_adequacy"
  )
  frame$pattern <- factor(frame$pattern, levels = c("consistent_adequate", "unstable_adequacy"))
  fit <- fit_scope(
    frame,
    paste0("sensitivity_", label),
    outcome = outcome,
    exposure = "pattern",
    include_mean_protein = include_mean_protein
  )
  if (is.null(fit)) return(NULL)
  prevalence_ratio <- fit$adjusted_pr
  margins <- fit$margins
  sensitivity_rows[[length(sensitivity_rows) + 1]] <<- data.frame(
    label = label,
    outcome = outcome,
    n = nrow(fit$data),
    key_group_n = sum(fit$data$pattern == "unstable_adequacy", na.rm = TRUE),
    key_group_events = sum(fit$data[[outcome]][fit$data$pattern == "unstable_adequacy"] == 1, na.rm = TRUE),
    adjusted_pr = prevalence_ratio[["estimate"]],
    lower = prevalence_ratio[["lower"]],
    upper = prevalence_ratio[["upper"]],
    p = prevalence_ratio[["p"]],
    adjusted_prevalence_difference = margins[["difference"]],
    pd_lower = margins[["lower"]],
    pd_upper = margins[["upper"]],
    direction_concordant = is.finite(prevalence_ratio[["estimate"]]) && prevalence_ratio[["estimate"]] > 1 &&
      is.finite(margins[["difference"]]) && margins[["difference"]] > 0,
    stringsAsFactors = FALSE
  )
}

run_sensitivity("primary_reestimated", dat)
run_sensitivity("protein_requirement_1_0_gkg", dat, threshold = 1.0)
run_sensitivity("ewgsop2_weakness", dat, outcome = "weak_ewgsop2")
run_sensitivity("fnih_bmi_adjusted_weakness", dat, outcome = "weak_fnih_bmi")
run_sensitivity(
  "plausible_energy_both_days",
  dat[dat$energy1_kcal >= 500 & dat$energy1_kcal <= 5000 & dat$energy2_kcal >= 500 & dat$energy2_kcal <= 5000, , drop = FALSE]
)
run_sensitivity("egfr_60_plus", dat[!is.na(dat$egfr) & dat$egfr >= 60, , drop = FALSE])
run_sensitivity("exclude_diabetes", dat[!is.na(dat$diabetes) & dat$diabetes == 0, , drop = FALSE])
run_sensitivity("age_65_plus", dat[!is.na(dat$age) & dat$age >= 65, , drop = FALSE])
run_sensitivity("adjusted_weight_in_obesity", dat, use_adjusted_weight = TRUE)
run_sensitivity("without_mean_protein_adjustment", dat, include_mean_protein = FALSE)
run_sensitivity("exclude_special_diet", dat[is.na(dat$special_diet) | dat$special_diet != 1, , drop = FALSE])

model_frame <- do.call(rbind, model_rows)
margin_frame <- do.call(rbind, margin_rows)
missing_frame <- do.call(rbind, missing_rows)
nonlinear_frame <- do.call(rbind, nonlinear_rows)
sensitivity_frame <- if (length(sensitivity_rows)) do.call(rbind, sensitivity_rows) else data.frame()
write.csv(model_frame, file.path(out_dir, "c_n07_models.csv"), row.names = FALSE)
write.csv(margin_frame, file.path(out_dir, "c_n07_adjusted_margins.csv"), row.names = FALSE)
write.csv(missing_frame, file.path(out_dir, "c_n07_missingness_mechanism.csv"), row.names = FALSE)
write.csv(nonlinear_frame, file.path(out_dir, "c_n07_nonlinearity.csv"), row.names = FALSE)
write.csv(sensitivity_frame, file.path(out_dir, "c_n07_sensitivities.csv"), row.names = FALSE)

prep <- fromJSON(file.path(out_dir, "c_n07_prep_status.json"))
primary_pr <- primary$adjusted_pr
primary_margin <- primary$margins
discovery_pr <- if (is.null(discovery)) rep(NA_real_, 4) else discovery$adjusted_pr
validation_pr <- if (is.null(validation)) rep(NA_real_, 4) else validation$adjusted_pr

sample_gate <- prep$actual_n >= 1500
retention_gate <- prep$complete_case_retention >= 0.80
key_group_gate <- prep$key_group_n >= 250 && prep$key_group_events >= 40
validation_precision_gate <- prep$validation_key_group_n >= 100 && prep$validation_key_group_events >= 15
clinical_gate <- (
  is.finite(primary_pr[["estimate"]]) && primary_pr[["estimate"]] >= 1.25 && primary_pr[["lower"]] >= 1.05
) || (
  is.finite(primary_margin[["difference"]]) && primary_margin[["difference"]] >= 0.04 && primary_margin[["lower"]] >= 0.01
)
temporal_gate <- all(is.finite(c(discovery_pr[["estimate"]], validation_pr[["estimate"]]))) &&
  discovery_pr[["estimate"]] > 1.10 && validation_pr[["estimate"]] > 1.10
sensitivity_gate <- nrow(sensitivity_frame) >= 8 &&
  sum(sensitivity_frame$direction_concordant, na.rm = TRUE) >= ceiling(0.75 * nrow(sensitivity_frame))
nonlinearity_completed <- nrow(nonlinear_frame) == 1 && is.finite(nonlinear_frame$p_nonlinear[1])

status <- list(
  candidate_code = "C-N07",
  data_analysis_completed = TRUE,
  prep = prep,
  primary = list(
    adjusted_prevalence_ratio = list(
      estimate = primary_pr[["estimate"]],
      lower = primary_pr[["lower"]],
      upper = primary_pr[["upper"]],
      p = primary_pr[["p"]]
    ),
    adjusted_prevalence_difference = list(
      estimate = primary_margin[["difference"]],
      lower = primary_margin[["lower"]],
      upper = primary_margin[["upper"]]
    ),
    predicted_prevalence_reference = primary_margin[["reference_prevalence"]],
    predicted_prevalence_exposed = primary_margin[["exposed_prevalence"]]
  ),
  discovery_adjusted_pr = discovery_pr[["estimate"]],
  validation_adjusted_pr = validation_pr[["estimate"]],
  nonlinearity = nonlinear_frame,
  sensitivity_concordant = if (nrow(sensitivity_frame)) sum(sensitivity_frame$direction_concordant, na.rm = TRUE) else 0,
  sensitivity_total = nrow(sensitivity_frame),
  gates = list(
    sample_gate = sample_gate,
    complete_case_retention_gate = retention_gate,
    key_group_event_gate = key_group_gate,
    validation_precision_gate = validation_precision_gate,
    clinical_magnitude_precision_gate = clinical_gate,
    temporal_replication_gate = temporal_gate,
    sensitivity_direction_gate = sensitivity_gate,
    nonlinearity_completed = nonlinearity_completed
  ),
  data_gate_pass = all(c(
    sample_gate,
    retention_gate,
    key_group_gate,
    validation_precision_gate,
    clinical_gate,
    temporal_gate,
    sensitivity_gate,
    nonlinearity_completed
  )),
  final_state = if (all(c(
    sample_gate,
    retention_gate,
    key_group_gate,
    validation_precision_gate,
    clinical_gate,
    temporal_gate,
    sensitivity_gate,
    nonlinearity_completed
  ))) "DATA_PASS_PENDING_FULL_TEXT_DEDUP" else "NO_GO_DATA_GATE"
)
write_json(status, file.path(out_dir, "c_n07_status.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
writeLines(
  c(
    "# Blind candidate C-N07 actual-data gate report",
    "",
    paste0("- Actual analytic N: ", prep$actual_n),
    paste0("- Complete-case N: ", prep$complete_case_n, " (", round(100 * prep$complete_case_retention, 1), "%)"),
    paste0("- Primary weakness events: ", prep$primary_events),
    paste0("- Key group N/events: ", prep$key_group_n, "/", prep$key_group_events),
    paste0("- Validation key group N/events: ", prep$validation_key_group_n, "/", prep$validation_key_group_events),
    paste0("- Adjusted PR: ", round(primary_pr[["estimate"]], 3), " (95% CI ", round(primary_pr[["lower"]], 3), " to ", round(primary_pr[["upper"]], 3), ")"),
    paste0("- Adjusted prevalence difference: ", round(100 * primary_margin[["difference"]], 2), " percentage points (95% CI ", round(100 * primary_margin[["lower"]], 2), " to ", round(100 * primary_margin[["upper"]], 2), ")"),
    paste0("- Discovery/validation adjusted PR: ", round(discovery_pr[["estimate"]], 3), "/", round(validation_pr[["estimate"]], 3)),
    paste0("- Concordant sensitivities: ", status$sensitivity_concordant, "/", status$sensitivity_total),
    paste0("- Final data state: ", status$final_state),
    "",
    "The scientific title remains blinded unless full-text deduplication also passes."
  ),
  file.path(out_dir, "c_n07_report.md")
)
cat(toJSON(status, pretty = TRUE, auto_unbox = TRUE), "\n")
