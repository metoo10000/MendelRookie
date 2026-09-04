options(stringsAsFactors = FALSE)
options(survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(jsonlite)
})

out_dir <- "nutrition_screen/output_c_n06"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dat <- read.csv(file.path(out_dir, "c_n06_analysis_core.csv"), na.strings = c("", "NA"))

numeric_vars <- c(
  "age", "sex", "race", "weight", "psu_u", "strata_u", "pir", "bmi", "diabetes", "hba1c",
  "b12_pg_ml", "mma_umol_l", "log_mma", "folate_nmol_l", "homocysteine", "creatinine_mg_dl",
  "egfr", "pn_sites", "pn_any", "pn_two_plus"
)
for (v in intersect(numeric_vars, names(dat))) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
dat$complete_case <- as.logical(dat$complete_case)
dat$cycle <- factor(dat$cycle, levels = c("A", "B", "C"))
dat$period <- factor(dat$period, levels = c("discovery", "validation"))
dat$sex <- factor(dat$sex)
dat$race <- factor(dat$race)
dat$smoking <- factor(dat$smoking, levels = c("never", "former", "current"))
dat$group <- factor(dat$group, levels = c("adequate", "discordant", "low_b12"))

make_design <- function(d) {
  svydesign(ids = ~psu_u, strata = ~strata_u, weights = ~weight, nest = TRUE, data = d)
}

safe_coef <- function(model, term, exponentiate = FALSE) {
  cf <- coef(model)
  vv <- vcov(model)
  if (!(term %in% names(cf)) || !is.finite(cf[[term]])) {
    return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_))
  }
  se <- sqrt(diag(vv))[[term]]
  est <- unname(cf[[term]])
  z <- c(est, est - 1.96 * se, est + 1.96 * se)
  if (exponentiate) z <- exp(z)
  p <- tryCatch(unname(summary(model)$coefficients[term, ncol(summary(model)$coefficients)]), error = function(e) NA_real_)
  c(estimate = z[1], lower = z[2], upper = z[3], p = p)
}

formula_for <- function(outcome, d, add_folate = FALSE) {
  rhs <- c("group", "age", "sex", "race", "bmi", "pir", "smoking", "diabetes", "hba1c", "egfr")
  if (length(unique(d$cycle[!is.na(d$cycle)])) > 1) rhs <- c(rhs, "cycle")
  if (add_folate) rhs <- c(rhs, "log_folate")
  as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
}

crude_formula_for <- function(outcome, d) {
  rhs <- "group"
  if (length(unique(d$cycle[!is.na(d$cycle)])) > 1) rhs <- c(rhs, "cycle")
  as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
}

model_rows <- list()
margin_rows <- list()
nonlinear_rows <- list()
missing_rows <- list()
sensitivity_rows <- list()

append_model <- function(scope, model_name, model, term = "groupdiscordant", exponentiate = TRUE) {
  z <- safe_coef(model, term, exponentiate = exponentiate)
  model_rows[[length(model_rows) + 1]] <<- data.frame(
    scope = scope,
    model = model_name,
    term = term,
    estimate = z[["estimate"]],
    lower = z[["lower"]],
    upper = z[["upper"]],
    p = z[["p"]],
    n = nrow(model$model),
    stringsAsFactors = FALSE
  )
  z
}

predictive_margin_difference <- function(model, d, exposed = "discordant", reference = "adequate") {
  beta <- coef(model)
  keep <- is.finite(beta)
  beta <- beta[keep]
  V <- vcov(model)[keep, keep, drop = FALSE]
  terms_no_y <- delete.response(terms(model))
  build <- function(level) {
    nd <- d
    nd$group <- factor(level, levels = levels(d$group))
    X <- model.matrix(terms_no_y, nd)
    X <- X[, names(beta), drop = FALSE]
    eta <- drop(X %*% beta)
    p <- plogis(eta)
    w <- d$weight
    ok <- is.finite(p) & is.finite(w) & w > 0
    p <- p[ok]
    X <- X[ok, , drop = FALSE]
    w <- w[ok]
    mean_p <- sum(w * p) / sum(w)
    grad <- colSums((w * p * (1 - p)) * X) / sum(w)
    list(mean = mean_p, grad = grad)
  }
  e <- build(exposed)
  r <- build(reference)
  diff <- e$mean - r$mean
  g <- e$grad - r$grad
  se <- sqrt(drop(t(g) %*% V %*% g))
  c(
    exposed_prevalence = e$mean,
    reference_prevalence = r$mean,
    difference = diff,
    lower = diff - 1.96 * se,
    upper = diff + 1.96 * se
  )
}

fit_scope <- function(d, scope, outcome = "pn_any", add_folate = FALSE) {
  needed <- c("group", "age", "sex", "race", "bmi", "pir", "smoking", "diabetes", "hba1c", "egfr", outcome, "weight", "psu_u", "strata_u")
  if (add_folate) needed <- c(needed, "folate_nmol_l")
  d <- d[complete.cases(d[, needed]), , drop = FALSE]
  if (nrow(d) < 300 || sum(d$group == "discordant", na.rm = TRUE) < 20 || sum(d[[outcome]] == 1, na.rm = TRUE) < 20) return(NULL)
  d$group <- droplevels(factor(d$group, levels = c("adequate", "discordant", "low_b12")))
  d$sex <- droplevels(factor(d$sex))
  d$race <- droplevels(factor(d$race))
  d$smoking <- droplevels(factor(d$smoking, levels = c("never", "former", "current")))
  d$cycle <- droplevels(factor(d$cycle, levels = c("A", "B", "C")))
  if (add_folate) d$log_folate <- log(pmax(d$folate_nmol_l, 0.1))
  des <- make_design(d)
  crude_pr <- svyglm(crude_formula_for(outcome, d), design = des, family = quasipoisson(link = "log"))
  adjusted_pr <- svyglm(formula_for(outcome, d, add_folate = add_folate), design = des, family = quasipoisson(link = "log"))
  adjusted_logit <- svyglm(formula_for(outcome, d, add_folate = add_folate), design = des, family = quasibinomial(link = "logit"))
  crude <- append_model(scope, paste0("crude_", outcome, "_prevalence_ratio"), crude_pr)
  adjusted <- append_model(scope, paste0("adjusted_", outcome, "_prevalence_ratio"), adjusted_pr)
  adjusted_or <- append_model(scope, paste0("adjusted_", outcome, "_odds_ratio"), adjusted_logit)
  margins <- predictive_margin_difference(adjusted_logit, d)
  margin_rows[[length(margin_rows) + 1]] <<- data.frame(
    scope = scope,
    outcome = outcome,
    exposed_prevalence = margins[["exposed_prevalence"]],
    reference_prevalence = margins[["reference_prevalence"]],
    difference = margins[["difference"]],
    lower = margins[["lower"]],
    upper = margins[["upper"]],
    n = nrow(d),
    key_group_n = sum(d$group == "discordant"),
    key_group_events = sum(d[[outcome]][d$group == "discordant"] == 1),
    stringsAsFactors = FALSE
  )
  list(data = d, design = des, crude = crude, adjusted_pr = adjusted, adjusted_or = adjusted_or, margins = margins)
}

primary <- fit_scope(dat, "overall")
discovery <- fit_scope(dat[dat$period == "discovery", , drop = FALSE], "discovery")
validation <- fit_scope(dat[dat$period == "validation", , drop = FALSE], "validation")
if (is.null(primary)) stop("Primary model could not be estimated")

# Missingness mechanism audit for full-model completeness.
miss_base <- dat[complete.cases(dat[, c("group", "age", "sex", "race", "weight", "psu_u", "strata_u")]), , drop = FALSE]
miss_base$complete_case_num <- as.integer(miss_base$complete_case)
miss_base$group <- factor(miss_base$group, levels = c("adequate", "discordant", "low_b12"))
miss_base$sex <- factor(miss_base$sex)
miss_base$race <- factor(miss_base$race)
miss_base$cycle <- factor(miss_base$cycle)
miss_design <- make_design(miss_base)
miss_model <- svyglm(complete_case_num ~ group + age + sex + race + cycle, design = miss_design, family = quasibinomial())
for (term in c("groupdiscordant", "grouplow_b12")) {
  z <- safe_coef(miss_model, term, exponentiate = TRUE)
  missing_rows[[length(missing_rows) + 1]] <- data.frame(
    term = term,
    odds_ratio = z[["estimate"]],
    lower = z[["lower"]],
    upper = z[["upper"]],
    p = z[["p"]],
    stringsAsFactors = FALSE
  )
}

# Four-knot restricted cubic spline of log(MMA) among serum-B12-normal participants.
weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}

rcs_basis <- function(x, knots) {
  pos3 <- function(z) pmax(z, 0)^3
  K <- length(knots)
  den <- (knots[K] - knots[1])^2
  h <- sapply(seq_len(K - 2), function(j) {
    (pos3(x - knots[j]) -
       pos3(x - knots[K - 1]) * (knots[K] - knots[j]) / (knots[K] - knots[K - 1]) +
       pos3(x - knots[K]) * (knots[K - 1] - knots[j]) / (knots[K] - knots[K - 1])) / den
  })
  colnames(h) <- paste0("rcs", seq_len(ncol(h)))
  h
}

nl <- primary$data[primary$data$b12_pg_ml >= 300, , drop = FALSE]
knots <- weighted_quantile(nl$log_mma, nl$weight, c(0.05, 0.35, 0.65, 0.95))
basis <- rcs_basis(nl$log_mma, knots)
nl$xlin <- nl$log_mma
nl$rcs1 <- basis[, 1]
nl$rcs2 <- basis[, 2]
nl$cycle <- droplevels(nl$cycle)
des_nl <- make_design(nl)
nl_formula <- pn_any ~ xlin + rcs1 + rcs2 + age + sex + race + bmi + pir + smoking + diabetes + hba1c + egfr + cycle
nl_model <- svyglm(nl_formula, design = des_nl, family = quasibinomial())
nl_test <- regTermTest(nl_model, ~rcs1 + rcs2)
nonlinear_rows[[1]] <- data.frame(
  outcome = "pn_any",
  n = nrow(nl),
  p_nonlinear = unname(nl_test$p),
  knots_mma_umol_l = paste(round(exp(knots), 4), collapse = "|"),
  stringsAsFactors = FALSE
)

# Frozen sensitivity analyses.
run_sensitivity <- function(label, d, b12_cut = 300, mma_cut = 0.40, outcome = "pn_any", add_folate = FALSE) {
  d <- d[!is.na(d$b12_pg_ml) & !is.na(d$mma_umol_l), , drop = FALSE]
  d$group <- ifelse(
    d$b12_pg_ml < b12_cut,
    "low_b12",
    ifelse(d$mma_umol_l > mma_cut, "discordant", "adequate")
  )
  d$group <- factor(d$group, levels = c("adequate", "discordant", "low_b12"))
  fit <- fit_scope(d, paste0("sensitivity_", label), outcome = outcome, add_folate = add_folate)
  if (is.null(fit)) return(NULL)
  z <- fit$adjusted_pr
  m <- fit$margins
  sensitivity_rows[[length(sensitivity_rows) + 1]] <<- data.frame(
    label = label,
    outcome = outcome,
    n = nrow(fit$data),
    key_group_n = sum(fit$data$group == "discordant"),
    key_group_events = sum(fit$data[[outcome]][fit$data$group == "discordant"] == 1),
    adjusted_pr = z[["estimate"]],
    lower = z[["lower"]],
    upper = z[["upper"]],
    p = z[["p"]],
    adjusted_prevalence_difference = m[["difference"]],
    pd_lower = m[["lower"]],
    pd_upper = m[["upper"]],
    direction_concordant = is.finite(z[["estimate"]]) && z[["estimate"]] > 1 && is.finite(m[["difference"]]) && m[["difference"]] > 0,
    stringsAsFactors = FALSE
  )
}

run_sensitivity("primary_reestimated", dat)
run_sensitivity("mma_cut_0_35", dat, mma_cut = 0.35)
run_sensitivity("mma_cut_0_45", dat, mma_cut = 0.45)
run_sensitivity("b12_cut_250", dat, b12_cut = 250)
run_sensitivity("egfr_90_plus", dat[!is.na(dat$egfr) & dat$egfr >= 90, , drop = FALSE])
run_sensitivity("exclude_diabetes", dat[!is.na(dat$diabetes) & dat$diabetes == 0, , drop = FALSE])
run_sensitivity("age_60_plus", dat[!is.na(dat$age) & dat$age >= 60, , drop = FALSE])
run_sensitivity("pn_two_plus_sites", dat, outcome = "pn_two_plus")
if (mean(!is.na(dat$folate_nmol_l)) >= 0.75) run_sensitivity("additional_folate_adjustment", dat, add_folate = TRUE)
for (cy in levels(dat$cycle)) {
  run_sensitivity(paste0("leave_out_cycle_", cy), dat[!is.na(dat$cycle) & dat$cycle != cy, , drop = FALSE])
}

model_df <- do.call(rbind, model_rows)
margin_df <- do.call(rbind, margin_rows)
nonlinear_df <- do.call(rbind, nonlinear_rows)
missing_df <- do.call(rbind, missing_rows)
sensitivity_df <- if (length(sensitivity_rows)) do.call(rbind, sensitivity_rows) else data.frame()
write.csv(model_df, file.path(out_dir, "c_n06_models.csv"), row.names = FALSE)
write.csv(margin_df, file.path(out_dir, "c_n06_adjusted_margins.csv"), row.names = FALSE)
write.csv(nonlinear_df, file.path(out_dir, "c_n06_nonlinearity.csv"), row.names = FALSE)
write.csv(missing_df, file.path(out_dir, "c_n06_missingness_mechanism.csv"), row.names = FALSE)
write.csv(sensitivity_df, file.path(out_dir, "c_n06_sensitivities.csv"), row.names = FALSE)

prep <- fromJSON(file.path(out_dir, "c_n06_prep_status.json"))
pr <- primary$adjusted_pr
margin <- primary$margins
dpr <- if (is.null(discovery)) rep(NA_real_, 4) else discovery$adjusted_pr
vpr <- if (is.null(validation)) rep(NA_real_, 4) else validation$adjusted_pr

sample_gate <- prep$actual_n >= 3000
retention_gate <- prep$complete_case_retention >= 0.80
key_group_gate <- prep$key_group_n >= 250 && prep$key_group_events >= 50
validation_precision_gate <- prep$validation_key_group_n >= 80 && prep$validation_key_group_events >= 15
clinical_gate <- (
  is.finite(pr[["estimate"]]) && pr[["estimate"]] >= 1.35 && pr[["lower"]] >= 1.10
) || (
  is.finite(margin[["difference"]]) && margin[["difference"]] >= 0.05 && margin[["lower"]] >= 0.01
)
temporal_gate <- all(is.finite(c(dpr[["estimate"]], vpr[["estimate"]]))) && dpr[["estimate"]] > 1.10 && vpr[["estimate"]] > 1.10
sensitivity_gate <- nrow(sensitivity_df) >= 6 && sum(sensitivity_df$direction_concordant, na.rm = TRUE) >= ceiling(0.75 * nrow(sensitivity_df))
nonlinearity_completed <- nrow(nonlinear_df) == 1 && is.finite(nonlinear_df$p_nonlinear[1])

status <- list(
  candidate_code = "C-N06",
  data_analysis_completed = TRUE,
  prep = prep,
  primary = list(
    adjusted_prevalence_ratio = list(estimate = pr[["estimate"]], lower = pr[["lower"]], upper = pr[["upper"]], p = pr[["p"]]),
    adjusted_prevalence_difference = list(estimate = margin[["difference"]], lower = margin[["lower"]], upper = margin[["upper"]]),
    predicted_prevalence_exposed = margin[["exposed_prevalence"]],
    predicted_prevalence_reference = margin[["reference_prevalence"]]
  ),
  discovery_adjusted_pr = dpr[["estimate"]],
  validation_adjusted_pr = vpr[["estimate"]],
  sensitivity_concordant = if (nrow(sensitivity_df)) sum(sensitivity_df$direction_concordant, na.rm = TRUE) else 0,
  sensitivity_total = nrow(sensitivity_df),
  nonlinearity = nonlinear_df,
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
  data_gate_pass = all(c(sample_gate, retention_gate, key_group_gate, validation_precision_gate, clinical_gate, temporal_gate, sensitivity_gate, nonlinearity_completed)),
  final_state = if (all(c(sample_gate, retention_gate, key_group_gate, validation_precision_gate, clinical_gate, temporal_gate, sensitivity_gate, nonlinearity_completed))) "DATA_PASS_PENDING_FULL_TEXT_DEDUP" else "NO_GO_DATA_GATE"
)
write_json(status, file.path(out_dir, "c_n06_status.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
writeLines(
  c(
    "# Blind candidate C-N06 actual-data gate report",
    "",
    paste0("- Actual analytic N: ", prep$actual_n),
    paste0("- Complete-case N: ", prep$complete_case_n, " (", round(100 * prep$complete_case_retention, 1), "%)"),
    paste0("- Primary outcome events: ", prep$primary_events),
    paste0("- Key group N/events: ", prep$key_group_n, "/", prep$key_group_events),
    paste0("- Validation key group N/events: ", prep$validation_key_group_n, "/", prep$validation_key_group_events),
    paste0("- Adjusted PR: ", round(pr[["estimate"]], 3), " (95% CI ", round(pr[["lower"]], 3), " to ", round(pr[["upper"]], 3), ")"),
    paste0("- Adjusted prevalence difference: ", round(100 * margin[["difference"]], 2), " percentage points (95% CI ", round(100 * margin[["lower"]], 2), " to ", round(100 * margin[["upper"]], 2), ")"),
    paste0("- Discovery/validation adjusted PR: ", round(dpr[["estimate"]], 3), "/", round(vpr[["estimate"]], 3)),
    paste0("- Concordant sensitivities: ", status$sensitivity_concordant, "/", status$sensitivity_total),
    paste0("- Final data state: ", status$final_state),
    "",
    "The scientific title remains blinded unless full-text deduplication also passes."
  ),
  file.path(out_dir, "c_n06_report.md")
)
cat(toJSON(status, pretty = TRUE, auto_unbox = TRUE), "\n")
