options(stringsAsFactors = FALSE)
options(survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(jsonlite)
})

out_dir <- "nutrition_screen/output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dat <- read.csv(file.path(out_dir, "c_n05_analysis_core.csv"), na.strings = c("", "NA"))

num_vars <- c("age", "race", "weight", "psu_u", "strata_u", "stfr", "ferritin", "hgb", "mcv", "crp_mg_l", "log_crp", "bmi", "pir", "diabetes", "anemia", "low_mcv", "log_stfr")
for (v in intersect(num_vars, names(dat))) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
dat$complete_case <- as.logical(dat$complete_case)
dat$cycle <- factor(dat$cycle, levels = c("C", "D", "E", "F"))
dat$period <- factor(dat$period, levels = c("discovery", "validation"))
dat$race <- factor(dat$race)
dat$smoking <- factor(dat$smoking, levels = c("never", "ever"))
dat$group <- factor(dat$group, levels = c("adequate", "discordant", "ferritin_low"))
dat$log_ferritin <- log(dat$ferritin)

coef_name <- "groupdiscordant"
model_rows <- list()
sens_rows <- list()
nonlinear_rows <- list()
missing_rows <- list()

safe_ci <- function(model, term) {
  cf <- coef(model)
  if (!(term %in% names(cf))) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  se <- sqrt(diag(vcov(model)))[term]
  est <- unname(cf[term])
  c(est, est - 1.96 * se, est + 1.96 * se, unname(summary(model)$coefficients[term, ncol(summary(model)$coefficients)]))
}

append_effect <- function(scope, model_name, model, term = coef_name, exponentiate = FALSE) {
  z <- safe_ci(model, term)
  if (exponentiate) z[1:3] <- exp(z[1:3])
  model_rows[[length(model_rows) + 1]] <<- data.frame(
    candidate_code = "C-N05", scope = scope, model = model_name, term = term,
    estimate = z[1], lower = z[2], upper = z[3], p = z[4], n = nrow(model$model),
    stringsAsFactors = FALSE
  )
  invisible(z)
}

make_design <- function(d) {
  svydesign(ids = ~psu_u, strata = ~strata_u, weights = ~weight, nest = TRUE, data = d)
}

full_formula_hgb <- hgb ~ group + age + race + bmi + log_crp + smoking + pir + cycle
full_formula_anemia <- anemia ~ group + age + race + bmi + log_crp + smoking + pir + cycle

fit_scope <- function(d, scope) {
  d <- d[d$complete_case & !is.na(d$group), , drop = FALSE]
  if (nrow(d) < 200 || length(unique(d$psu_u)) < 4 || sum(d$group == "discordant") < 20) return(NULL)
  d$group <- droplevels(factor(d$group, levels = c("adequate", "discordant", "ferritin_low")))
  des <- make_design(d)
  crude_h <- svyglm(hgb ~ group + cycle, design = des, family = gaussian())
  crude_a <- svyglm(anemia ~ group + cycle, design = des, family = quasibinomial())
  adj_h <- svyglm(full_formula_hgb, design = des, family = gaussian())
  adj_a <- svyglm(full_formula_anemia, design = des, family = quasibinomial())
  append_effect(scope, "crude_hgb_difference_g_dl", crude_h)
  append_effect(scope, "crude_anemia_or", crude_a, exponentiate = TRUE)
  zh <- append_effect(scope, "adjusted_hgb_difference_g_dl", adj_h)
  za <- append_effect(scope, "adjusted_anemia_or", adj_a, exponentiate = TRUE)
  list(data = d, design = des, adjusted_hgb = zh, adjusted_anemia = za, hgb_model = adj_h, anemia_model = adj_a)
}

primary <- fit_scope(dat, "overall")
discovery <- fit_scope(dat[dat$period == "discovery", ], "discovery")
validation <- fit_scope(dat[dat$period == "validation", ], "validation")
if (is.null(primary)) stop("Primary model could not be estimated")

# Missingness mechanism: whether full-adjustment variables are jointly observed.
des_all <- make_design(dat[!is.na(dat$group) & !is.na(dat$age) & !is.na(dat$race), ])
miss_fit <- svyglm(I(complete_case) ~ group + age + race + cycle, design = des_all, family = quasibinomial())
for (term in c("groupdiscordant", "groupferritin_low")) {
  z <- safe_ci(miss_fit, term)
  missing_rows[[length(missing_rows) + 1]] <- data.frame(
    term = term, odds_ratio = exp(z[1]), lower = exp(z[2]), upper = exp(z[3]), p = z[4], stringsAsFactors = FALSE
  )
}

# Restricted cubic spline nonlinearity among ferritin-normal women.
weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  ord <- order(x); x <- x[ord]; w <- w[ord]
  cw <- cumsum(w) / sum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}

rcs_terms <- function(x, knots) {
  k <- knots
  pos3 <- function(z) pmax(z, 0)^3
  den <- (k[length(k)] - k[1])^2
  h <- sapply(seq_len(length(k) - 2), function(j) {
    (pos3(x - k[j]) - pos3(x - k[length(k) - 1]) * (k[length(k)] - k[j]) / (k[length(k)] - k[length(k) - 1]) +
       pos3(x - k[length(k)]) * (k[length(k) - 1] - k[j]) / (k[length(k)] - k[length(k) - 1])) / den
  })
  colnames(h) <- paste0("rcs", seq_len(ncol(h)))
  h
}

nl <- primary$data[primary$data$ferritin >= 15, , drop = FALSE]
knots <- weighted_quantile(nl$log_stfr, nl$weight, c(0.05, 0.35, 0.65, 0.95))
basis <- rcs_terms(nl$log_stfr, knots)
nl$xlin <- nl$log_stfr
nl$rcs1 <- basis[, 1]
nl$rcs2 <- basis[, 2]
des_nl <- make_design(nl)
nl_h <- svyglm(hgb ~ xlin + rcs1 + rcs2 + age + race + bmi + log_crp + smoking + pir + cycle, design = des_nl, family = gaussian())
nl_a <- svyglm(anemia ~ xlin + rcs1 + rcs2 + age + race + bmi + log_crp + smoking + pir + cycle, design = des_nl, family = quasibinomial())
test_h <- regTermTest(nl_h, ~rcs1 + rcs2)
test_a <- regTermTest(nl_a, ~rcs1 + rcs2)
nonlinear_rows[[1]] <- data.frame(outcome = "hemoglobin", n = nrow(nl), p_nonlinear = unname(test_h$p), knots = paste(round(exp(knots), 3), collapse = "|"), stringsAsFactors = FALSE)
nonlinear_rows[[2]] <- data.frame(outcome = "anemia", n = nrow(nl), p_nonlinear = unname(test_a$p), knots = paste(round(exp(knots), 3), collapse = "|"), stringsAsFactors = FALSE)

# Sensitivity helper: thresholds are frozen before each model is fitted.
run_sensitivity <- function(label, d, ferritin_cut = 15, stfr_cut = 5.33, add_log_ferritin = FALSE) {
  d <- d[d$complete_case, , drop = FALSE]
  d$group2 <- ifelse(d$ferritin < ferritin_cut, "ferritin_low", ifelse(d$stfr > stfr_cut, "discordant", "adequate"))
  d$group2 <- factor(d$group2, levels = c("adequate", "discordant", "ferritin_low"))
  if (nrow(d) < 200 || sum(d$group2 == "discordant") < 20) return(NULL)
  des <- make_design(d)
  if (add_log_ferritin) {
    mh <- svyglm(hgb ~ group2 + log_ferritin + age + race + bmi + log_crp + smoking + pir + cycle, design = des, family = gaussian())
    ma <- svyglm(anemia ~ group2 + log_ferritin + age + race + bmi + log_crp + smoking + pir + cycle, design = des, family = quasibinomial())
  } else {
    mh <- svyglm(hgb ~ group2 + age + race + bmi + log_crp + smoking + pir + cycle, design = des, family = gaussian())
    ma <- svyglm(anemia ~ group2 + age + race + bmi + log_crp + smoking + pir + cycle, design = des, family = quasibinomial())
  }
  zh <- safe_ci(mh, "group2discordant")
  za <- safe_ci(ma, "group2discordant"); za[1:3] <- exp(za[1:3])
  sens_rows[[length(sens_rows) + 1]] <<- data.frame(
    label = label, n = nrow(d), key_group_n = sum(d$group2 == "discordant"), key_group_anemia = sum(d$anemia[d$group2 == "discordant"]),
    hgb_difference = zh[1], hgb_lower = zh[2], hgb_upper = zh[3], hgb_p = zh[4],
    anemia_or = za[1], anemia_lower = za[2], anemia_upper = za[3], anemia_p = za[4],
    direction_concordant = is.finite(zh[1]) && is.finite(za[1]) && zh[1] < 0 && za[1] > 1,
    stringsAsFactors = FALSE
  )
}

run_sensitivity("primary_reestimated", dat)
run_sensitivity("ferritin_cut_30", dat, ferritin_cut = 30)
run_sensitivity("exclude_crp_gt_10_mg_l", dat[is.na(dat$crp_mg_l) | dat$crp_mg_l <= 10, ])
run_sensitivity("exclude_diabetes", dat[dat$diabetes == 0, ])
run_sensitivity("cycles_2005_2010", dat[dat$cycle != "C", ])
run_sensitivity("stfr_cut_5_0", dat, stfr_cut = 5.0)
run_sensitivity("stfr_cut_5_5", dat, stfr_cut = 5.5)
run_sensitivity("adjust_for_log_ferritin", dat, add_log_ferritin = TRUE)
for (cy in levels(dat$cycle)) run_sensitivity(paste0("leave_out_cycle_", cy), dat[dat$cycle != cy, ])

model_df <- do.call(rbind, model_rows)
sens_df <- do.call(rbind, sens_rows)
nonlinear_df <- do.call(rbind, nonlinear_rows)
missing_df <- do.call(rbind, missing_rows)
write.csv(model_df, file.path(out_dir, "c_n05_models.csv"), row.names = FALSE)
write.csv(sens_df, file.path(out_dir, "c_n05_sensitivities.csv"), row.names = FALSE)
write.csv(nonlinear_df, file.path(out_dir, "c_n05_nonlinearity.csv"), row.names = FALSE)
write.csv(missing_df, file.path(out_dir, "c_n05_missingness_mechanism.csv"), row.names = FALSE)

prep <- fromJSON(file.path(out_dir, "c_n05_prep_status.json"))
ph <- primary$adjusted_hgb
pa <- primary$adjusted_anemia
dh <- if (is.null(discovery)) rep(NA_real_, 4) else discovery$adjusted_hgb
da <- if (is.null(discovery)) rep(NA_real_, 4) else discovery$adjusted_anemia
vh <- if (is.null(validation)) rep(NA_real_, 4) else validation$adjusted_hgb
va <- if (is.null(validation)) rep(NA_real_, 4) else validation$adjusted_anemia

sample_gate <- prep$core_n >= 3000
retention_gate <- prep$complete_case_retention >= 0.80
key_group_gate <- prep$key_group_n >= 250 && prep$key_group_anemia_events >= 30
clinical_gate <- (is.finite(ph[1]) && ph[1] <= -0.50 && ph[3] <= -0.20) || (is.finite(pa[1]) && pa[1] >= 1.50 && pa[2] >= 1.20)
temporal_gate <- all(is.finite(c(dh[1], da[1], vh[1], va[1]))) && dh[1] < 0 && da[1] > 1 && vh[1] < 0 && va[1] > 1
sens_gate <- sum(sens_df$direction_concordant, na.rm = TRUE) >= ceiling(0.70 * nrow(sens_df))
nonlinear_completed <- nrow(nonlinear_df) == 2 && all(is.finite(nonlinear_df$p_nonlinear))

status <- list(
  candidate_code = "C-N05",
  data_analysis_completed = TRUE,
  prep = prep,
  primary = list(
    adjusted_hgb_difference = list(estimate = ph[1], lower = ph[2], upper = ph[3], p = ph[4]),
    adjusted_anemia_or = list(estimate = pa[1], lower = pa[2], upper = pa[3], p = pa[4])
  ),
  discovery = list(hgb_difference = dh[1], anemia_or = da[1]),
  validation = list(hgb_difference = vh[1], anemia_or = va[1]),
  nonlinearity = nonlinear_df,
  sensitivity_concordant = sum(sens_df$direction_concordant, na.rm = TRUE),
  sensitivity_total = nrow(sens_df),
  gates = list(
    sample_gate = sample_gate,
    complete_case_retention_gate = retention_gate,
    key_group_event_gate = key_group_gate,
    clinical_magnitude_and_precision_gate = clinical_gate,
    temporal_replication_gate = temporal_gate,
    sensitivity_direction_gate = sens_gate,
    nonlinearity_completed = nonlinear_completed
  ),
  data_gate_pass = all(c(sample_gate, retention_gate, key_group_gate, clinical_gate, temporal_gate, sens_gate, nonlinear_completed)),
  final_state = if (all(c(sample_gate, retention_gate, key_group_gate, clinical_gate, temporal_gate, sens_gate, nonlinear_completed))) "DATA_PASS_PENDING_FULL_TEXT_DEDUP" else "NO_GO_DATA_GATE"
)
write_json(status, file.path(out_dir, "c_n05_status.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")

report <- c(
  "# Blind nutrition candidate C-N05: actual data gate report",
  "",
  paste0("- Actual classified cohort: ", prep$core_n),
  paste0("- Complete-case full model: ", prep$complete_case_n, " (", round(100 * prep$complete_case_retention, 1), "%)"),
  paste0("- Key group: ", prep$key_group_n, "; anemia events: ", prep$key_group_anemia_events),
  paste0("- Adjusted hemoglobin difference: ", round(ph[1], 3), " g/dL (95% CI ", round(ph[2], 3), " to ", round(ph[3], 3), ")"),
  paste0("- Adjusted anemia OR: ", round(pa[1], 3), " (95% CI ", round(pa[2], 3), " to ", round(pa[3], 3), ")"),
  paste0("- Discovery/validation hemoglobin differences: ", round(dh[1], 3), " / ", round(vh[1], 3)),
  paste0("- Discovery/validation anemia ORs: ", round(da[1], 3), " / ", round(va[1], 3)),
  paste0("- Concordant sensitivities: ", sum(sens_df$direction_concordant, na.rm = TRUE), "/", nrow(sens_df)),
  paste0("- Data gate: ", status$final_state),
  "",
  "No title is released unless the candidate also passes full-text deduplication."
)
writeLines(report, file.path(out_dir, "c_n05_report.md"))
cat(toJSON(status, pretty = TRUE, auto_unbox = TRUE), "\n")
