from pathlib import Path

path = Path("nutrition_stage/analyze_c_n12.R")
text = path.read_text(encoding="utf-8")
old = '''  rhs <- c("group2", "age", "sex", "race", "bmi", "smoking", "hba1c", "egfr", "folate_nmol_l", "pir", "albumin_g_dl")
  if (length(unique(d$cycle)) > 1) rhs <- c(rhs, "cycle")
  if (!is.null(extra_adjust)) rhs <- c(rhs, extra_adjust)
  des <- make_design(d)
'''
new = '''  d <- droplevels(d)
  rhs <- c("group2", "age", "sex", "race", "bmi", "smoking", "hba1c", "egfr", "folate_nmol_l", "pir", "albumin_g_dl")
  if (length(unique(d$cycle[!is.na(d$cycle)])) > 1) rhs <- c(rhs, "cycle")
  if (!is.null(extra_adjust)) rhs <- c(rhs, extra_adjust)
  for (v in intersect(c("sex", "race", "smoking", "cycle"), rhs)) {
    if (length(unique(d[[v]][!is.na(d[[v]])])) < 2) rhs <- setdiff(rhs, v)
  }
  des <- make_design(d)
'''
if old not in text:
    raise SystemExit("Frozen patch target not found; refusing broad rewrite")
path.write_text(text.replace(old, new, 1), encoding="utf-8")
print("Applied invariant-factor technical patch; no cohort, threshold, outcome, or gate changed")
