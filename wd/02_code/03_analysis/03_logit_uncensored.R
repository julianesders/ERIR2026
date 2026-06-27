# ─────────────────────────────────────────────────────────────────────────────
# 03_logit_uncensored.R  — Absorbing-treatment logit robustness
#
# Plain logit with an ABSORBING treatment indicator (= 1 from first_treat
# onward) on the full AGS8 panel (no risk-set censoring). Completely
# independent of 02_hazard.R — different outcome, different sample logic.
#
# Treatment at AGS8 level, FE at AGS2 level: no collinearity.
# Absorbing dummy: treat_absorb = 1{year >= first_treat}.
#
# Specification grid — state Green vote share only.
# Direct (3 specs):
#   (1) sk + bev + chg + state_gruene + log_dens
#   (2) (1) + pers
#   (3) (1) + county_funded
# Broad (2 specs):
#   (1) sk + bev + chg + state_gruene + log_dens
#   (2) (1) + pers
#   [county-funded excluded: near-collinear with broad treat_absorb]
#
# Combined output tables (03_output/03_logit_uncensored/):
#   logit_coef.{tex,csv}   Raw coefficients — Direct (1)-(3) | Broad (1)-(2)
#   logit_ame.{tex,csv}    Average marginal effects — same layout
#   logit_diag_{direct,broad}.csv   Treated obs per state/year
# ─────────────────────────────────────────────────────────────────────────────

options("marginaleffects_safe" = FALSE)

library(data.table)
library(fixest)
library(marginaleffects)

# -- Paths --------------------------------------------------------------------

argv      <- commandArgs(trailingOnly = FALSE)
self_flag <- grep("--file=", argv, value = TRUE)
self <- if (length(self_flag)) {
  normalizePath(sub("--file=", "", self_flag))
} else if (
  requireNamespace("rstudioapi", quietly = TRUE) &&
  rstudioapi::isAvailable()
) {
  normalizePath(rstudioapi::getSourceEditorContext()$path)
} else {
  stop("Cannot determine script path. Run as: Rscript 03_logit_uncensored.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "03_output", "03_logit_uncensored")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))

# -- I/O ----------------------------------------------------------------------

.read_frame <- function(p) fread(p,
  colClasses = list(character = c("AGS8", "AGS5", "AGS2", "treat_type")))

# -- Regressor construction ---------------------------------------------------

BASE_COLS <- c("log_pop_dens", "log_steuerkraft_L1",
               "bev_stock_p100k_L1", "ev_chargepoints_p100k_L1",
               "state_gruene_L1")

zfit <- function(d) {
  bev_cap <- quantile(d$bev_stock_p100k_L1,       0.99, na.rm = TRUE)
  chg_cap <- quantile(d$ev_chargepoints_p100k_L1, 0.99, na.rm = TRUE)
  d[, bev_win := pmin(pmax(bev_stock_p100k_L1,       0), bev_cap)]
  d[, chg_win := pmin(pmax(ev_chargepoints_p100k_L1, 0), chg_cap)]
  d[, log_dens_z     := z(log_pop_dens)]
  d[, sk_z           := z(log_steuerkraft_L1)]
  d[, bev_z          := z(sqrt(bev_win))]
  d[, chg_z          := z(sqrt(chg_win))]
  d[, state_gruene_z := z(state_gruene_L1)]
  d[, pers_z         := z(log(n_vze_personal_L1))]
  d
}

# -- Labels -------------------------------------------------------------------

ET_DICT <- c(dict,
  sk_z           = "Log tax capacity (z, L1)",
  bev_z          = "$\\sqrt{\\cdot}$ BEV stock p100k (z, L1)",
  chg_z          = "$\\sqrt{\\cdot}$ EV chargers p100k (z, L1)",
  pers_z         = "Log Personnel FTE p100k (z, L1)",
  state_gruene_z = "State Green share (z, L1)",
  kreis_funded   = "County-funded (strict past)",
  log_dens_z     = "Log pop. density (z)"
)

# -- Model fitting ------------------------------------------------------------

.fit_models <- function(ph, first_treat_col, label) {
  cat(sprintf("\n=== Absorbing logit | treat: %s ===\n", first_treat_col))

  ph[, treat_absorb := as.integer(
    !is.na(get(first_treat_col)) & year >= get(first_treat_col)
  )]
  ph <- ph[complete.cases(ph[, ..BASE_COLS])]
  ph <- zfit(ph)

  n_treated <- uniqueN(ph[treat_absorb == 1L, AGS8])
  cat(sprintf(
    "  N = %d obs | %d AGS8 | %d treated AGS8 | absorb rate = %.2f%%\n",
    nrow(ph), uniqueN(ph$AGS8), n_treated,
    100 * mean(ph$treat_absorb)
  ))

  diag_dt <- rbindlist(list(
    ph[, .(n_obs = .N, n_treated_obs = sum(treat_absorb)),
       by = .(grp = AGS2)][, dim := "AGS2"],
    ph[, .(n_obs = .N, n_treated_obs = sum(treat_absorb)),
       by = .(grp = as.character(year))][, dim := "year"]
  ))
  fwrite(diag_dt, file.path(out_dir, sprintf("logit_diag_%s.csv", label)))

  fe <- "| year + AGS2"
  if (label == "direct") {
    fs <- list(
      `D(1)` = as.formula(sprintf(
        "treat_absorb ~ sk_z + bev_z + chg_z + state_gruene_z + log_dens_z %s",
        fe)),
      `D(2)` = as.formula(sprintf(
        "treat_absorb ~ sk_z + bev_z + chg_z + state_gruene_z + pers_z + log_dens_z %s",
        fe)),
      `D(3)` = as.formula(sprintf(
        "treat_absorb ~ sk_z + bev_z + chg_z + state_gruene_z + kreis_funded + log_dens_z %s",
        fe))
    )
  } else {
    fs <- list(
      `B(1)` = as.formula(sprintf(
        "treat_absorb ~ sk_z + bev_z + chg_z + state_gruene_z + log_dens_z %s",
        fe)),
      `B(2)` = as.formula(sprintf(
        "treat_absorb ~ sk_z + bev_z + chg_z + state_gruene_z + pers_z + log_dens_z %s",
        fe))
    )
  }

  models <- lapply(fs, function(f)
    feglm(f, data = ph, family = binomial("logit"), cluster = ~AGS5))
  names(models) <- names(fs)

  cat("  Computing AMEs...\n")
  ame_list <- lapply(models, function(m)
    tryCatch(as.data.frame(avg_slopes(m, vcov = ~AGS5)),
             error = function(e) { cat("  AME failed\n"); NULL }))
  names(ame_list) <- names(models)

  ph2 <- ph[!is.na(pers_z)]
  if (label == "direct") {
    ev <- list(`D(1)` = sum(ph$treat_absorb),
               `D(2)` = sum(ph2$treat_absorb),
               `D(3)` = sum(ph$treat_absorb))
  } else {
    ev <- list(`B(1)` = sum(ph$treat_absorb),
               `B(2)` = sum(ph2$treat_absorb))
  }

  list(models = models, ame_list = ame_list, ev = ev, ph = ph)
}

# -- Combined table builders --------------------------------------------------

.stars <- function(p) ifelse(is.na(p), "",
  ifelse(p < 0.01, "$^{***}$", ifelse(p < 0.05, "$^{**}$", ifelse(p < 0.1, "$^{*}$", ""))))

.coef_cell <- function(models, nm, term) {
  ct <- coeftable(models[[nm]])
  if (!term %in% rownames(ct)) return(c(est = "", se = ""))
  b <- ct[term, "Estimate"]; s <- ct[term, "Std. Error"]; p <- ct[term, ncol(ct)]
  c(est = sprintf("%.3f%s", b, .stars(p)),
    se  = sprintf("(%.3f)", s))
}

.ame_cell <- function(ame_list, nm, term) {
  am <- ame_list[[nm]]
  if (is.null(am)) return(c(est = "", se = ""))
  row <- am[am$term == term, , drop = FALSE]
  if (nrow(row) == 0L) return(c(est = "", se = ""))
  c(est = sprintf("%.4f%s", row$estimate[1L], .stars(row$p.value[1L])),
    se  = sprintf("(%.4f)", row$std.error[1L]))
}

build_combined_tables <- function(res_d, res_b) {
  models_all  <- c(res_d$models,  res_b$models)
  ame_all     <- c(res_d$ame_list, res_b$ame_list)
  ev_all      <- c(res_d$ev,      res_b$ev)
  col_keys  <- names(models_all)          # internal lookup keys
  col_names <- col_keys                   # keep for backwards-compat references
  n_d   <- length(res_d$models)
  n_b   <- length(res_b$models)
  n_all <- n_d + n_b
  # Sequential display labels (1)-(n_all), each centred with \SetCell[c=1]{c}
  col_labels <- sprintf("\\SetCell[c=1]{c} (%d)", seq_len(n_all))

  fmt_n <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)

  stats_all <- setNames(lapply(col_names, function(nm)
    list(n   = nobs(models_all[[nm]]),
         pr2 = tryCatch(as.numeric(r2(models_all[[nm]], "pr2")),
                        error = function(e) NA_real_))),
    col_names)

  write_estimates_csv(models_all,
    file.path(out_dir, "logit_coef.csv"))

  # colspec: label + n_all data columns
  colspec_str <- sprintf("l *{%d}{l}", n_all)

  # Spanning header: "Direct treatment" over n_d cols, "Broad treatment" over n_b cols
  span_d <- sprintf("\\SetCell[c=%d]{c} Direct treatment", n_d)
  span_b <- sprintf("\\SetCell[c=%d]{c} Broad treatment",  n_b)
  empty_d <- paste(rep("", n_d - 1L), collapse = " & ")
  empty_b <- paste(rep("", n_b - 1L), collapse = " & ")
  span_row <- paste0(
    " & ", span_d,
    if (n_d > 1L) paste0(" & ", empty_d) else "",
    " & ", span_b,
    if (n_b > 1L) paste0(" & ", empty_b) else "",
    " \\\\"
  )
  col_label_row <- paste0(
    " & ", paste(col_labels, collapse = " & "), " \\\\"
  )
  header_rows <- c(span_row, "\\hline", col_label_row)

  terms_ordered <- c("sk_z", "bev_z", "chg_z", "state_gruene_z",
                     "pers_z", "kreis_funded", "log_dens_z")

  .body_row_coef <- function(term) {
    lab <- ET_DICT[term]; if (is.na(lab)) lab <- term
    r1 <- c(paste0("\\quad ", lab),
             vapply(col_names, function(nm)
               .coef_cell(models_all, nm, term)["est"], character(1)))
    r2 <- c("",
             vapply(col_names, function(nm)
               .coef_cell(models_all, nm, term)["se"],  character(1)))
    c(paste(r1, collapse = " & "), " \\\\",
      paste(r2, collapse = " & "), " \\\\")
  }

  .body_row_ame <- function(term) {
    lab <- ET_DICT[term]; if (is.na(lab)) lab <- term
    r1 <- c(paste0("\\quad ", lab),
             vapply(col_names, function(nm)
               .ame_cell(ame_all, nm, term)["est"], character(1)))
    r2 <- c("",
             vapply(col_names, function(nm)
               .ame_cell(ame_all, nm, term)["se"],  character(1)))
    c(paste(r1, collapse = " & "), " \\\\",
      paste(r2, collapse = " & "), " \\\\")
  }

  body_coef <- c(
    .body_row_coef("sk_z"),   "",
    .body_row_coef("bev_z"),
    .body_row_coef("chg_z"),  "",
    .body_row_coef("state_gruene_z"), "",
    .body_row_coef("pers_z"), "",
    .body_row_coef("kreis_funded"),
    "\\hline",
    .body_row_coef("log_dens_z")
  )

  body_ame <- c(
    .body_row_ame("sk_z"),   "",
    .body_row_ame("bev_z"),
    .body_row_ame("chg_z"),  "",
    .body_row_ame("state_gruene_z"), "",
    .body_row_ame("pers_z"), "",
    .body_row_ame("kreis_funded"),
    "\\hline",
    .body_row_ame("log_dens_z")
  )

  footer_rows <- c(
    paste(c("Year/State FE", rep("Yes", n_all)), collapse = " & "),
    " \\\\",
    paste(c("Observations",
            vapply(col_names, function(nm)
              fmt_n(stats_all[[nm]]$n), character(1))),
          collapse = " & "),
    " \\\\",
    paste(c("Treated AGS8-years",
            vapply(col_names, function(nm)
              fmt_n(ev_all[[nm]]), character(1))),
          collapse = " & "),
    " \\\\",
    paste(c("Pseudo-$R^2$",
            vapply(col_names, function(nm)
              sprintf("%.3f", stats_all[[nm]]$pr2), character(1))),
          collapse = " & "),
    " \\\\"
  )

  note_coef <- paste0(
    "Raw coefficients of uncensored logit model of funding onset, periods 2016-2023. ",
    "Treatment type is specified at the top of each column. County-funded indicator excluded in ",
    "columns (4)-(5) due to collinearity with broad treatment. AGS5 and year FE included. ",
    "AGS5-clustered SEs in parentheses. $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$."
  )
  note_ame <- paste0(
    "Uncensored logit model of funding onset, periods 2016-2023. ",
    "Treatment type is specified at the top of each column. County-funded indicator excluded in ",
    "columns (4)-(5) due to collinearity with broad treatment. Average marginal effects (AME) ",
    "at sample averages reported with AGS5-clustered SEs in parentheses. ",
    "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$."
  )

  write_longtblr(
    stem        = file.path(out_dir, "logit_coef"),
    caption     = "Absorbing-Treatment Logit: Raw Coefficients (Appendix)",
    label       = "tab:logit_coef",
    note        = note_coef,
    colspec     = colspec_str,
    header_rows = header_rows,
    body_rows   = body_coef,
    footer_rows = footer_rows,
    env         = "tblr"
  )
  cat("  Coef table: logit_coef.tex\n")

  write_longtblr(
    stem        = file.path(out_dir, "logit_ame"),
    caption     = "Absorbing-Treatment Logit: Average Marginal Effects (Appendix)",
    label       = "tab:logit_ame",
    note        = note_ame,
    colspec     = colspec_str,
    header_rows = header_rows,
    body_rows   = body_ame,
    footer_rows = footer_rows,
    env         = "tblr"
  )
  cat("  AME table: logit_ame.tex\n")
}

# -- Run ----------------------------------------------------------------------

ph <- .read_frame(file.path(data_final, "frame_logit_full.csv"))

res_direct <- .fit_models(copy(ph), "first_treat_direct", "direct")
res_broad  <- .fit_models(copy(ph), "first_treat_broad",  "broad")

build_combined_tables(res_direct, res_broad)

cat(sprintf("\nAbsorbing logit outputs -> %s\n", out_dir))