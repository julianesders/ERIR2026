# ─────────────────────────────────────────────────────────────────────────────
# 02_hazard.R   — Part (i): drivers of EMK onset
#
# Discrete-time hazard on the AGS8 risk set (year >= 2015).
# All specs: year + AGS2 FE, SEs clustered on AGS5.
# Primary event: onset_direct (frame_hazard.csv).
#
# 3-column specification grid (log_dens_z universal, listed last):
#   (1) sk + bev + chg + state_gruene + log_dens    ← baseline
#   (2) (1) + pers                                  ← no-Stadtstaaten
#   (3) (1) + kreis_funded
#
# Tables produced in 03_output/02_hazard/:
#   hazard_hr_ame.{tex,csv}              HR + AME, 3 specs (main text)
#   hazard_coef.{tex,csv}                Raw cloglog coefficients (appendix)
#   hazard_logit_hr_ame.{tex,csv}        Logit link OR + AME (appendix)
#   hazard_logit_coef.{tex,csv}          Logit link raw coefficients (appendix)
#   hazard_broad_hr_ame.{tex,csv}        Broad treatment HR + AME, 2 specs (appendix)
#   hazard_broad_coef.{tex,csv}          Broad treatment raw coef (appendix)
#   tab_hazard_diag_{main,appendix}.csv  events per state/year
# ─────────────────────────────────────────────────────────────────────────────

options("marginaleffects_safe" = FALSE)

library(data.table)
library(fixest)
library(marginaleffects)
library(ggplot2)

# -- Paths ────────────────────────────────────────────────────────────────────

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
  stop("Cannot determine script path. Run as: Rscript 02_hazard.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "03_output", "02_hazard")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))

# -- I/O ----------------------------------------------------------------------

.read_frame <- function(p) fread(p,
  colClasses = list(character = c("AGS8", "AGS5", "AGS2", "treat_type")))

# -- Regressor construction ---------------------------------------------------
# BEV stock and charger covariates: winsorised at the 99th pct, then
# square-root transformed before z-scoring. Sqrt is the variance-stabilising
# transform for Poisson-distributed count data; sqrt(0) = 0 preserves the
# zero mass.

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

# -- etable display settings --------------------------------------------------

ET_DICT <- c(dict,
  sk_z           = "Log tax capacity p.c. (z, L1)",
  bev_z          = "$\\sqrt{\\cdot}$ BEV stock p100k (z, L1)",
  chg_z          = "$\\sqrt{\\cdot}$ EV chargers p100k (z, L1)",
  pers_z         = "Log Personnel p.c. (z, L1)",
  state_gruene_z = "State Green share (z, L1)",
  kreis_funded   = "County funded",
  log_dens_z     = "Log pop. density (z)"
)

# Groups for the 3-spec main model (state Grüne only)
HAZARD_GROUPS <- list(
  "Fiscal capacity"         = "sk_z",
  "EV ecosystem"            = c("bev_z", "chg_z"),
  "Political"               = "state_gruene_z",
  "Administrative capacity" = "pers_z",
  "County coverage"         = "kreis_funded",
  .hline                    = NULL,
  "Controls"                = "log_dens_z"
)

# -- Diagnostics helper -------------------------------------------------------

.diag <- function(ph, event_col, label) {
  diag_dt <- rbindlist(list(
    ph[, .(n_obs = .N, n_events = sum(get(event_col))),
       by = .(grp = AGS2)][, dim := "AGS2"],
    ph[, .(n_obs = .N, n_events = sum(get(event_col))),
       by = .(grp = as.character(year))][, dim := "year"]
  ))
  fwrite(diag_dt,
         file.path(out_dir, sprintf("tab_hazard_diag_%s.csv", label)))
  zero <- diag_dt[dim == "AGS2" & n_events == 0L]
  if (nrow(zero))
    cat(sprintf("  NOTE: %d AGS2 with 0 events\n", nrow(zero)))
}

# -- Logit coef+AME table (same 3-spec structure as main cloglog) ------------

run_logit_main <- function(ph, event_col, stem_hr_ame, stem_coef) {
  cat(sprintf("\n=== Logit hazard table | event: %s ===\n", event_col))
  ph <- ph[complete.cases(ph[, ..BASE_COLS])]
  ph <- zfit(ph)
  cat(sprintf("  N = %d obs | %d AGS8 | %d events (%.2f%%)\n",
              nrow(ph), uniqueN(ph$AGS8),
              sum(ph[[event_col]]), 100 * mean(ph[[event_col]])))
  .diag(ph, event_col, basename(stem_coef))

  fe  <- "| year + AGS2"
  fs3 <- list(
    `(1)` = as.formula(sprintf(
      "%s ~ sk_z + bev_z + chg_z + state_gruene_z + log_dens_z %s",
      event_col, fe)),
    `(2)` = as.formula(sprintf(
      "%s ~ sk_z + bev_z + chg_z + state_gruene_z + pers_z + log_dens_z %s",
      event_col, fe)),
    `(3)` = as.formula(sprintf(
      "%s ~ sk_z + bev_z + chg_z + state_gruene_z + kreis_funded + log_dens_z %s",
      event_col, fe))
  )
  models    <- lapply(fs3, function(f)
    feglm(f, data = ph, family = binomial("logit"), cluster = ~AGS5))
  names(models) <- names(fs3)
  col_names <- names(models)

  cat("  Computing AMEs...\n")
  ame_list <- lapply(models, function(m)
    tryCatch(as.data.frame(avg_slopes(m, vcov = ~AGS5)),
             error = function(e) { cat("  AME failed\n"); NULL }))
  names(ame_list) <- col_names

  ph2 <- ph[!is.na(pers_z)]
  ev  <- list(`(1)` = sum(ph[[event_col]]),
              `(2)` = sum(ph2[[event_col]]),
              `(3)` = sum(ph[[event_col]]))
  stats <- setNames(lapply(col_names, function(nm)
    list(n   = nobs(models[[nm]]),
         pr2 = tryCatch(as.numeric(r2(models[[nm]], "pr2")),
                        error = function(e) NA_real_))),
    col_names)

  write_estimates_csv(models, paste0(stem_coef, ".csv"))

  .stars <- function(p) ifelse(is.na(p), "",
    ifelse(p < 0.01, "$^{***}$", ifelse(p < 0.05, "$^{**}$", ifelse(p < 0.1, "$^{*}$", ""))))

  .or_cell <- function(nm, term) {
    ct <- coeftable(models[[nm]])
    if (!term %in% rownames(ct)) return(c(est = "", ci = ""))
    b <- ct[term, "Estimate"]; s <- ct[term, "Std. Error"]
    c(est = sprintf("%.2f", exp(b)),
      ci  = sprintf("[%.2f,\\,%.2f]", exp(b - 1.96 * s), exp(b + 1.96 * s)))
  }

  .ame_cell <- function(nm, term) {
    am <- ame_list[[nm]]
    if (is.null(am)) return(c(est = "", se = ""))
    row <- am[am$term == term, , drop = FALSE]
    if (nrow(row) == 0L) return(c(est = "", se = ""))
    c(est = sprintf("%.4f%s", row$estimate[1L], .stars(row$p.value[1L])),
      se  = sprintf("(%.4f)", row$std.error[1L]))
  }

  .crc <- function(term) {
    lab <- ET_DICT[term]; if (is.na(lab)) lab <- term
    r1  <- c(paste0("\\quad ", lab),
             unlist(lapply(col_names, function(nm)
               c(.or_cell(nm, term)["est"], .ame_cell(nm, term)["est"]))))
    r2  <- c("",
             unlist(lapply(col_names, function(nm)
               c(.or_cell(nm, term)["ci"], .ame_cell(nm, term)["se"]))))
    c(paste(r1, collapse = " & "), " \\\\",
      paste(r2, collapse = " & "), " \\\\")
  }

  fmt_n <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)
  .s2 <- function(v, a = "r") paste0("\\SetCell[c=2]{", a, "} ", v)

  header_c <- c(
    paste0(" & \\SetCell[c=2]{c} (1) & & \\SetCell[c=2]{c} (2)",
           " & & \\SetCell[c=2]{c} (3) & \\\\"),
    "\\hline",
    " & OR & AME & OR & AME & OR & AME \\\\"
  )
  body_c <- c(
    .crc("sk_z"),
    "",
    .crc("bev_z"),
    .crc("chg_z"),
    "",
    .crc("state_gruene_z"),
    "",
    .crc("pers_z"),
    "",
    .crc("kreis_funded"),
    "\\hline",
    .crc("log_dens_z")
  )
  footer_c <- c(
    paste(c("Year/State FE",
            .s2("Yes","c"),"", .s2("Yes","c"),"", .s2("Yes","c"),""),
          collapse = " & "), " \\\\",
    paste(c("Observations",
            .s2(fmt_n(stats[["(1)"]]$n),"c"),"",
            .s2(fmt_n(stats[["(2)"]]$n),"c"),"",
            .s2(fmt_n(stats[["(3)"]]$n),"c"),""),
          collapse = " & "), " \\\\",
    paste(c("Onset events",
            .s2(fmt_n(ev[["(1)"]]),"c"),"",
            .s2(fmt_n(ev[["(2)"]]),"c"),"",
            .s2(fmt_n(ev[["(3)"]]),"c"),""),
          collapse = " & "), " \\\\",
    paste(c("Pseudo-$R^2$",
            .s2(sprintf("%.3f", stats[["(1)"]]$pr2),"c"),"",
            .s2(sprintf("%.3f", stats[["(2)"]]$pr2),"c"),"",
            .s2(sprintf("%.3f", stats[["(3)"]]$pr2),"c"),""),
          collapse = " & "), " \\\\"
  )
  write_longtblr(
    stem        = stem_hr_ame,
    caption     = "Discrete-Time Hazard: Logit Link, Odds Ratios and Average Marginal Effects",
    label       = "tab:hazard_logit_hr_ame",
    note        = paste0(
      "Discrete-time hazard (logit link) of direct funding onset, AGS8 risk set 2016-2023. ",
      "AGS5 and year FE included. Odds ratios (OR) are computed as $\\exp(\\hat{\\beta})$, ",
      "95\\% CIs in brackets. Average marginal effects (AME) at sample averages reported ",
      "with AGS5-clustered SEs in parentheses. $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$."),
    colspec     = "l r l r l r l",
    header_rows = header_c,
    body_rows   = body_c,
    footer_rows = footer_c,
    resize      = TRUE
  )
  cat(sprintf("  OR/AME table: %s.tex\n", stem_hr_ame))

  ev_counts <- setNames(lapply(col_names, function(nm) ev[[nm]]), col_names)
  write_coef_longtblr(
    models       = models,
    stem         = stem_coef,
    caption      = "Discrete-Time Hazard: Logit Link, Raw Coefficients",
    label        = "tab:hazard_logit_coef",
    note         = paste0(
      "Raw coefficients of discrete-time hazard (logit link) of direct funding onset, ",
      "AGS8 risk set 2016-2023. AGS5 and year FE included. AGS5-clustered SEs in parentheses. ",
      "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$."),
    groups       = HAZARD_GROUPS,
    var_labels   = ET_DICT,
    event_counts = ev_counts
  )
  cat(sprintf("  Raw coef table: %s.tex\n", stem_coef))

  invisible(list(models = models, ph = ph, ame_list = ame_list))
}

# -- Combined HR/AME main table + raw-coef appendix table --------------------

run_hazard_main <- function(ph, event_col, stem_hr_ame, stem_coef) {
  cat(sprintf("\n=== Main hazard table | event: %s ===\n", event_col))
  ph <- ph[complete.cases(ph[, ..BASE_COLS])]
  ph <- zfit(ph)
  cat(sprintf("  N = %d obs | %d AGS8 | %d events (%.2f%%)\n",
              nrow(ph), uniqueN(ph$AGS8),
              sum(ph[[event_col]]), 100 * mean(ph[[event_col]])))
  .diag(ph, event_col, basename(stem_coef))

  fe  <- "| year + AGS2"
  fs3 <- list(
    `(1)` = as.formula(sprintf(
      "%s ~ sk_z + bev_z + chg_z + state_gruene_z + log_dens_z %s",
      event_col, fe)),
    `(2)` = as.formula(sprintf(
      "%s ~ sk_z + bev_z + chg_z + state_gruene_z + pers_z + log_dens_z %s",
      event_col, fe)),
    `(3)` = as.formula(sprintf(
      "%s ~ sk_z + bev_z + chg_z + state_gruene_z + kreis_funded + log_dens_z %s",
      event_col, fe))
  )
  models    <- lapply(fs3, function(f)
    feglm(f, data = ph, family = binomial("cloglog"), cluster = ~AGS5))
  names(models) <- names(fs3)
  col_names <- names(models)

  cat("  Computing AMEs...\n")
  ame_list <- lapply(models, function(m)
    tryCatch(as.data.frame(avg_slopes(m, vcov = ~AGS5)),
             error = function(e) { cat("  AME failed\n"); NULL }))
  names(ame_list) <- col_names

  ph2 <- ph[!is.na(pers_z)]
  ev  <- list(`(1)` = sum(ph[[event_col]]),
              `(2)` = sum(ph2[[event_col]]),
              `(3)` = sum(ph[[event_col]]))
  stats <- setNames(lapply(col_names, function(nm)
    list(n   = nobs(models[[nm]]),
         ev  = ev[[nm]],
         pr2 = tryCatch(as.numeric(r2(models[[nm]], "pr2")),
                        error = function(e) NA_real_))),
    col_names)

  write_estimates_csv(models, paste0(stem_coef, ".csv"))

  # ── A. HR/AME combined table (7-col: label + 3*(HR, AME)) ─────────────────
  N_C    <- 7L
  .stars <- function(p) ifelse(is.na(p), "",
    ifelse(p < 0.01, "$^{***}$", ifelse(p < 0.05, "$^{**}$", ifelse(p < 0.1, "$^{*}$", ""))))

  .hr_cell <- function(nm, term) {
    ct <- coeftable(models[[nm]])
    if (!term %in% rownames(ct)) return(c(est = "", ci = ""))
    b <- ct[term, "Estimate"]; s <- ct[term, "Std. Error"]
    c(est = sprintf("%.2f", exp(b)),
      ci  = sprintf("[%.2f,\\,%.2f]", exp(b - 1.96 * s), exp(b + 1.96 * s)))
  }

  .ame_cell <- function(nm, term) {
    am <- ame_list[[nm]]
    if (is.null(am)) return(c(est = "", se = ""))
    row <- am[am$term == term, , drop = FALSE]
    if (nrow(row) == 0L) return(c(est = "", se = ""))
    c(est = sprintf("%.4f%s", row$estimate[1L], .stars(row$p.value[1L])),
      se  = sprintf("(%.4f)", row$std.error[1L]))
  }

  .crc <- function(term) {
    lab <- ET_DICT[term]; if (is.na(lab)) lab <- term
    r1  <- c(paste0("\\quad ", lab),
             unlist(lapply(col_names, function(nm)
               c(.hr_cell(nm, term)["est"], .ame_cell(nm, term)["est"]))))
    r2  <- c("",
             unlist(lapply(col_names, function(nm)
               c(.hr_cell(nm, term)["ci"], .ame_cell(nm, term)["se"]))))
    c(paste(r1, collapse = " & "), " \\\\",
      paste(r2, collapse = " & "), " \\\\")
  }

  fmt_n <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)
  .s2 <- function(v, a = "r") paste0("\\SetCell[c=2]{", a, "} ", v)

  header_c <- c(
    paste0(" & \\SetCell[c=2]{c} (1) & & \\SetCell[c=2]{c} (2)",
           " & & \\SetCell[c=2]{c} (3) & \\\\"),
    "\\hline",
    " & HR & AME & HR & AME & HR & AME \\\\"
  )

  body_c <- c(
    .crc("sk_z"),
    "",
    .crc("bev_z"),
    .crc("chg_z"),
    "",
    .crc("state_gruene_z"),
    "",
    .crc("pers_z"),
    "",
    .crc("kreis_funded"),
    "\\hline",
    .crc("log_dens_z")
  )

  footer_c <- c(
    paste(c("Year/State FE",
            .s2("Yes","c"),"", .s2("Yes","c"),"", .s2("Yes","c"),""),
          collapse = " & "), " \\\\",
    paste(c("Observations",
            .s2(fmt_n(stats[["(1)"]]$n),"c"),"",
            .s2(fmt_n(stats[["(2)"]]$n),"c"),"",
            .s2(fmt_n(stats[["(3)"]]$n),"c"),""),
          collapse = " & "), " \\\\",
    paste(c("Onset events",
            .s2(fmt_n(ev[["(1)"]]),"c"),"",
            .s2(fmt_n(ev[["(2)"]]),"c"),"",
            .s2(fmt_n(ev[["(3)"]]),"c"),""),
          collapse = " & "), " \\\\",
    paste(c("Pseudo-$R^2$",
            .s2(sprintf("%.3f", stats[["(1)"]]$pr2),"c"),"",
            .s2(sprintf("%.3f", stats[["(2)"]]$pr2),"c"),"",
            .s2(sprintf("%.3f", stats[["(3)"]]$pr2),"c"),""),
          collapse = " & "), " \\\\"
  )

  note_c <- paste0(
    "Discrete-time hazard (cloglog link) of direct funding onset, AGS8 risk set 2016-2023. ",
    "AGS5 and year FE included. Hazard ratios (HR) computed as $\\exp(\\hat{\\beta})$, ",
    "95\\% CIs in brackets. Average marginal effects (AME) at sample averages reported ",
    "with AGS5-clustered SEs in parentheses. $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$."
  )

  write_longtblr(
    stem        = stem_hr_ame,
    caption     = "Discrete-Time Hazard: Hazard Ratios and Average Marginal Effects",
    label       = "tab:hazard_hr_ame",
    note        = note_c,
    colspec     = "l r l r l r l",
    header_rows = header_c,
    body_rows   = body_c,
    footer_rows = footer_c,
    resize      = TRUE
  )
  cat(sprintf("  HR/AME table: %s.tex\n", stem_hr_ame))

  # ── B. Raw cloglog coefficients (appendix) ─────────────────────────────────
  note_coef <- paste0(
    "Raw coefficients of discrete-time hazard (cloglog link) of direct funding onset, ",
    "AGS8 risk set 2016-2023. AGS5 and year FE included. AGS5-clustered SEs in parentheses. ",
    "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$."
  )
  ev_counts <- setNames(lapply(col_names, function(nm) ev[[nm]]), col_names)
  write_coef_longtblr(
    models       = models,
    stem         = stem_coef,
    caption      = "Discrete-Time Hazard: Raw Coefficients",
    label        = "tab:hazard_coef",
    note         = note_coef,
    groups       = HAZARD_GROUPS,
    var_labels   = ET_DICT,
    event_counts = ev_counts
  )
  cat(sprintf("  Raw coef table: %s.tex\n", stem_coef))

  invisible(list(models = models, ph = ph, ame_list = ame_list))
}

# -- Broad-treatment HR/AME + coef (2 specs, no county-funded) ----------------

HAZARD_GROUPS_BROAD <- list(
  "Fiscal capacity"         = "sk_z",
  "EV ecosystem"            = c("bev_z", "chg_z"),
  "Political"               = "state_gruene_z",
  "Administrative capacity" = "pers_z",
  .hline                    = NULL,
  "Controls"                = "log_dens_z"
)

run_hazard_broad <- function(ph_broad, stem_hr_ame, stem_coef) {
  cat(sprintf("\n=== Broad hazard HR/AME | event: onset_broad ===\n"))
  ph_brd <- ph_broad[complete.cases(ph_broad[, ..BASE_COLS])]
  ph_brd <- zfit(ph_brd)
  cat(sprintf("  N = %d | %d AGS8 | %d onset_broad (%.2f%%)\n",
              nrow(ph_brd), uniqueN(ph_brd$AGS8),
              sum(ph_brd$onset_broad), 100 * mean(ph_brd$onset_broad)))

  fe  <- "| year + AGS2"
  fs2 <- list(
    `(1)` = as.formula(sprintf(
      "onset_broad ~ sk_z + bev_z + chg_z + state_gruene_z + log_dens_z %s",
      fe)),
    `(2)` = as.formula(sprintf(
      "onset_broad ~ sk_z + bev_z + chg_z + state_gruene_z + pers_z + log_dens_z %s",
      fe))
  )
  models    <- lapply(fs2, function(f)
    feglm(f, data = ph_brd, family = binomial("cloglog"), cluster = ~AGS5))
  names(models) <- names(fs2)
  col_names <- names(models)

  cat("  Computing AMEs...\n")
  ame_list <- lapply(models, function(m)
    tryCatch(as.data.frame(avg_slopes(m, vcov = ~AGS5)),
             error = function(e) { cat("  AME failed\n"); NULL }))
  names(ame_list) <- col_names

  ph2_brd <- ph_brd[!is.na(pers_z)]
  ev  <- list(`(1)` = sum(ph_brd$onset_broad),
              `(2)` = sum(ph2_brd$onset_broad))
  stats <- setNames(lapply(col_names, function(nm)
    list(n   = nobs(models[[nm]]),
         pr2 = tryCatch(as.numeric(r2(models[[nm]], "pr2")),
                        error = function(e) NA_real_))),
    col_names)

  write_estimates_csv(models, paste0(stem_coef, ".csv"))

  .stars <- function(p) ifelse(is.na(p), "",
    ifelse(p < 0.01, "$^{***}$", ifelse(p < 0.05, "$^{**}$", ifelse(p < 0.1, "$^{*}$", ""))))

  .hr_cell <- function(nm, term) {
    ct <- coeftable(models[[nm]])
    if (!term %in% rownames(ct)) return(c(est = "", ci = ""))
    b <- ct[term, "Estimate"]; s <- ct[term, "Std. Error"]
    c(est = sprintf("%.2f", exp(b)),
      ci  = sprintf("[%.2f,\\,%.2f]", exp(b - 1.96 * s), exp(b + 1.96 * s)))
  }

  .ame_cell <- function(nm, term) {
    am <- ame_list[[nm]]
    if (is.null(am)) return(c(est = "", se = ""))
    row <- am[am$term == term, , drop = FALSE]
    if (nrow(row) == 0L) return(c(est = "", se = ""))
    c(est = sprintf("%.4f%s", row$estimate[1L], .stars(row$p.value[1L])),
      se  = sprintf("(%.4f)", row$std.error[1L]))
  }

  .crc <- function(term) {
    lab <- ET_DICT[term]; if (is.na(lab)) lab <- term
    r1  <- c(paste0("\\quad ", lab),
             unlist(lapply(col_names, function(nm)
               c(.hr_cell(nm, term)["est"], .ame_cell(nm, term)["est"]))))
    r2  <- c("",
             unlist(lapply(col_names, function(nm)
               c(.hr_cell(nm, term)["ci"], .ame_cell(nm, term)["se"]))))
    c(paste(r1, collapse = " & "), " \\\\",
      paste(r2, collapse = " & "), " \\\\")
  }

  fmt_n <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)
  .s2 <- function(v, a = "r") paste0("\\SetCell[c=2]{", a, "} ", v)

  header_c <- c(
    paste0(" & \\SetCell[c=2]{c} (1) & & \\SetCell[c=2]{c} (2) & \\\\"),
    "\\hline",
    " & HR & AME & HR & AME \\\\"
  )
  body_c <- c(
    .crc("sk_z"),
    "",
    .crc("bev_z"),
    .crc("chg_z"),
    "",
    .crc("state_gruene_z"),
    "",
    .crc("pers_z"),
    "\\hline",
    .crc("log_dens_z")
  )
  footer_c <- c(
    paste(c("Year/State FE",
            .s2("Yes","c"),"", .s2("Yes","c"),""),
          collapse = " & "), " \\\\",
    paste(c("Observations",
            .s2(fmt_n(stats[["(1)"]]$n),"c"),"",
            .s2(fmt_n(stats[["(2)"]]$n),"c"),""),
          collapse = " & "), " \\\\",
    paste(c("Onset events",
            .s2(fmt_n(ev[["(1)"]]),"c"),"",
            .s2(fmt_n(ev[["(2)"]]),"c"),""),
          collapse = " & "), " \\\\",
    paste(c("Pseudo-$R^2$",
            .s2(sprintf("%.3f", stats[["(1)"]]$pr2),"c"),"",
            .s2(sprintf("%.3f", stats[["(2)"]]$pr2),"c"),""),
          collapse = " & "), " \\\\"
  )
  write_longtblr(
    stem        = stem_hr_ame,
    caption     = "Discrete-Time Hazard: Broad Treatment, Hazard Ratios and Average Marginal Effects",
    label       = "tab:hazard_broad_hr_ame",
    note        = paste0(
      "Discrete-time hazard (cloglog link) of broad funding onset, AGS8 risk set 2016-2023. ",
      "AGS5 and year FE included. County-funded indicator excluded due to collinearity with ",
      "broad treatment. Hazard ratios (HR) computed as $\\exp(\\hat{\\beta})$, ",
      "95\\% CIs in brackets. Average marginal effects (AME) at sample averages reported ",
      "with AGS5-clustered SEs in parentheses. $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$."),
    colspec     = "l r l r l",
    header_rows = header_c,
    body_rows   = body_c,
    footer_rows = footer_c,
    resize      = TRUE
  )
  cat(sprintf("  Broad HR/AME: %s.tex\n", stem_hr_ame))

  ev_counts <- setNames(lapply(col_names, function(nm) ev[[nm]]), col_names)
  write_coef_longtblr(
    models       = models,
    stem         = stem_coef,
    caption      = "Discrete-Time Hazard: Broad Treatment, Raw Coefficients",
    label        = "tab:hazard_broad_coef",
    note         = paste0(
      "Raw coefficients of discrete-time hazard (cloglog link) of broad funding onset, ",
      "AGS8 risk set 2016-2023. County-funded indicator excluded due to collinearity with ",
      "broad treatment. AGS5 and year FE included. AGS5-clustered SEs in parentheses. ",
      "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$."),
    groups       = HAZARD_GROUPS_BROAD,
    var_labels   = ET_DICT,
    event_counts = ev_counts
  )
  cat(sprintf("  Broad coef: %s.tex\n", stem_coef))

  invisible(list(models = models, ph = ph_brd, ame_list = ame_list))
}

# -- Run ----------------------------------------------------------------------

ph_direct   <- .read_frame(file.path(data_final, "frame_hazard.csv"))
ph_broad    <- .read_frame(file.path(data_final, "frame_hazard_cov.csv"))

# Main paper table: HR/AME (main text) + raw coef (appendix)
run_hazard_main(
  ph_direct, "onset_direct",
  stem_hr_ame = file.path(out_dir, "hazard_hr_ame"),
  stem_coef   = file.path(out_dir, "hazard_coef")
)

# Appendix logit (risk-set censored): OR/AME + raw coef
run_logit_main(
  ph_direct, "onset_direct",
  stem_hr_ame = file.path(out_dir, "hazard_logit_hr_ame"),
  stem_coef   = file.path(out_dir, "hazard_logit_coef")
)

# Broad treatment (cloglog): HR/AME + raw coef
run_hazard_broad(
  ph_broad,
  stem_hr_ame = file.path(out_dir, "hazard_broad_hr_ame"),
  stem_coef   = file.path(out_dir, "hazard_broad_coef")
)

cat(sprintf("\nHazard outputs -> %s\n", out_dir))
