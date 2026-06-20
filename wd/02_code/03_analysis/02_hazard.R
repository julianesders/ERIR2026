# ─────────────────────────────────────────────────────────────────────────────
# 02_hazard.R   — Part (i): drivers of EMK onset
#
# Discrete-time hazard on the AGS8 risk set (year >= 2015).
# All specs: year + AGS2 FE, SEs clustered on AGS2.
# Primary event: onset_direct (frame_hazard.csv).
#
# 5-column specification grid (log_dens_z universal, listed last):
#   (1) sk + bev + chg + state_gruene + log_dens    ← baseline
#   (2) (1) + pers                                  ← no-Stadtstaaten
#   (3) sk + eco + pers + state_gruene + log_dens
#   (4) sk + bev + chg + fed_gruene + log_dens
#   (5) sk + bev + chg + kreis_funded + log_dens
#
# Tables produced in 04_results/02_hazard/:
#   tab_hazard_main.{tex,csv}              Table 1: cloglog, direct (main)
#   tab_hazard_appendix.{tex,csv}          Table 2: logit, direct, censored
#   tab_hazard_ame_main.{tex,csv}          AME for Table 1
#   tab_hazard_robust_uncensored.{tex,csv} logit on frame_logit_full
#   tab_hazard_robust_penalized.{tex,csv}  brglm2 penalised logit
#   tab_hazard_robust_broad.{tex,csv}      cloglog, broad treatment
#   tab_hazard_diag_{main,appendix}.csv    events per AGS2/year
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
out_dir    <- file.path(root, "04_results", "02_hazard")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))

# -- I/O ----------------------------------------------------------------------

.read_frame <- function(p) fread(p,
  colClasses = list(character = c("AGS8", "AGS5", "AGS2", "treat_type")))

# -- Regressor construction ---------------------------------------------------
# BEV stock and charger covariates are winsorised at 99th pct on the
# estimation sample before log1p; this mirrors WINSOR_Q = 0.99 applied to
# the outcome variables in 00_prep_analysis.R.

# fed_gruene_L1 has ~13k NAs and is used only in col (4); excluding it from
# the shared filter keeps the common sample at ~75k instead of ~62k.
# feglm handles per-spec NAs automatically.
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
  d[, bev_z          := z(log1p(bev_win))]
  d[, chg_z          := z(log1p(chg_win))]
  d[, state_gruene_z := z(state_gruene_L1)]
  d[, fed_gruene_z   := z(fed_gruene_L1)]
  d[, eco_z          := z(eco_index_L1)]
  d[, pers_z         := z(log1p(pmax(n_vze_personal_L1, 0)))]
  d
}

# -- Formula grid -------------------------------------------------------------

# 5-column grid: log_dens_z is a universal control, listed last in every spec.
# Old cols (1) [no density] and (2) [+density] collapsed to a single baseline
# after density was promoted to a universal control in all specs.
make_formulas <- function(y) {
  fe <- "| year + AGS2"
  list(
    `(1)` = as.formula(sprintf(
      "%s ~ sk_z + bev_z + chg_z + state_gruene_z + log_dens_z %s", y, fe)),
    `(2)` = as.formula(sprintf(
      "%s ~ sk_z + bev_z + chg_z + pers_z + state_gruene_z + log_dens_z %s",
      y, fe)),
    `(3)` = as.formula(sprintf(
      "%s ~ sk_z + eco_z + pers_z + state_gruene_z + log_dens_z %s", y, fe)),
    `(4)` = as.formula(sprintf(
      "%s ~ sk_z + bev_z + chg_z + fed_gruene_z + log_dens_z %s", y, fe)),
    `(5)` = as.formula(sprintf(
      "%s ~ sk_z + bev_z + chg_z + kreis_funded + log_dens_z %s", y, fe))
  )
}

# -- etable display settings --------------------------------------------------

ET_DICT <- c(dict,
  sk_z           = "log1p Steuerkraft p.c. (z, L1)",
  bev_z          = "BEV stock p100k (z, L1, 99\\%)",
  chg_z          = "EV chargers p100k (z, L1, 99\\%)",
  eco_z          = "EV ecosystem index (z, L1)",
  pers_z         = "log1p Personnel p.c. (z, L1)",
  state_gruene_z = "Gr\\\"une vote share, state (z, L1)",
  fed_gruene_z   = "Gr\\\"une vote share, federal (z, L1)",
  kreis_funded   = "Kreis funded",
  log_dens_z     = "Log pop. density (z)"
)

ET_KEEP <- c("sk_z", "bev_z", "chg_z", "eco_z", "pers_z",
             "state_gruene_z", "fed_gruene_z", "kreis_funded", "log_dens_z")

et_base <- list(
  dict         = ET_DICT,
  keep_raw     = ET_KEEP,
  digits       = 3,
  se.below     = TRUE,
  signif.code  = c("***" = 0.01, "**" = 0.05, "*" = 0.1),
  fitstat      = ~ n + pr2
)

.etable_both <- function(models, stem, note, ...) {
  args <- c(list(models), et_base, list(notes = note), list(...))
  do.call(etable, args)                                   # console
  do.call(etable, c(args, list(                           # tex file
    file    = paste0(stem, ".tex"),
    replace = TRUE
  )))
  write_estimates_csv(models, paste0(stem, ".csv"))       # csv file
}

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

# -- Main table runner --------------------------------------------------------

run_table <- function(ph, event_col, link, stem, title_tag, note_extra = "") {
  cat(sprintf("\n=== %s | event: %s | link: %s ===\n",
              stem, event_col, link$family))
  ph <- ph[complete.cases(ph[, ..BASE_COLS])]
  ph <- zfit(ph)
  cat(sprintf("  N = %d obs | %d AGS8 | %d events (%.2f%%)\n",
              nrow(ph), uniqueN(ph$AGS8),
              sum(ph[[event_col]]), 100 * mean(ph[[event_col]])))
  .diag(ph, event_col, basename(stem))
  fs <- make_formulas(event_col)
  models <- lapply(fs, function(f)
    feglm(f, data = ph, family = link, cluster = ~AGS2))
  names(models) <- names(fs)
  note <- sprintf(
    paste0(
      "Discrete-time hazard on the AGS8 risk set (year \\geq 2015). ",
      "Event: \\texttt{%s}. ",
      "Cols (3)--(4): include log1p municipal personnel p.c. (z, L1). ",
      "BEV stock and EV charger covariates winsorised at the 99th pct ",
      "before log1p transformation. ",
      "All columns: year + AGS2 fixed effects; SEs clustered on AGS2. ",
      "%s",
      "$^{*}$ $p<0.1$, $^{**}$ $p<0.05$, $^{***}$ $p<0.01$."
    ),
    event_col,
    if (nchar(note_extra)) paste0(note_extra, " ") else ""
  )
  .etable_both(models, file.path(out_dir, stem), note,
               title = title_tag)
  invisible(list(models = models, ph = ph))
}

# -- AME table ----------------------------------------------------------------

run_ame <- function(models, stem) {
  cat(sprintf("\nComputing AMEs for %s...\n", stem))
  ame_rows <- rbindlist(lapply(seq_along(models), function(i) {
    nm <- names(models)[i]
    am <- tryCatch(
      avg_slopes(models[[i]], vcov = FALSE),
      error = function(e) { cat("  AME failed:", nm, "\n"); NULL }
    )
    if (is.null(am)) return(NULL)
    dt <- as.data.table(am)
    dt[, .(model = nm, term, estimate)]
  }), fill = TRUE)

  if (nrow(ame_rows) == 0L) return(invisible(NULL))
  fwrite(ame_rows, file.path(out_dir, paste0(stem, ".csv")))

  # Wide tex: point AMEs only (vcov=FALSE means no SEs available)
  ame_rows[, cell := sprintf("%.4f", estimate)]
  term_order <- intersect(ET_KEEP, unique(ame_rows$term))
  col_order  <- names(models)
  est_wide   <- dcast(ame_rows[term %in% term_order],
                      term ~ model, value.var = "cell", fill = "")
  est_wide   <- est_wide[match(term_order, term)]

  n_cols <- length(col_order)
  tex_lines <- c(
    sprintf("\\begin{tabular}{l*{%d}{r}}", n_cols),
    "\\hline\\hline",
    paste0(c("", col_order), collapse = " & "), " \\\\",
    "\\hline",
    vapply(seq_along(term_order), function(i) {
      lab <- ET_DICT[term_order[i]]
      if (is.na(lab)) lab <- term_order[i]
      row <- unlist(est_wide[i, col_order, with = FALSE])
      paste0(paste(c(lab, row), collapse = " & "), " \\\\")
    }, character(1L)),
    "\\hline",
    sprintf(
      "\\multicolumn{%d}{l}{\\small Point AMEs (cloglog, vcov=FALSE; SEs not available with FE).}",
      n_cols + 1L
    ),
    "\\end{tabular}"
  )
  writeLines(tex_lines, file.path(out_dir, paste0(stem, ".tex")))
  cat("AME table written.\n")
  invisible(ame_rows)
}

# -- Robustness ---------------------------------------------------------------

run_robustness <- function(ph_full, ph_broad) {
  cat("\n========== Robustness ==========\n")

  # (a) Uncensored logit
  cat("\n--- (a) Uncensored logit (frame_logit_full) ---\n")
  ph_u <- ph_full[complete.cases(ph_full[, ..BASE_COLS])]
  ph_u <- zfit(ph_u)
  cat(sprintf("  N = %d | %d AGS8 | %d onset_direct (%.2f%%)\n",
              nrow(ph_u), uniqueN(ph_u$AGS8),
              sum(ph_u$onset_direct), 100 * mean(ph_u$onset_direct)))
  fs_u <- make_formulas("onset_direct")
  m_unc <- lapply(fs_u, function(f)
    feglm(f, data = ph_u, family = binomial("logit"), cluster = ~AGS2))
  names(m_unc) <- names(fs_u)
  note_unc <- paste0(
    "Robustness: uncensored logit. Post-onset AGS8-years are retained; ",
    "the onset indicator is mechanically zero after first onset. ",
    "Otherwise identical to Table 1. year + AGS2 FE; SEs clustered on AGS2. ",
    "$^{*}$ $p<0.1$, $^{**}$ $p<0.05$, $^{***}$ $p<0.01$."
  )
  .etable_both(m_unc,
               file.path(out_dir, "tab_hazard_robust_uncensored"),
               note_unc,
               title = "Robustness: uncensored logit (onset\\_direct)")

  # (b) Penalised logit (brglm2), all 6 specs
  if (requireNamespace("brglm2", quietly = TRUE)) {
    cat("\n--- (b) Penalised logit (brglm2), all 6 specs ---\n")
    ph_br <- copy(ph_full[complete.cases(ph_full[, ..BASE_COLS])])
    ph_br <- zfit(ph_br)
    ph_br[, year_f := factor(year)][, ags2_f := factor(AGS2)]

    spec_vars <- list(
      `(1)` = "sk_z + bev_z + chg_z + state_gruene_z + log_dens_z",
      `(2)` = "sk_z + bev_z + chg_z + pers_z + state_gruene_z + log_dens_z",
      `(3)` = "sk_z + eco_z + pers_z + state_gruene_z + log_dens_z",
      `(4)` = "sk_z + bev_z + chg_z + fed_gruene_z + log_dens_z",
      `(5)` = "sk_z + bev_z + chg_z + kreis_funded + log_dens_z"
    )
    br_rows <- rbindlist(lapply(names(spec_vars), function(nm) {
      f_br <- as.formula(sprintf(
        "onset_direct ~ %s + year_f + ags2_f", spec_vars[[nm]]))
      fit <- tryCatch(
        glm(f_br, data = as.data.frame(ph_br),
            family = binomial("logit"),
            method = brglm2::brglm_fit),
        error = function(e) {
          cat(sprintf("  brglm2 failed spec %s: %s\n", nm, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(fit)) return(NULL)
      co   <- coef(fit); se <- sqrt(diag(vcov(fit)))
      keep <- !grepl("^(year_f|ags2_f|\\(Intercept\\))", names(co))
      data.table(model = nm, term = names(co)[keep],
                 estimate = co[keep], se = se[keep], n = nobs(fit))
    }), fill = TRUE)

    if (nrow(br_rows) > 0L) {
      fwrite(br_rows,
             file.path(out_dir, "tab_hazard_robust_penalized.csv"))
      # Pivot to wide for tex display
      br_rows[, cell := sprintf("%.3f", estimate)]
      br_wide <- dcast(br_rows, term ~ model, value.var = "cell", fill = "")
      se_rows <- copy(br_rows)[, cell := sprintf("(%.3f)", se)]
      se_wide <- dcast(se_rows, term ~ model, value.var = "cell", fill = "")
      term_ord <- intersect(ET_KEEP, br_wide$term)
      br_wide  <- br_wide[match(term_ord, term)]
      se_wide  <- se_wide[match(term_ord, term)]
      col_ord  <- names(spec_vars)
      n_c <- length(col_ord)
      tex_br <- c(
        sprintf("\\begin{tabular}{l*{%d}{r}}", n_c),
        "\\hline\\hline",
        paste0(c("", col_ord), collapse = " & "), " \\\\",
        "\\hline",
        do.call(c, lapply(seq_along(term_ord), function(i) {
          lab <- ET_DICT[term_ord[i]]; if (is.na(lab)) lab <- term_ord[i]
          e_r <- unlist(br_wide[i, col_ord, with = FALSE])
          s_r <- unlist(se_wide[i, col_ord, with = FALSE])
          c(paste0(c(lab, e_r), collapse=" & "), " \\\\",
            paste0(c("",  s_r), collapse=" & "), " \\\\")
        })),
        "\\hline",
        sprintf(
          "\\multicolumn{%d}{l}{\\small Penalised logit (brglm2, Firth correction). No cluster correction.}",
          n_c + 1L),
        "\\end{tabular}"
      )
      writeLines(tex_br,
                 file.path(out_dir, "tab_hazard_robust_penalized.tex"))
      cat("Penalised logit table written.\n")
    }
  } else {
    cat("brglm2 not installed — skipping penalised logit.\n")
  }

  # (c) Broad treatment (cloglog, censored)
  cat("\n--- (c) Broad treatment cloglog (frame_hazard_cov) ---\n")
  ph_brd <- ph_broad[complete.cases(ph_broad[, ..BASE_COLS])]
  ph_brd <- zfit(ph_brd)
  cat(sprintf("  N = %d | %d AGS8 | %d onset_broad (%.2f%%)\n",
              nrow(ph_brd), uniqueN(ph_brd$AGS8),
              sum(ph_brd$onset_broad), 100 * mean(ph_brd$onset_broad)))
  fs_brd <- make_formulas("onset_broad")
  m_brd <- lapply(fs_brd, function(f)
    feglm(f, data = ph_brd, family = binomial("cloglog"), cluster = ~AGS2))
  names(m_brd) <- names(fs_brd)
  note_brd <- paste0(
    "Robustness: broad treatment (direct OR Kreis-broadcast coverage). ",
    "Otherwise identical to Table 1. ",
    "year + AGS2 FE; SEs clustered on AGS2. ",
    "$^{*}$ $p<0.1$, $^{**}$ $p<0.05$, $^{***}$ $p<0.01$."
  )
  .etable_both(m_brd,
               file.path(out_dir, "tab_hazard_robust_broad"),
               note_brd,
               title = "Robustness: broad treatment cloglog (onset\\_broad)")
}

# -- Run ----------------------------------------------------------------------

ph_direct  <- .read_frame(file.path(data_final, "frame_hazard.csv"))
ph_broad   <- .read_frame(file.path(data_final, "frame_hazard_cov.csv"))
ph_full    <- .read_frame(file.path(data_final, "frame_logit_full.csv"))
ph_cov_full <- .read_frame(file.path(data_final, "frame_logit_cov_full.csv"))

# Table 1: cloglog, direct, censored
res_main <- run_table(
  ph_direct, "onset_direct", binomial("cloglog"),
  stem      = "tab_hazard_main",
  title_tag = "Discrete-time hazard (cloglog), direct treatment"
)

# AME for Table 1
run_ame(res_main$models, "tab_hazard_ame_main")

# Table 2 (appendix): logit, direct, censored
run_table(
  ph_direct, "onset_direct", binomial("logit"),
  stem      = "tab_hazard_appendix",
  title_tag = "Appendix: logit (censored), direct treatment",
  note_extra = "Appendix table: logit link; risk-set censoring retained."
)

# Table 3+: robustness
run_robustness(ph_full, ph_broad)

cat(sprintf("\nHazard outputs -> %s\n", out_dir))
