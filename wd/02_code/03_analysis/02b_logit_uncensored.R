# ─────────────────────────────────────────────────────────────────────────────
# 02b_logit_uncensored.R  — Absorbing-treatment logit robustness
#
# Plain logit with an ABSORBING treatment indicator (= 1 from first_treat
# onward) on the full AGS8 panel (no risk-set censoring). Completely
# independent of 02_hazard.R — different outcome, different sample logic.
#
# Treatment at AGS8 level, FE at AGS2 level: no collinearity.
# Absorbing dummy defined as: treat_absorb = 1{year >= first_treat}.
#
# Same 5-column spec grid and AGS2-clustered SEs as the main hazard table.
# Run for both direct and broad treatment definitions.
#
# Outputs (04_results/02b_logit_uncensored/):
#   tab_absorb_main_{direct,broad}.{tex,csv}   Table: absorbing logit coefs
#   tab_absorb_ame_{direct,broad}.{tex,csv}    AME table
#   tab_absorb_diag_{direct,broad}.csv         Events per AGS2/year
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
  stop("Cannot determine script path. Run as: Rscript 02b_logit_uncensored.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "04_results", "02b_logit_uncensored")
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
  d[, bev_z          := z(log1p(bev_win))]
  d[, chg_z          := z(log1p(chg_win))]
  d[, state_gruene_z := z(state_gruene_L1)]
  d[, fed_gruene_z   := z(fed_gruene_L1)]
  d[, eco_z          := z(eco_index_L1)]
  d[, pers_z         := z(log1p(pmax(n_vze_personal_L1, 0)))]
  d
}

# -- Formula grid (same 5-col layout as 02_hazard.R) -------------------------

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

# -- etable display -----------------------------------------------------------

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
  dict        = ET_DICT,
  keep_raw    = ET_KEEP,
  digits      = 3,
  se.below    = TRUE,
  signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.1),
  fitstat     = ~ n + pr2
)

.etable_both <- function(models, stem, note, ...) {
  args <- c(list(models), et_base, list(notes = note), list(...))
  do.call(etable, args)
  do.call(etable, c(args, list(
    file    = paste0(stem, ".tex"),
    replace = TRUE
  )))
  write_estimates_csv(models, paste0(stem, ".csv"))
}

# -- AME table ----------------------------------------------------------------

.run_ame <- function(models, stem) {
  cat("Computing AMEs...\n")
  ame_rows <- rbindlist(lapply(seq_along(models), function(i) {
    nm <- names(models)[i]
    am <- tryCatch(avg_slopes(models[[i]], vcov = FALSE),
                   error = function(e) { cat("  AME failed:", nm, "\n"); NULL })
    if (is.null(am)) return(NULL)
    as.data.table(am)[, .(model = nm, term, estimate)]
  }), fill = TRUE)
  if (nrow(ame_rows) == 0L) return(invisible(NULL))
  fwrite(ame_rows, paste0(stem, ".csv"))

  ame_rows[, cell := sprintf("%.4f", estimate)]
  term_order <- intersect(ET_KEEP, unique(ame_rows$term))
  col_order  <- names(models)
  est_wide   <- dcast(ame_rows[term %in% term_order],
                      term ~ model, value.var = "cell", fill = "")
  est_wide   <- est_wide[match(term_order, term)]
  n_c <- length(col_order)
  tex <- c(
    sprintf("\\begin{tabular}{l*{%d}{r}}", n_c), "\\hline\\hline",
    paste0(c("", col_order), collapse = " & "), " \\\\", "\\hline",
    vapply(seq_along(term_order), function(i) {
      lab <- ET_DICT[term_order[i]]
      if (is.na(lab)) lab <- term_order[i]
      paste0(paste(c(lab, unlist(est_wide[i, col_order, with = FALSE])),
                   collapse = " & "), " \\\\")
    }, character(1L)),
    "\\hline",
    sprintf("\\multicolumn{%d}{l}{\\small Point AMEs (logit, vcov=FALSE).}",
            n_c + 1L),
    "\\end{tabular}"
  )
  writeLines(tex, paste0(stem, ".tex"))
  cat("AME table written.\n")
  invisible(ame_rows)
}

# -- Core routine -------------------------------------------------------------

run_absorb_logit <- function(ph, first_treat_col, label) {
  cat(sprintf(
    "\n=== Absorbing logit | treat: %s | label: %s ===\n",
    first_treat_col, label
  ))

  # Build absorbing dummy: 1 from first_treat year onward
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

  # Diagnostics
  diag_dt <- rbindlist(list(
    ph[, .(n_obs = .N, n_treated_obs = sum(treat_absorb)),
       by = .(grp = AGS2)][, dim := "AGS2"],
    ph[, .(n_obs = .N, n_treated_obs = sum(treat_absorb)),
       by = .(grp = as.character(year))][, dim := "year"]
  ))
  fwrite(diag_dt,
         file.path(out_dir, sprintf("tab_absorb_diag_%s.csv", label)))

  fs <- make_formulas("treat_absorb")
  models <- lapply(fs, function(f)
    feglm(f, data = ph, family = binomial("logit"), cluster = ~AGS2))
  names(models) <- names(fs)

  note <- sprintf(
    paste0(
      "Absorbing-treatment logit on the full AGS8 panel (year $\\geq$ 2015, ",
      "no risk-set censoring). Outcome: \\texttt{treat\\_absorb} $= 1$ ",
      "for all years $\\geq$ \\texttt{%s} (absorbing), 0 otherwise. ",
      "Treatment at AGS8 level; FE at AGS2 level — no collinearity. ",
      "BEV stock and EV charger covariates winsorised at 99th pct before ",
      "log1p. Year + AGS2 FE; SEs clustered on AGS2. ",
      "$^{*}$ $p<0.1$, $^{**}$ $p<0.05$, $^{***}$ $p<0.01$."
    ),
    first_treat_col
  )

  .etable_both(
    models,
    file.path(out_dir, sprintf("tab_absorb_main_%s", label)),
    note,
    title = sprintf("Absorbing-treatment logit (%s)", label)
  )
  .run_ame(models,
           file.path(out_dir, sprintf("tab_absorb_ame_%s", label)))

  invisible(models)
}

# -- Run ----------------------------------------------------------------------

ph <- .read_frame(file.path(data_final, "frame_logit_full.csv"))

run_absorb_logit(ph, "first_treat_direct", "direct")
run_absorb_logit(ph, "first_treat_broad",  "broad")

cat(sprintf("\nAbsorbing logit outputs -> %s\n", out_dir))
