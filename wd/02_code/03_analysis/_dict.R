# Shared variable / label conventions for analysis scripts.
# Sourced from each NN_*.R via `source(file.path(code_dir, "03_analysis", "_dict.R"))`.

# Variable labels for etable / modelsummary
dict <- c(
  # IDs / FEs
  AGS8 = "AGS8 (Municipality)",
  AGS5 = "AGS5 (County)",
  AGS2 = "AGS2 (State)",
  year = "Year",

  # Hazard channels (lagged, z-scored on estimation sample)
  log_dens_z      = "Log pop. density (z)",
  log_dens        = "Log pop. density",
  sk_z            = "log1p Tax capacity (z, L1)",
  kk_z            = "Log Purchasing power (z, L1)",
  bev_z           = "log1p BEV stock p100k (z, L1)",
  chg_z           = "log1p Charging pts p100k (z, L1)",
  pers_z          = "log1p Personnel FTE p100k (z, L1)",
  muni_gruene_z   = "Muni Green share (z, L1)",
  state_gruene_z  = "State Green share (z, L1)",
  fed_gruene_z    = "Fed Green share (z, L1)",
  kreis_funded    = "County-funded (strict past)",

  # Baseline snapshot z-scores (DiD covariates)
  kk_base_z          = "Baseline log purchasing power (z)",
  sk_base_z          = "Baseline log tax capacity (z)",
  dens_base_z        = "Baseline log pop. density (z)",
  green_base_z       = "Baseline muni Green share (z)",
  state_green_base_z = "Baseline state Green share (z)",
  bev_base_z         = "Baseline log1p BEV stock p100k (z)",
  chg_base_z         = "Baseline log1p Charging pts p100k (z)",

  # Events
  onset_direct = "Direct-treatment onset",
  onset_broad  = "Coverage-event onset",

  # Outcomes
  bev_neuzulassungen_p100k = "BEV new registrations (p100k)",
  bev_corporate_p100k      = "BEV new registrations, corporate (p100k)",
  bev_private_p100k        = "BEV new registrations, private (p100k)",
  bev_stock_p100k          = "BEV stock (p100k)",
  ice_neuzulassungen_p100k = "ICE new registrations (p100k) — placebo"
)

# Outcome -> human label
OUTCOME_LABELS <- c(
  bev_neuzulassungen_p100k = "BEV new registrations per 100k population",
  bev_corporate_p100k      = "BEV new registrations (corporate) per 100k",
  bev_private_p100k        = "BEV new registrations (private) per 100k",
  bev_stock_p100k          = "BEV stock per 100k population",
  ice_neuzulassungen_p100k = "ICE new registrations per 100k (placebo)"
)

# Event-study horizon for plots / aggte. ES_MAX is the *display* default; each
# downstream script should compute its own data-driven cap (last outcome year
# minus earliest reachable cohort) via `es_max_data_driven()` below and pass it
# to aggte() / cap displayed ribbons.
ES_MIN <- -4L
ES_MAX <-  7L

# Stadtstaaten (n_vze_personal conflates municipal / Länder roles)
STADTSTAATEN <- c("02", "04", "11")

# BJS never-treated coding. didimputation v0.5.1 docs say `0`; we lock it here
# so 03/04/06 cannot drift apart. The A1 guard in 03_did_main.R verifies the
# untreated sample is non-degenerate; if it is, switch this to NA_real_ in one
# place and rerun.
BJS_NEVER <- 0  # numeric so fifelse stays double; flip to NA_real_ if A1 trips

# Canonical CS-dr xformla — used by 03, 04, 05, 06. Baseline (z'd) covariates.
# Note: bev_base_z and chg_base_z are dropped because (a) baseline 2014-16
# BEV/charge counts are mass-zero across the AGS8 cross-section, so they
# collapse to a near-constant column inside small (g,t) cells and trigger
# singular-matrix errors in the DR fit; (b) baseline outcome is mechanically
# related to the outcome family and is the wrong control for a flow.
XFORMLA_CS <- ~ sk_base_z + state_green_base_z + dens_base_z

# Data-driven upper horizon: last outcome year minus earliest reachable cohort.
# Pass the data.table and a non-NA outcome column.
es_max_data_driven <- function(dat, yname, gname_col = "gname_cs") {
  yr_max  <- dat[!is.na(get(yname)), max(year)]
  g_min   <- dat[get(gname_col) > 0L, min(get(gname_col))]
  if (!is.finite(yr_max) || !is.finite(g_min)) return(ES_MAX)
  max(1L, yr_max - g_min)
}

# Standard z helper, used after every sample restriction so each model's
# regressors are standardised on its own estimation sample.
z <- function(x) as.numeric(scale(x))

# Population-weighted quantile cuts. Returns the q-th sample quantiles of x
# weighted by w (NA-safe).
wq <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x  <- x[ok]; w <- w[ok]
  o  <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1L]], numeric(1))
}

# CSV+TeX twin writer for fixest etable. Writes a .tex and a .csv with the
# same coefficient/SE numbers, plus a manifest row.
write_estimates_csv <- function(models, file) {
  rows <- lapply(seq_along(models), function(i) {
    m  <- models[[i]]
    nm <- names(models)[i]
    co <- coef(m); se <- sqrt(diag(vcov(m)))
    data.table(
      model = nm,
      term  = names(co),
      estimate = as.numeric(co),
      se       = as.numeric(se[names(co)]),
      n        = nobs(m)
    )
  })
  dt <- data.table::rbindlist(rows)
  dt[, ci_lo := estimate - 1.96 * se]
  dt[, ci_hi := estimate + 1.96 * se]
  data.table::fwrite(dt, file)
  invisible(dt)
}

# Resolve project root from the running script path. Drop-in for the
# argv/self_flag block at the top of every analysis script.
resolve_root <- function() {
  argv      <- commandArgs(trailingOnly = FALSE)
  self_flag <- grep("--file=", argv, value = TRUE)
  self <- if (length(self_flag)) {
    normalizePath(sub("--file=", "", self_flag))
  } else if (
    requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()
  ) {
    normalizePath(rstudioapi::getSourceEditorContext()$path)
  } else {
    stop("Cannot determine script path. Run as: Rscript <script>.R")
  }
  dirname(dirname(dirname(self)))
}

# Standard output directory for a script: 04_results/<script_stem>/
results_dir <- function(root, stem) {
  d <- file.path(root, "04_results", stem)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

# ── Table writers ──────────────────────────────────────────────────────────────

# Core longtblr writer (tabularray package). All regression/effect and
# descriptive tables use this format. `note{}` (empty key) produces an
# unnumbered footnote at the bottom of the float.
write_longtblr <- function(stem, caption, label, note, colspec,
                            header_rows, body_rows, footer_rows,
                            rowsep = "-3pt") {
  tex <- c(
    "\\begin{longtblr}[",
    paste0("    caption = {", caption, "},"),
    paste0("    label = {", label, "},"),
    paste0("    note{} = {\\small ", note, "},"),
    "]{",
    paste0("  colspec = {", colspec, "},"),
    paste0("  rowsep  = ", rowsep, ","),
    "}",
    "\\hline\\hline",
    header_rows,
    "\\hline",
    body_rows,
    "\\hline",
    footer_rows,
    "\\hline\\hline",
    "\\end{longtblr}"
  )
  writeLines(tex, paste0(stem, ".tex"))
  invisible(tex)
}

# Coefficient table writer for fixest models. Builds body from a `groups` list
# (names = section labels used as blank-line separators, values = char vector
# of term names; NULL value inserts a \\hline separator).
# event_counts: named list model_name -> integer, for "Onset events" footer row.
write_coef_longtblr <- function(
  models, stem, caption, label, note, groups, var_labels,
  event_counts = NULL, show_pr2 = TRUE
) {
  col_names <- names(models)
  n_spec    <- length(col_names)

  .stars <- function(p) ifelse(is.na(p), "",
    ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", ""))))

  .cell <- function(nm, term) {
    ct <- tryCatch(coeftable(models[[nm]]), error = function(e) NULL)
    if (is.null(ct) || !term %in% rownames(ct)) return(c(est = "", se = ""))
    b <- ct[term, "Estimate"]; s <- ct[term, "Std. Error"]; p <- ct[term, ncol(ct)]
    c(est = sprintf("%.3f%s", b, .stars(p)), se = sprintf("(%.3f)", s))
  }

  .coef_rows <- function(term) {
    lab   <- var_labels[term]; if (is.na(lab)) lab <- term
    cells <- lapply(col_names, .cell, term = term)
    r1    <- c(paste0("\\quad ", lab), vapply(cells, `[[`, character(1), "est"))
    r2    <- c("",                     vapply(cells, `[[`, character(1), "se"))
    c(paste(r1, collapse = " & "), " \\\\",
      paste(r2, collapse = " & "), " \\\\")
  }

  body_rows <- c()
  for (i in seq_along(groups)) {
    terms <- groups[[i]]
    if (is.null(terms)) {
      body_rows <- c(body_rows, "\\hline")
    } else {
      if (i > 1L && !is.null(groups[[i - 1L]])) body_rows <- c(body_rows, "")
      body_rows <- c(body_rows, unlist(lapply(terms, .coef_rows)))
    }
  }

  header_rows <- paste0(paste(c("", col_names), collapse = " & "), " \\\\")

  fmt_n <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)
  n_vals <- setNames(
    vapply(col_names, function(nm)
      tryCatch(nobs(models[[nm]]), error = function(e) NA_integer_), integer(1)),
    col_names)
  pr2_vals <- if (show_pr2)
    setNames(vapply(col_names, function(nm)
      tryCatch(as.numeric(r2(models[[nm]], "pr2")), error = function(e) NA_real_),
      numeric(1)), col_names)
  else NULL

  footer_rows <- c(
    paste0(paste(c("Year/State FE", rep("Yes", n_spec)), collapse = " & "), " \\\\"),
    paste0(paste(c("Observations",  fmt_n(n_vals)),       collapse = " & "), " \\\\")
  )
  if (!is.null(event_counts)) {
    ev_vals <- vapply(col_names, function(nm) fmt_n(event_counts[[nm]]), character(1))
    footer_rows <- c(footer_rows,
      paste0(paste(c("Onset events", ev_vals), collapse = " & "), " \\\\"))
  }
  if (!is.null(pr2_vals)) {
    footer_rows <- c(footer_rows,
      paste0(paste(c("Pseudo-$R^2$", sprintf("%.3f", pr2_vals)), collapse = " & "), " \\\\"))
  }

  write_longtblr(stem = stem, caption = caption, label = label, note = note,
                 colspec = sprintf("l *{%d}{r}", n_spec),
                 header_rows = header_rows,
                 body_rows   = body_rows,
                 footer_rows = footer_rows)
  write_estimates_csv(models, paste0(stem, ".csv"))
  invisible(models)
}
