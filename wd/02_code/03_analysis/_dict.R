# Shared variable / label conventions for analysis scripts.
# Sourced from each NN_*.R via `source(file.path(code_dir, "03_analysis", "_dict.R"))`.

# Variable labels for etable / modelsummary
dict <- c(
  # IDs / FEs
  AGS8 = "AGS8 (Gemeinde)",
  AGS5 = "AGS5 (Kreis)",
  AGS2 = "AGS2 (Bundesland)",
  year = "Year",

  # Hazard channels (lagged, z-scored on estimation sample)
  log_dens_z      = "Log pop. density (z)",
  log_dens        = "Log pop. density",
  sk_z            = "log1p Steuerkraft (z, L1)",
  kk_z            = "Log Kaufkraft (z, L1)",
  bev_z           = "log1p BEV stock p100k (z, L1)",
  chg_z           = "log1p Charging pts p100k (z, L1)",
  pers_z          = "log1p Personnel VZE p100k (z, L1)",
  muni_gruene_z   = "Muni Grüne share (z, L1)",
  state_gruene_z  = "State Grüne share (z, L1)",
  fed_gruene_z    = "Fed Grüne share (z, L1)",
  eco_index_L1    = "EV ecosystem index (PCA, L1)",
  kreis_funded    = "Kreis-funded (strict past)",

  # Baseline snapshot z-scores (DiD covariates)
  kk_base_z    = "Baseline log Kaufkraft (z)",
  sk_base_z    = "Baseline log Steuerkraft (z)",
  dens_base_z  = "Baseline log pop. density (z)",
  green_base_z = "Baseline muni Grüne share (z)",
  bev_base_z   = "Baseline log1p BEV stock p100k (z)",
  chg_base_z   = "Baseline log1p Charging pts p100k (z)",

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
XFORMLA_CS <- ~ kk_base_z + sk_base_z + dens_base_z + green_base_z

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
