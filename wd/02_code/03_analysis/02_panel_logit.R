library(data.table)
library(fixest)

# -- Paths ---------------------------------------------------------------------

argv      <- commandArgs(trailingOnly = FALSE)
self_flag <- grep("--file=", argv, value = TRUE)
self <- if (length(self_flag)) {
  normalizePath(sub("--file=", "", self_flag))
} else if (
  requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()
) {
  normalizePath(rstudioapi::getSourceEditorContext()$path)
} else {
  stop("Cannot determine script path. Run as: Rscript 02_panel_logit.R")
}
root       <- dirname(dirname(dirname(self)))
data_int   <- file.path(root, "01_data", "02_intermediate")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "04_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -- Load and prepare ----------------------------------------------------------

panel <- fread(
  file.path(data_final, "emk_inkar_panel_ags8.csv"),
  colClasses = list(character = c("AGS8", "AGS5", "AGS2"))
)
panel <- panel[year >= 2015]

# Risk-set: units exit after their first onset year.
# Never-funded units are right-censored (all-zero rows throughout).
panel[,
  first_treat := if (any(emk_absorbing == 1L))
    min(year[emk_absorbing == 1L])
  else
    .Machine$integer.max,
  by = AGS8
]
ph <- panel[year <= first_treat][, first_treat := NULL]
panel[, first_treat := NULL]

# -- Design diagnostic: direct AGS8 assignment vs Kreis broadcast -------------

emk <- fread(
  file.path(data_int, "emk", "emk_ags_matched.csv"),
  colClasses = list(character = c("AGS8", "AGS5")),
  na.strings  = c("", "NA")
)
n_direct <- sum(!is.na(emk$AGS8))
cat(sprintf(
  "Projects: %d/%d directly AGS8-assigned (%d%%), %d/%d Kreis-broadcast (%d%%)\n",
  n_direct,              nrow(emk), round(100 * n_direct / nrow(emk)),
  nrow(emk) - n_direct,  nrow(emk), round(100 * (1 - n_direct / nrow(emk)))
))
cat(sprintf(
  "Risk set: %d obs | %d AGS8 | %d onsets | onset rate %.1f%%\n",
  nrow(ph), ph[, uniqueN(AGS8)], ph[, sum(emk_absorbing)],
  100 * ph[, mean(emk_absorbing)]
))

# -- Separation check under AGS5 FE -------------------------------------------

n_ags5_zero <- ph[, .(ever = max(emk_absorbing)), by = AGS5][ever == 0L, .N]
cat(sprintf(
  "AGS5 groups with no treated AGS8: %d of %d (dropped under AGS5 FE)\n",
  n_ags5_zero, ph[, uniqueN(AGS5)]
))

# -- Standardise continuous regressors on the estimation sample ---------------
# sk_z/sk_sq_z: build the square from the z-scored level so the polynomial
# is internally consistent. eco_index already standardised via PCA; no rescale.

z <- function(x) as.numeric(scale(x))

ph[, log_dens_z := z(log_pop_dens)]
ph[, pendler_z  := z(q_pendlersaldo)]
ph[, sk_z       := z(q_gest_bev_L1)]
ph[, sk_sq_z    := sk_z^2]
ph[, bev_z      := z(bev_stock_p100k_L1)]
ph[, chg_z      := z(ev_chargepoints_p100k_L1)]

sk_mean <- mean(ph$q_gest_bev_L1, na.rm = TRUE)
sk_sd   <- sd(ph$q_gest_bev_L1,   na.rm = TRUE)

# -- Stadtstaaten exclusion for personnel spec --------------------------------
# Hamburg (02), Bremen (04), Berlin (11): personnel conflates municipal/state
# roles, incomparable to Flachenlander. Baseline also run on this restricted
# sample so the only difference from the personnel model is the added regressor.

STADTSTAATEN <- c("02", "04", "11")
ph_ns <- ph[!(AGS2 %in% STADTSTAATEN)]
ph_ns[, sk_z    := z(q_gest_bev_L1)]   # re-z-score on restricted sample
ph_ns[, sk_sq_z := sk_z^2]
ph_ns[, pers_z  := z(n_vze_personal_L1)]
ph_ns[, log_dens_z := z(log_pop_dens)]
ph_ns[, pendler_z  := z(q_pendlersaldo)]
ph_ns[, bev_z      := z(bev_stock_p100k_L1)]
ph_ns[, chg_z      := z(ev_chargepoints_p100k_L1)]

sk_mean_ns <- mean(ph_ns$q_gest_bev_L1, na.rm = TRUE)
sk_sd_ns   <- sd(ph_ns$q_gest_bev_L1,   na.rm = TRUE)

cat(sprintf(
  "No-Stadtstaaten: %d obs dropped | %d obs | %d onsets\n",
  nrow(ph) - nrow(ph_ns), nrow(ph_ns), ph_ns[, sum(emk_absorbing)]
))

# -- NA coverage diagnostic ---------------------------------------------------
# Fail-fast: report coverage of every source variable before estimation so that
# all-NA predictors (e.g. eco_index when KBA Bestand is absent) are diagnosed
# without scrolling through individual feglm error messages.

cat("\nSource variable coverage in risk set (before z-scoring):\n")
src_vars <- c(
  "log_pop_dens", "q_pendlersaldo", "muni_gruene_L1",
  "eco_index_L1", "q_gest_bev_L1",
  "bev_stock_p100k_L1", "ev_chargepoints_p100k_L1", "n_vze_personal_L1"
)
for (v in src_vars) {
  if (v %in% names(ph)) {
    pct <- 100 * mean(!is.na(ph[[v]]))
    cat(sprintf("  %-30s: %5.1f%% non-NA\n", v, pct))
  } else {
    cat(sprintf("  %-30s: MISSING COLUMN -- panel may need regeneration\n", v))
  }
}
cat("\n")

# -- Specifications -----------------------------------------------------------
# year FE: discrete-time baseline hazard (non-parametric).
# AGS2 FE: identifies off cross-Kreis variation (16 Bundeslander).
# AGS5 FE: within-Kreis variation; absorbs any AGS5-constant regressors.
#   n_vze_personal is measured at AGS5 -> collinear with AGS5 FE, AGS2-only.
# Clustering on AGS5: coarsest level at which treatment is assigned.

f_base_a2 <- emk_absorbing ~
  log_dens_z + pendler_z + muni_gruene_L1 +
  eco_index_L1 + sk_z + sk_sq_z | year + AGS2

f_base_a5 <- emk_absorbing ~
  log_dens_z + pendler_z + muni_gruene_L1 +
  eco_index_L1 + sk_z + sk_sq_z | year + AGS5

f_comp_a2 <- emk_absorbing ~
  log_dens_z + pendler_z + muni_gruene_L1 +
  bev_z + chg_z + sk_z + sk_sq_z | year + AGS2

f_comp_a5 <- emk_absorbing ~
  log_dens_z + pendler_z + muni_gruene_L1 +
  bev_z + chg_z + sk_z + sk_sq_z | year + AGS5

f_pers_a2 <- emk_absorbing ~
  log_dens_z + pendler_z + muni_gruene_L1 +
  eco_index_L1 + sk_z + sk_sq_z + pers_z | year + AGS2

# -- Estimation ----------------------------------------------------------------

cll <- binomial("cloglog")   # primary: discrete-time proportional hazards
lgt <- binomial("logit")     # robustness

m_cll_base_a2 <- feglm(f_base_a2, data = ph,    family = cll, cluster = ~AGS5)
m_cll_base_a5 <- feglm(f_base_a5, data = ph,    family = cll, cluster = ~AGS5)
m_cll_comp_a2 <- feglm(f_comp_a2, data = ph,    family = cll, cluster = ~AGS5)
m_cll_comp_a5 <- feglm(f_comp_a5, data = ph,    family = cll, cluster = ~AGS5)
m_cll_base_ns <- feglm(f_base_a2, data = ph_ns, family = cll, cluster = ~AGS5)
m_cll_pers_a2 <- feglm(f_pers_a2, data = ph_ns, family = cll, cluster = ~AGS5)
m_lgt_base_a2 <- feglm(f_base_a2, data = ph,    family = lgt, cluster = ~AGS5)
m_lgt_base_a5 <- feglm(f_base_a5, data = ph,    family = lgt, cluster = ~AGS5)

# -- Diagnostics ---------------------------------------------------------------

cat("\nObservations used per model:\n")
for (pair in list(
  list(m_cll_base_a2, "cloglog base      AGS2"),
  list(m_cll_base_a5, "cloglog base      AGS5"),
  list(m_cll_comp_a2, "cloglog comp      AGS2"),
  list(m_cll_comp_a5, "cloglog comp      AGS5"),
  list(m_cll_base_ns, "cloglog base      AGS2 (no-Stadtstaaten)"),
  list(m_cll_pers_a2, "cloglog base+pers AGS2 (no-Stadtstaaten)"),
  list(m_lgt_base_a2, "logit   base      AGS2"),
  list(m_lgt_base_a5, "logit   base      AGS5")
)) {
  cat(sprintf("  %-45s n = %d\n", pair[[2]], nobs(pair[[1]])))
}

# -- Steuerkraft turning point ------------------------------------------------
# tp_z = -b_sk / (2 * b_sk_sq); back-convert to original units.
# If outside data range, the hump hypothesis is not supported in-sample.

tp <- function(m, mu, sig) {
  b <- coef(m)
  if (!all(c("sk_z", "sk_sq_z") %in% names(b))) return(NA_real_)
  (-b["sk_z"] / (2 * b["sk_sq_z"])) * sig + mu
}

cat(sprintf("\nSteuerkraft turning points (EUR/capita):\n"))
cat(sprintf("  cloglog base AGS2 (full):             %.0f\n",
  tp(m_cll_base_a2, sk_mean,    sk_sd)))
cat(sprintf("  cloglog base AGS2 (no-Stadtstaaten):  %.0f\n",
  tp(m_cll_base_ns, sk_mean_ns, sk_sd_ns)))
cat(sprintf("  Data range in estimation sample:      [%.0f, %.0f]\n",
  quantile(ph$q_gest_bev_L1, 0.01, na.rm = TRUE),
  quantile(ph$q_gest_bev_L1, 0.99, na.rm = TRUE)))

# -- Output -------------------------------------------------------------------
# Column groups: cloglog AGS2 | cloglog AGS5 | logit robustness (AGS2+AGS5).
# * = no-Stadtstaaten sample. Personnel spec is AGS2-only (AGS5 FE would
#   absorb n_vze_personal which is constant within Kreis).

dict <- c(
  log_dens_z   = "Log pop. density (z)",
  pendler_z    = "Commuter balance (z)",
  muni_gruene_L1 = "Mun. Gruene vote share (L1)",
  eco_index_L1 = "EV ecosystem index (L1)",
  sk_z         = "Steuerkraft (z, L1)",
  sk_sq_z      = "Steuerkraft sq. (z, L1)",
  bev_z        = "BEV stock p100k (z, L1)",
  chg_z        = "Charging pts p100k (z, L1)",
  pers_z       = "Personnel VZE (z, L1)"
)

models  <- list(m_cll_base_a2, m_cll_comp_a2, m_cll_base_ns, m_cll_pers_a2,
                m_cll_base_a5, m_cll_comp_a5,
                m_lgt_base_a2, m_lgt_base_a5)
headers <- c("Base", "Comp", "Base*", "Base*+Pers",
             "Base", "Comp", "Base", "Base")

do.call(etable, c(models, list(headers = headers, depvar = FALSE,
                               digits = 4, dict = dict)))

do.call(etable, c(models, list(headers = headers, depvar = FALSE,
                               digits = 4, dict = dict,
                               file    = file.path(out_dir, "tab_treatment_onset.tex"),
                               replace = TRUE)))

cat(sprintf("\nTeX table written to: %s\n",
  file.path(out_dir, "tab_treatment_onset.tex")))
