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
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "04_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -- Load and prepare ----------------------------------------------------------

panel <- fread(
  file.path(data_final, "emk_inkar_panel_ags8.csv"),
  colClasses = list(character = c("AGS8", "AGS5", "AGS2"))
)
panel <- panel[year >= 2015]

# Risk set: units exit after their first onset year; never-funded are censored.
panel[,
  first_treat := if (any(emk_absorbing == 1L))
    min(year[emk_absorbing == 1L])
  else
    .Machine$integer.max,
  by = AGS8
]
ph <- panel[year <= first_treat][, first_treat := NULL]
panel[, first_treat := NULL]

cat(sprintf(
  "Risk set: %d obs | %d AGS8 | %d onsets | onset rate %.1f%%\n",
  nrow(ph), ph[, uniqueN(AGS8)], ph[, sum(emk_absorbing)],
  100 * ph[, mean(emk_absorbing)]
))

# -- Standardise continuous regressors on the full estimation sample -----------
# eco_index_L1 already standardised via PCA; no rescaling needed.

z <- function(x) as.numeric(scale(x))

ph[, log_dens_z := z(log_pop_dens)]
ph[, sk_z       := z(log_steuerkraft_L1)]
ph[, bev_z      := z(log1p(bev_stock_p100k_L1))]
ph[, chg_z      := z(log1p(ev_chargepoints_p100k_L1))]


# -- Stadtstaaten exclusion for personnel specs --------------------------------
# Hamburg (02), Bremen (04), Berlin (11): n_vze_personal conflates
# municipal and Länder roles; personnel specs use the no-Stadtstaaten sample.

STADTSTAATEN <- c("02", "04", "11")
ph_ns <- ph[!(AGS2 %in% STADTSTAATEN)]
ph_ns[, sk_z       := z(log_steuerkraft_L1)]   # re-z-score on restricted sample
ph_ns[, log_dens_z := z(log_pop_dens)]
ph_ns[, bev_z      := z(log1p(bev_stock_p100k_L1))]
ph_ns[, chg_z      := z(log1p(ev_chargepoints_p100k_L1))]
ph_ns[, pers_z     := z(log1p(n_vze_personal_L1))]

cat(sprintf(
  "No-Stadtstaaten: %d obs dropped | %d obs | %d AGS8 | %d onsets\n",
  nrow(ph) - nrow(ph_ns), nrow(ph_ns),
  ph_ns[, uniqueN(AGS8)], ph_ns[, sum(emk_absorbing)]
))

# -- Formulas ------------------------------------------------------------------
# Two ecosystem variants:
#   eco  — PCA composite index (eco_index_L1)
#   comp — two separate indicators (bev_z, chg_z)
# Personnel variant: no-Stadtstaaten sample only (see above)
# All specs: year + AGS2 FE; clustering at AGS5.

f_eco_a2 <- emk_absorbing ~
  log_dens_z + state_gruene_L1 +
  eco_index_L1 + sk_z | year + AGS2

f_comp_a2 <- emk_absorbing ~
  log_dens_z + state_gruene_L1 +
  bev_z + chg_z + sk_z | year + AGS2

f_eco_pers <- emk_absorbing ~
  log_dens_z + state_gruene_L1 +
  eco_index_L1 + sk_z + pers_z | year + AGS2

f_comp_pers <- emk_absorbing ~
  log_dens_z + state_gruene_L1 +
  bev_z + chg_z + sk_z + pers_z | year + AGS2

# -- Estimation ----------------------------------------------------------------

cll <- binomial("cloglog")   # primary: discrete-time proportional hazards
lgt <- binomial("logit")     # robustness

cat("\nEstimating cloglog models...\n")
m_cll_eco_a2    <- feglm(f_eco_a2,    data = ph,    family = cll, cluster = ~AGS5)
m_cll_comp_a2   <- feglm(f_comp_a2,   data = ph,    family = cll, cluster = ~AGS5)
m_cll_eco_pers  <- feglm(f_eco_pers,  data = ph_ns, family = cll, cluster = ~AGS5)
m_cll_comp_pers <- feglm(f_comp_pers, data = ph_ns, family = cll, cluster = ~AGS5)

cat("Estimating logit models...\n")
m_lgt_eco_a2    <- feglm(f_eco_a2,    data = ph,    family = lgt, cluster = ~AGS5)
m_lgt_comp_a2   <- feglm(f_comp_a2,   data = ph,    family = lgt, cluster = ~AGS5)
m_lgt_eco_pers  <- feglm(f_eco_pers,  data = ph_ns, family = lgt, cluster = ~AGS5)
m_lgt_comp_pers <- feglm(f_comp_pers, data = ph_ns, family = lgt, cluster = ~AGS5)

# -- Observations per model ----------------------------------------------------

cat("\nObservations used per model:\n")
for (s in list(
  list(m_cll_eco_a2,    "cloglog  eco        AGS2   full sample"),
  list(m_cll_comp_a2,   "cloglog  comp       AGS2   full sample"),
  list(m_cll_eco_pers,  "cloglog  eco+pers   AGS2   no-Stadtstaaten"),
  list(m_cll_comp_pers, "cloglog  comp+pers  AGS2   no-Stadtstaaten"),
  list(m_lgt_eco_a2,    "logit    eco        AGS2   full sample"),
  list(m_lgt_comp_a2,   "logit    comp       AGS2   full sample"),
  list(m_lgt_eco_pers,  "logit    eco+pers   AGS2   no-Stadtstaaten"),
  list(m_lgt_comp_pers, "logit    comp+pers  AGS2   no-Stadtstaaten")
)) cat(sprintf("  %-48s n = %d\n", s[[2]], nobs(s[[1]])))

# -- Output tables -------------------------------------------------------------

dict <- c(
  log_dens_z      = "Log pop. density (z)",
  state_gruene_L1 = "State Gruene vote share (L1)",
  eco_index_L1    = "EV ecosystem index PCA (L1)",
  sk_z            = "log1p Steuerkraft (z, L1)",
  bev_z           = "log(1+BEV stock p100k) (z, L1)",
  chg_z           = "log(1+Charging pts p100k) (z, L1)",
  pers_z          = "log(1+Municipal personnel VZE p100k) (z, L1)"
)

tab_note <- paste0(
  "Discrete-time hazard on AGS8 risk set (year >= 2015). ",
  "Cols (3)-(4): Hamburg, Bremen, Berlin excluded (personnel conflates municipal/Länder roles). ",
  "year + AGS2 FE throughout. year FE = non-parametric baseline hazard. ",
  "SEs clustered at AGS5."
)

col_labels <- c(
  "(1) Eco", "(2) Comp",
  "(3) Eco+P", "(4) Comp+P"
)

cll_models <- list(
  m_cll_eco_a2, m_cll_comp_a2,
  m_cll_eco_pers, m_cll_comp_pers
)
lgt_models <- list(
  m_lgt_eco_a2, m_lgt_comp_a2,
  m_lgt_eco_pers, m_lgt_comp_pers
)

cat("\n=== Cloglog: discrete-time hazard (primary) ===\n")
do.call(etable, c(
  setNames(cll_models, col_labels),
  list(depvar = FALSE, digits = 4, dict = dict, notes = tab_note)
))

cat("\n=== Logit: robustness ===\n")
do.call(etable, c(
  setNames(lgt_models, col_labels),
  list(depvar = FALSE, digits = 4, dict = dict, notes = tab_note)
))

do.call(etable, c(
  setNames(cll_models, col_labels),
  list(
    depvar  = FALSE, digits = 4, dict = dict, notes = tab_note,
    title   = "Discrete-time hazard (cloglog): drivers of EMK treatment onset",
    file    = file.path(out_dir, "tab_onset_cloglog.tex"),
    replace = TRUE
  )
))

do.call(etable, c(
  setNames(lgt_models, col_labels),
  list(
    depvar  = FALSE, digits = 4, dict = dict, notes = tab_note,
    title   = "Logit robustness: drivers of EMK treatment onset",
    file    = file.path(out_dir, "tab_onset_logit.tex"),
    replace = TRUE
  )
))

cat(sprintf("\nTeX tables written to: %s\n", out_dir))
