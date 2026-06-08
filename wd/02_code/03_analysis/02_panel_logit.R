library(data.table)
library(fixest)

# -- Paths ---------------------------------------------------------------------

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

# Filtered panel for hazard logit: retain only pre-treatment years and the year
# of first treatment, dropping the absorbing post-treatment spell.
panel[,
  first_treat := if (any(emk_absorbing == 1L))
    min(year[emk_absorbing == 1L])
  else
    .Machine$integer.max,
  by = AGS8
]
panel_hazard <- panel[year <= first_treat][, first_treat := NULL]
panel[, first_treat := NULL]

cat(sprintf(
  "Full panel:   %d obs | %d AGS8\n",
  nrow(panel), panel[, uniqueN(AGS8)]
))
cat(sprintf(
  "Hazard panel: %d obs | %d AGS8 | %d treatment onsets\n",
  nrow(panel_hazard), panel_hazard[, uniqueN(AGS8)],
  panel_hazard[, sum(emk_absorbing)]
))

# -- Specifications ------------------------------------------------------------

# (1) Baseline
fml1 <- emk_absorbing ~ log_pop_dens + q_pendlersaldo +
  muni_gruene_L1 + eco_index_L1 + q_gest_bev_L1 + steuerkraft_sq_L1 |
  year + AGS2

# (2) Baseline + lagged personnel
fml2 <- emk_absorbing ~ log_pop_dens + q_pendlersaldo +
  muni_gruene_L1 + eco_index_L1 + q_gest_bev_L1 + steuerkraft_sq_L1 +
  n_vze_personal_L1 |
  year + AGS2

# (3) EV ecosystem components instead of index
fml3 <- emk_absorbing ~ log_pop_dens + q_pendlersaldo +
  muni_gruene_L1 + bev_stock_p100k_L1 + ev_chargepoints_p100k_L1 +
  q_gest_bev_L1 + steuerkraft_sq_L1 |
  year + AGS2

# -- Estimation ----------------------------------------------------------------

m_lpm_1    <- feols(fml1, data = panel,        cluster = ~AGS2)
m_lpm_2    <- feols(fml2, data = panel,        cluster = ~AGS2)
m_lpm_3    <- feols(fml3, data = panel,        cluster = ~AGS2)

m_logit_1  <- feglm(fml1, data = panel,        family = binomial("logit"), cluster = ~AGS2)
m_logit_2  <- feglm(fml2, data = panel,        family = binomial("logit"), cluster = ~AGS2)
m_logit_3  <- feglm(fml3, data = panel,        family = binomial("logit"), cluster = ~AGS2)

m_hazard_1 <- feglm(fml1, data = panel_hazard, family = binomial("logit"), cluster = ~AGS2)
m_hazard_2 <- feglm(fml2, data = panel_hazard, family = binomial("logit"), cluster = ~AGS2)
m_hazard_3 <- feglm(fml3, data = panel_hazard, family = binomial("logit"), cluster = ~AGS2)

# -- Output --------------------------------------------------------------------

etable(
  m_lpm_1,    m_lpm_2,    m_lpm_3,
  m_logit_1,  m_logit_2,  m_logit_3,
  m_hazard_1, m_hazard_2, m_hazard_3,
  headers = c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"),
  depvar  = FALSE,
  digits  = 4
)

etable(
  m_lpm_1,    m_lpm_2,    m_lpm_3,
  m_logit_1,  m_logit_2,  m_logit_3,
  m_hazard_1, m_hazard_2, m_hazard_3,
  headers = c("(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)"),
  depvar  = FALSE,
  digits  = 4,
  file    = file.path(out_dir, "tab_treatment_onset.tex"),
  replace = TRUE
)

cat(sprintf(
  "\nTeX table written to: %s\n",
  file.path(out_dir, "tab_treatment_onset.tex")
))
