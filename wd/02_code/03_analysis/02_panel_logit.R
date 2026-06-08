library(data.table)
library(fixest)

# ── Paths ──────────────────────────────────────────────────────────────────────

argv      <- commandArgs(trailingOnly = FALSE)
self_flag <- grep("--file=", argv, value = TRUE)
if (!length(self_flag)) stop("Run as: Rscript 02_panel_logit.R")
root       <- dirname(dirname(dirname(normalizePath(sub("--file=", "", self_flag)))))
DATA_FINAL <- file.path(root, "01_data", "03_final")
OUT_DIR    <- file.path(root, "04_results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load and prepare ───────────────────────────────────────────────────────────

panel <- fread(
  file.path(DATA_FINAL, "emk_inkar_panel_ags8.csv"),
  colClasses = list(character = c("AGS8", "AGS5", "AGS2"))
)

panel <- panel[year > 2015]

# Filtered panel for hazard logit: retain only pre-treatment years and the year
# of first treatment, dropping the absorbing post-treatment spell.
panel[, first_treat := fifelse(
  any(emk_absorbing == 1L),
  min(year[emk_absorbing == 1L]),
  .Machine$integer.max
), by = AGS8]
panel_hazard <- panel[year <= first_treat][, first_treat := NULL]
panel[, first_treat := NULL]

cat(sprintf("Full panel:   %d obs | %d AGS8\n", nrow(panel), panel[, uniqueN(AGS8)]))
cat(sprintf("Hazard panel: %d obs | %d AGS8 | %d treatment onsets\n",
  nrow(panel_hazard), panel_hazard[, uniqueN(AGS8)], panel_hazard[, sum(emk_absorbing)]))

# ── Specification ──────────────────────────────────────────────────────────────

fml <- emk_absorbing ~ log_pop_dens + q_pendlersaldo +
  muni_gruene_L1 + eco_index_L1 + q_gest_bev_L1 + steuerkraft_sq_L1 |
  year + AGS2

# ── Estimation ─────────────────────────────────────────────────────────────────

m_lpm    <- feols(fml, data = panel,        cluster = ~AGS2)
m_logit  <- feglm(fml, data = panel,        family = binomial("logit"), cluster = ~AGS2)
m_hazard <- feglm(fml, data = panel_hazard, family = binomial("logit"), cluster = ~AGS2)

etable(
  m_lpm, m_logit, m_hazard,
  headers = c("LPM", "Logit", "Hazard Logit"),
  depvar  = FALSE,
  digits  = 4
)

etable(
  m_lpm, m_logit, m_hazard,
  headers = c("LPM", "Logit", "Hazard Logit"),
  depvar  = FALSE,
  digits  = 4,
  file    = file.path(OUT_DIR, "tab_treatment_onset.tex"),
  replace = TRUE
)

cat(sprintf("\nTeX table written to: %s\n", file.path(OUT_DIR, "tab_treatment_onset.tex")))
