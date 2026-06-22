# ───────────────────────────────────────────────────────────────────────────────
# 00_prep_analysis.R
#
# Build estimation frames (.csv) consumed by every downstream analysis script.
#
# Outputs (01_data/03_final/):
#   frame_hazard.csv         — risk set on DIRECT onsets; broad units stay in
#                              risk set with kreis_funded switching on (year >= 2015).
#   frame_hazard_cov.csv     — risk set on BROAD (coverage) onsets (appendix).
#   frame_logit_full.csv     — full panel (no censoring), DIRECT onset flag;
#                              counterfactual to frame_hazard for plain logit.
#   frame_logit_cov_full.csv — full panel, BROAD onset flag.
#   frame_did_broad.csv      — full panel for staggered DiD on broad treatment.
#   frame_did_direct.csv     — broadcast-only units dropped entirely.
#
# Also writes a baseline-snapshot table (earliest available 2014–2016 value per
# AGS8 per covariate) used as a time-invariant covariate set for CS-dr and for
# tercile stratification in heterogeneity.
#
# Outputs (04_results/00_prep_analysis/):
#   tab_cohorts.{tex,csv}   — onsets per year × {broad, direct}, pre/post counts
#   prep_log.txt            — provenance log: n, dropped, base-year choices
# ───────────────────────────────────────────────────────────────────────────────

library(data.table)
library(fixest)        # etable

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
  stop("Cannot determine script path. Run as: Rscript 00_prep_analysis.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "04_results", "00_prep_analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))

# -- Constants -----------------------------------------------------------------
# Latest data year for first pass (plan: toggle via MAX_COHORT once policy-regime
# scope is decided). NA disables the cohort cap.
MAX_COHORT  <- NA_integer_
BASE_WINDOW <- 2014:2016
START_YEAR  <- 2015L
MIN_YEAR    <- 2010L     # drop pre-MIN_YEAR rows; trim panel before save

# -- Load panel + neighbors ----------------------------------------------------

panel <- fread(
  file.path(data_final, "emk_inkar_panel_ags8.csv"),
  colClasses = list(character = c("AGS8", "AGS5", "AGS2"))
)
nbrs <- fread(
  file.path(data_final, "spatial_neighbors_ags8.csv"),
  colClasses = list(character = "AGS8")
)
panel <- merge(panel, nbrs, by = c("AGS8", "year"), all.x = TRUE)

# Defensive zfill
panel[, AGS8 := sprintf("%08d", as.integer(AGS8))]
panel[, AGS5 := sprintf("%05d", as.integer(AGS5))]
panel[, AGS2 := sprintf("%02d", as.integer(AGS2))]

# Stadtstaaten flag (kept as column, not a separate frame)
panel[, ns_flag := !(AGS2 %in% STADTSTAATEN)]

# Year filter: trim pre-MIN_YEAR rows. Baseline snapshot (BASE_WINDOW = 2014–
# 2016) and hazard z-scoring (year >= START_YEAR = 2015) are both unaffected.
# The DiD frames otherwise hold many years of NA-outcome rows that downstream
# scripts would drop per call; this filter just does it once at the source.
.n_pre <- nrow(panel)
panel <- panel[year >= MIN_YEAR]
cat(sprintf("MIN_YEAR=%d filter: -%d rows; panel years now %d-%d\n",
            MIN_YEAR, .n_pre - nrow(panel),
            min(panel$year), max(panel$year)))

# Early-treatment hygiene: EMK projects can begin before MIN_YEAR (e.g. some
# 2008-09 projects). Those units would otherwise appear "always-treated" in
# the post-MIN_YEAR frame, which CS/BJS handle poorly. NA-out their treatment
# year so they are reclassified as untreated controls. Print counts so it's
# transparent how many units are affected.
.early_direct <- panel[first_treat_direct < MIN_YEAR, uniqueN(AGS8)]
.early_broad  <- panel[first_treat_broad  < MIN_YEAR, uniqueN(AGS8)]
if (.early_direct > 0L || .early_broad > 0L) {
  cat(sprintf("Pre-MIN_YEAR treatment: %d AGS8 with first_treat_direct < %d, ",
              .early_direct, MIN_YEAR))
  cat(sprintf("%d with first_treat_broad < %d. Reclassified as untreated.\n",
              .early_broad, MIN_YEAR))
  panel[first_treat_direct < MIN_YEAR, first_treat_direct := NA_integer_]
  panel[first_treat_broad  < MIN_YEAR, first_treat_broad  := NA_integer_]
  panel[is.na(first_treat_direct) & is.na(first_treat_broad),
        treat_type := "never"]
}

# Data-quality filter: drop very small Gemeinden. KBA registers vehicles at
# the holder's HQ AGS8; leasing companies / Großkunden-Halter incorporated in
# tiny Gemeinden create extreme per-100k spikes (e.g. AGS8 01054108, pop 304,
# 143 BEVs in 2022 -> 47,000 per 100k). UNIT-LEVEL filter: drop AGS8 whose
# minimum population in the panel is below threshold, so the panel doesn't
# acquire churning unit-years (which would silently censor the hazard and
# change which (g,t) cells exist in CS).
POP_MIN <- 500L
.n_pre <- nrow(panel)
.small_units <- panel[, .(min_pop = min(xbev, na.rm = TRUE)), by = AGS8][
  min_pop < POP_MIN, AGS8]
.dropped_treat_broad  <- uniqueN(
  panel[AGS8 %in% .small_units & !is.na(first_treat_broad),  AGS8])
.dropped_treat_direct <- uniqueN(
  panel[AGS8 %in% .small_units & !is.na(first_treat_direct), AGS8])
panel <- panel[!(AGS8 %in% .small_units)]
cat(sprintf(
  "POP_MIN=%d unit-level filter: -%d rows (%d AGS8 dropped; %d broad-treated, %d direct-treated)\n",
  POP_MIN, .n_pre - nrow(panel), length(.small_units),
  .dropped_treat_broad, .dropped_treat_direct))

# Winsorize per-100k outcomes at the 99th percentile to cap remaining
# fleet-registration artifacts at larger Gemeinden.
WINSOR_Q <- 0.99
.winsor_cols <- c("bev_neuzulassungen_p100k", "bev_corporate_p100k",
                  "bev_private_p100k", "bev_stock_p100k",
                  "ice_neuzulassungen_p100k")
for (.c in .winsor_cols) {
  .cap <- quantile(panel[[.c]], WINSOR_Q, na.rm = TRUE)
  .n_capped <- panel[get(.c) > .cap, .N]
  panel[get(.c) > .cap, (.c) := .cap]
  cat(sprintf("Winsorized %-26s at %.1f (q=%.3f, %d cells capped)\n",
              .c, .cap, WINSOR_Q, .n_capped))
}

# log1p twins of the per-100k outcomes for proportional-effect specifications.
# log1p handles the floor of zeros (many small Gemeinden have 0 new BEVs in
# early years); the level reads as "additional registrations per 100k", the
# log1p reads as approximate percent change in the rate.
for (.c in .winsor_cols) {
  panel[, (paste0("log1p_", .c)) := log1p(pmax(get(.c), 0))]
}

# Optional cohort cap: censor units whose first-treat exceeds MAX_COHORT.
if (!is.na(MAX_COHORT)) {
  panel[first_treat_direct > MAX_COHORT, first_treat_direct := NA_integer_]
  panel[first_treat_broad  > MAX_COHORT, first_treat_broad  := NA_integer_]
}

cat(sprintf("Panel: %d obs | %d AGS8 | years %d–%d\n",
            nrow(panel), uniqueN(panel$AGS8),
            min(panel$year), max(panel$year)))
cat(sprintf("  ever-broad: %d | ever-direct: %d | never: %d\n",
            uniqueN(panel[!is.na(first_treat_broad),  AGS8]),
            uniqueN(panel[!is.na(first_treat_direct), AGS8]),
            uniqueN(panel[treat_type == "never",      AGS8])))

# -- Baseline snapshot per AGS8 ------------------------------------------------
# Per-variable earliest available year in BASE_WINDOW (Kaufkraft coverage may
# start late, so we take the per-variable earliest year). Records which year
# was used per variable in the provenance log.

base_vars <- c(
  kk_base          = "log_kaufkraft",
  sk_base          = "log_steuerkraft",
  dens_base        = "log_pop_dens",
  bev_base         = "bev_stock_p100k",
  chg_base         = "ev_chargepoints_p100k",
  green_base       = "muni_gruene",
  state_green_base = "state_gruene",
  pop_base         = "xbev"
)

snap <- panel[year %in% BASE_WINDOW, c("AGS8", "year", unname(base_vars)),
              with = FALSE]
base_dt <- snap[, .(AGS8 = unique(AGS8))]
year_log <- list()
for (nm in names(base_vars)) {
  v  <- base_vars[[nm]]
  d  <- snap[!is.na(get(v)), .(AGS8, year, val = get(v))]
  setorder(d, AGS8, year)
  d  <- d[, .SD[1L], by = AGS8]
  yr_col <- paste0(nm, "_year")
  setnames(d, c("year", "val"), c(yr_col, nm))
  base_dt <- merge(base_dt, d, by = "AGS8", all.x = TRUE)
  year_log[[nm]] <- d[, .N, by = c(yr_col)][order(get(yr_col))]
}

# log1p transforms for the BEV / charging baselines (snapshots, not yet z'd)
base_dt[, bev_base   := log1p(pmax(bev_base,   0))]
base_dt[, chg_base   := log1p(pmax(chg_base,   0))]

# Population-weighted and unweighted tercile / quintile cuts on baseline KK & SK
.qcut <- function(v, w = NULL, n = 5L) {
  if (is.null(w)) {
    cut(v,
        breaks = quantile(v, probs = seq(0, 1, length.out = n + 1L), na.rm = TRUE),
        include.lowest = TRUE, labels = FALSE)
  } else {
    cuts <- wq(v, w, probs = seq(0, 1, length.out = n + 1L)[-c(1L, n + 1L)])
    cuts <- c(-Inf, cuts, Inf)
    cut(v, breaks = cuts, include.lowest = TRUE, labels = FALSE)
  }
}

for (nm in c("kk_base", "sk_base")) {
  base_dt[, paste0(nm, "_q5")      := .qcut(get(nm), n = 5L)]
  base_dt[, paste0(nm, "_terc")    := .qcut(get(nm), n = 3L)]
  base_dt[, paste0(nm, "_q5_pw")   := .qcut(get(nm), w = pop_base, n = 5L)]
  base_dt[, paste0(nm, "_terc_pw") := .qcut(get(nm), w = pop_base, n = 3L)]
}

# z-scored baselines (on the AGS8 cross-section); kept as time-invariant covs
for (nm in c("kk_base", "sk_base", "dens_base",
             "bev_base", "chg_base",
             "green_base", "state_green_base")) {
  base_dt[, paste0(nm, "_z") := z(get(nm))]
}

cat(sprintf("Baseline snapshot: %d AGS8 with at least one base var\n", nrow(base_dt)))

# B5: NA diagnostics on the baseline z-covariates used by CS-dr. If any of
# these are systematically missing in a Bundesland-shaped pattern, the CS
# sample silently loses the units and the headline ATT is on a non-random
# subset. Print counts; if green_base is the binding loser, swap in a fallback
# from state_gruene baseline.
cat("\nBaseline covariate NA counts (relevant to CS sample):\n")
.diag_cols <- c("kk_base_z", "sk_base_z", "dens_base_z",
                "green_base_z", "state_green_base_z",
                "bev_base_z", "chg_base_z")
for (.c in .diag_cols) {
  cat(sprintf("  %-20s NA=%d (%.1f%% of AGS8)\n",
              .c, sum(is.na(base_dt[[.c]])),
              100 * mean(is.na(base_dt[[.c]]))))
}
.cov_keep <- complete.cases(
  base_dt[, c("sk_base_z", "state_green_base_z", "dens_base_z"),
          with = FALSE])
cat(sprintf(
  "CS xformla complete cases (sk+state_green+dens): %d / %d AGS8\n",
  sum(.cov_keep), nrow(base_dt)))

# green_base fallback: per AGS8, where muni_gruene is missing for all base
# years but state_gruene is present, use the state-level baseline as a fallback.
# Document the swap so it's traceable.
if (any(is.na(base_dt$green_base))) {
  state_snap <- panel[year %in% BASE_WINDOW & !is.na(state_gruene),
                      .(state_gruene_base = state_gruene[1L]),
                      by = AGS8]
  base_dt <- merge(base_dt, state_snap, by = "AGS8", all.x = TRUE)
  .n_swap <- base_dt[is.na(green_base) & !is.na(state_gruene_base), .N]
  base_dt[is.na(green_base) & !is.na(state_gruene_base),
          green_base := state_gruene_base]
  base_dt[, green_base_z := z(green_base)]
  base_dt[, state_gruene_base := NULL]
  cat(sprintf("green_base: filled %d AGS8 from state_gruene baseline; ",
              .n_swap))
  cat(sprintf("remaining NA = %d\n", sum(is.na(base_dt$green_base))))
}

fwrite(base_dt, file.path(out_dir, "baseline_snapshot.csv"))

# Merge baseline into panel (constant within AGS8)
panel <- merge(panel, base_dt, by = "AGS8", all.x = TRUE)

# -- Hazard frame (DIRECT events, AGS8 risk set) ------------------------------
# Risk set: year >= 2015; drop unit-years AFTER the unit's direct onset.
# Broadcast-only units (treat_type == "broadcast_only") remain in the risk set
# (they never have a direct event); kreis_funded switches on past the AGS5
# Kreis project's start_year.

ph <- panel[year >= START_YEAR]
ph <- ph[is.na(first_treat_direct) | year <= first_treat_direct]
ph[, onset_direct := as.integer(!is.na(first_treat_direct) & year == first_treat_direct)]
ph[, kreis_funded := as.integer(!is.na(kreis_funded_year) & kreis_funded_year < year)]

cat(sprintf(
  "frame_hazard: %d obs | %d AGS8 | %d direct onsets | rate %.2f%%\n",
  nrow(ph), uniqueN(ph$AGS8), ph[, sum(onset_direct)],
  100 * ph[, mean(onset_direct)]
))
fwrite(ph, file.path(data_final, "frame_hazard.csv"))

# Coverage-event hazard frame (appendix; "determinants of coverage")
ph_cov <- panel[year >= START_YEAR]
ph_cov <- ph_cov[is.na(first_treat_broad) | year <= first_treat_broad]
ph_cov[, onset_broad := as.integer(!is.na(first_treat_broad) & year == first_treat_broad)]
ph_cov[, kreis_funded := as.integer(!is.na(kreis_funded_year) & kreis_funded_year < year)]
cat(sprintf(
  "frame_hazard_cov: %d obs | %d AGS8 | %d coverage onsets | rate %.2f%%\n",
  nrow(ph_cov), uniqueN(ph_cov$AGS8), ph_cov[, sum(onset_broad)],
  100 * ph_cov[, mean(onset_broad)]
))
fwrite(ph_cov, file.path(data_final, "frame_hazard_cov.csv"))

# -- Logit-counterfactual frames (full panel, no risk-set censoring) ----------
# Same year >= START_YEAR floor and all other filters; post-onset AGS8-years
# are retained. Onset stays a year-of-event indicator (mechanically zero for
# post-onset years). Consumed by 02b_logit_uncensored.R as a plain-logit
# counterfactual to the discrete-time hazard.

logit_full <- panel[year >= START_YEAR]
logit_full[, onset_direct := as.integer(!is.na(first_treat_direct) & year == first_treat_direct)]
logit_full[, kreis_funded := as.integer(!is.na(kreis_funded_year) & kreis_funded_year < year)]
cat(sprintf(
  "frame_logit_full: %d obs | %d AGS8 | %d direct onsets | rate %.2f%%\n",
  nrow(logit_full), uniqueN(logit_full$AGS8), logit_full[, sum(onset_direct)],
  100 * logit_full[, mean(onset_direct)]
))
fwrite(logit_full, file.path(data_final, "frame_logit_full.csv"))

logit_cov_full <- panel[year >= START_YEAR]
logit_cov_full[, onset_broad  := as.integer(!is.na(first_treat_broad)  & year == first_treat_broad)]
logit_cov_full[, kreis_funded := as.integer(!is.na(kreis_funded_year) & kreis_funded_year < year)]
cat(sprintf(
  "frame_logit_cov_full: %d obs | %d AGS8 | %d coverage onsets | rate %.2f%%\n",
  nrow(logit_cov_full), uniqueN(logit_cov_full$AGS8), logit_cov_full[, sum(onset_broad)],
  100 * logit_cov_full[, mean(onset_broad)]
))
fwrite(logit_cov_full, file.path(data_final, "frame_logit_cov_full.csv"))

# -- DiD frames ----------------------------------------------------------------
# `did` package: gname = 0 for never-treated
# `didimputation` (BJS) package: gname = Inf for never-treated
# Integer ags8_id required by both.

panel[, ags8_id := .GRP, by = AGS8]

frame_did_broad <- copy(panel)
frame_did_broad[, gname_cs  := fifelse(is.na(first_treat_broad), 0L,
                                       as.integer(first_treat_broad))]
# BJS never-treated coding lives in `BJS_NEVER` (_dict.R) so 03/04/06 cannot
# drift apart. The A1 guard at the top of 03_did_main.R verifies that the
# untreated estimation sample is non-degenerate; if it is, flip `BJS_NEVER`.
frame_did_broad[, gname_bjs := fifelse(is.na(first_treat_broad), BJS_NEVER,
                                       as.numeric(first_treat_broad))]
fwrite(frame_did_broad, file.path(data_final, "frame_did_broad.csv"))

frame_did_direct <- panel[treat_type != "broadcast_only"]
frame_did_direct[, gname_cs  := fifelse(is.na(first_treat_direct), 0L,
                                        as.integer(first_treat_direct))]
frame_did_direct[, gname_bjs := fifelse(is.na(first_treat_direct), BJS_NEVER,
                                        as.numeric(first_treat_direct))]
fwrite(frame_did_direct, file.path(data_final, "frame_did_direct.csv"))

cat(sprintf(
  "frame_did_broad : %d obs | %d AGS8 | %d ever-treated\n",
  nrow(frame_did_broad), uniqueN(frame_did_broad$AGS8),
  uniqueN(frame_did_broad[gname_cs > 0, AGS8])
))
cat(sprintf(
  "frame_did_direct: %d obs | %d AGS8 | %d ever-treated (broadcast-only dropped)\n",
  nrow(frame_did_direct), uniqueN(frame_did_direct$AGS8),
  uniqueN(frame_did_direct[gname_cs > 0, AGS8])
))

# -- Cohort table --------------------------------------------------------------
# Onsets per year × {broad, direct}; available pre/post given data window.

yr_max <- panel[, max(year)]
yr_min <- panel[, min(year)]

cohort_tab <- rbindlist(list(
  data.table(
    frame  = "broad",
    cohort = sort(unique(panel[!is.na(first_treat_broad),  first_treat_broad])),
    key    = "cohort"
  ),
  data.table(
    frame  = "direct",
    cohort = sort(unique(panel[!is.na(first_treat_direct), first_treat_direct])),
    key    = "cohort"
  )
), fill = TRUE)

cohort_tab[, n_units := vapply(seq_len(.N), function(i) {
  fr <- frame[i]; ck <- cohort[i]
  if (fr == "broad")  uniqueN(panel[first_treat_broad  == ck, AGS8])
  else                uniqueN(panel[first_treat_direct == ck, AGS8])
}, integer(1))]
cohort_tab[, post_avail := yr_max - cohort]
cohort_tab[, pre_avail  := cohort - yr_min]
setcolorder(cohort_tab, c("frame", "cohort", "n_units", "pre_avail", "post_avail"))

cat("\n=== Cohort table ===\n")
print(cohort_tab)
fwrite(cohort_tab, file.path(out_dir, "tab_cohorts.csv"))

# tex via fixest internal helper (no fixest model -> manual)
tex_path <- file.path(out_dir, "tab_cohorts.tex")
tex_lines <- c(
  "\\begin{tabular}{llrrr}",
  "\\hline",
  "Frame & Cohort & N units & Pre-periods & Post-periods \\\\",
  "\\hline",
  cohort_tab[, sprintf("%s & %d & %d & %d & %d \\\\",
                       frame, cohort, n_units, pre_avail, post_avail)],
  "\\hline",
  "\\end{tabular}"
)
writeLines(tex_lines, tex_path)

zero_post <- cohort_tab[post_avail <= 0L]
if (nrow(zero_post)) {
  cat(sprintf(
    "WARNING: %d cohort(s) have 0 post-periods; `did` will drop them:\n",
    nrow(zero_post)))
  print(zero_post)
}

# -- Provenance log ------------------------------------------------------------

log_path <- file.path(out_dir, "prep_log.txt")
con <- file(log_path, open = "wt")
writeLines(c(
  sprintf("Run: %s", format(Sys.time())),
  sprintf("MAX_COHORT: %s", ifelse(is.na(MAX_COHORT), "NA (no cap)", MAX_COHORT)),
  sprintf("Panel obs: %d | AGS8: %d | years %d–%d",
          nrow(panel), uniqueN(panel$AGS8), yr_min, yr_max),
  "",
  "Baseline-window source year used per variable (counts):"
), con)
for (nm in names(year_log)) {
  writeLines(sprintf("  %s:", nm), con)
  write.table(year_log[[nm]], con, row.names = FALSE, quote = FALSE)
}
close(con)
cat(sprintf("\nLogs / cohort table -> %s\n", out_dir))
cat(sprintf("Frames saved        -> %s/frame_*.csv\n", data_final))
