# ─────────────────────────────────────────────────────────────────────────────
# 04_did_main.R   — Staggered DiD estimates (Callaway-Sant'Anna)
#
# Estimator : did::att_gt; conditional doubly-robust (est_method = "dr") is the
#             PRIMARY spec; unconditional (xformla = NULL) is a comparison column.
# Bootstrap : multiplier, B = 2000, clustered on AGS5 (county)
# ES bands  : simultaneous (cband = TRUE in aggte), pointwise fallback
# Outcome   : bev_neuzulassungen_p100k (level; winsorised at 99th pct)
# Covariates (conditional spec, XFORMLA_CS from _dict.R):
#   sk_base_z + state_green_base_z + dens_base_z
#   (Tax capacity | State Green share | Log pop. density; baseline 2014-16)
#
# Conditional CS is the headline (selection-on-observables / conditional
# parallel trends); the unconditional column is retained only to show
# sensitivity to covariate adjustment. The shared estimation + output machinery
# lives in `_did_helpers.R` and is reused by 07_spillovers.R.
#
# Cohorts with < COHORT_MIN treated units are dropped upstream in
# 00_prep_analysis.R; the frames loaded here are already post-drop.
#
# Output sections (each: <stem>.{png,tex,csv}; conditional primary + uncond)
# ────────────────
# A+B. Combined direct + broad — tab_es_main_combined.*, es_main_combined.png
# C.   Robustness: not-yet-treated (direct, ever-treated pool only)
#                  — es_robust_notyet.*
# D.   ATT summary — est_att_main.csv  (one row per estimated cell)
#
# Secondary outcome ES (corporate/private BEV, ICE placebo) -> 06_heterogeneity.R
#
# All outputs -> 03_output/04_did_main/
# ─────────────────────────────────────────────────────────────────────────────

library(data.table)
library(did)
library(ggplot2)

# ── Paths ─────────────────────────────────────────────────────────────────────

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
  stop("Cannot determine script path. Run as: Rscript 04_did_main.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "03_output", "04_did_main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))
source(file.path(code_dir, "03_analysis", "_did_helpers.R"))
set.seed(202L)

# ── Load frames ───────────────────────────────────────────────────────────────

frames <- list(
  direct = .read_frame(file.path(data_final, "frame_did_direct.csv")),
  broad  = .read_frame(file.path(data_final, "frame_did_broad.csv"))
)

for (nm in names(frames)) {
  d <- frames[[nm]]
  stopifnot(is.integer(d$ags8_id), is.integer(d$gname_cs))
  if (!"state_green_base_z" %in% names(d))
    stop(nm, " frame missing state_green_base_z — re-run 00_prep_analysis.R")
}

# ── A1 guard ──────────────────────────────────────────────────────────────────
# Validates the ESTIMATION sample: the frames here are post-COHORT_MIN-drop, so
# this confirms the never-treated control pool survived the upstream filtering.

for (nm in names(frames)) {
  d            <- frames[[nm]]
  never_u      <- uniqueN(d$ags8_id[d$gname_cs == 0L])
  untreated_n  <- sum(d$gname_cs == 0L)
  cat(sprintf("CS [%s]: never_units=%d; untreated_obs=%d (%.0f%%)\n",
              nm, never_u, untreated_n,
              100 * untreated_n / nrow(d)))
  if (!(untreated_n > 0.5 * nrow(d) && never_u > 1000L))
    stop(sprintf("A1 guard FAILED on frame '%s' (untreated_obs=%d of %d, ",
                 nm, untreated_n, nrow(d)),
         sprintf("never_units=%d). Re-check 00_prep_analysis.R.", never_u))
}

# ── Table note strings (03-specific) ───────────────────────────────────────────

.NOTE_TAIL <- paste0(
  "Outcome: BEV new registrations per 100,000 population ",
  "(level; winsorised at 99th percentile). ",
  "Bootstrap SEs clustered at county (AGS5), $B=", CS_BITERS, "$. ",
  "Graph CI bands: simultaneous (95\\%) where available; pointwise (95\\%) ",
  "fallback when the simultaneous critical value is unavailable. ",
  "Table: pointwise SEs; ",
  "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$ (two-sided). ",
  "Conditional spec covariates (z-scored, earliest year 2014--16): ",
  "tax capacity, state Green vote share, log pop.\\ density. ",
  "Pre-treatment test: the package's joint Wald $\\chi^2$ of $H_0$ that all ",
  "pre-treatment $\\text{ATT}(g,t)=0$ (AGS5-clustered multiplier bootstrap; ",
  "$df=$ number of pre-treatment $(g,t)$ cells), shown as ``--'' if unavailable."
)

NOTE_MAIN <- paste0(
  "Callaway--Sant'Anna (2021) CS estimator. ",
  "The headline column is the conditional spec; ",
  "the unconditional column is shown for comparison. ", .NOTE_TAIL
)

NOTE_NOTYET <- paste0(
  "\\textcite{callaway2021difference} doubly robust estimates for direct treatment ",
  "using not-yet-treated as control group. The dependent variable is BEV registration ",
  "per 100k inhabitants. Conditional specifications control for baseline tax capacity, ",
  "population density, and state Green vote share. SEs are AGS5-clustered derived from ",
  "a multiplier boostrap ($B=", CS_BITERS, "$). $p$-value of a joint pre-test for ",
  "$H_0:\\text{ATT}(g,t)=0$ for all $g,t$ in the last line. ",
  "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$."
)

# ═══════════════════════════════════════════════════════════════════════════════
# A + B  Main specs — direct and broad frames (conditional dr primary | uncond)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n=== A. Direct frame ===\n")
r_dir_dr  <- run_spec(frames$direct, OUTCOME_BEV, XFORMLA_CS, "nevertreated",
                      "dr",  "direct/cond_dr")
r_dir_unc <- run_spec(frames$direct, OUTCOME_BEV, NULL,       "nevertreated",
                      "reg", "direct/uncond")

cat("\n=== B. Broad frame ===\n")
r_brd_dr  <- run_spec(frames$broad, OUTCOME_BEV, XFORMLA_CS, "nevertreated",
                      "dr",  "broad/cond_dr")
r_brd_unc <- run_spec(frames$broad, OUTCOME_BEV, NULL,       "nevertreated",
                      "reg", "broad/uncond")

# ═══════════════════════════════════════════════════════════════════════════════
# C  Robustness — not-yet-treated (direct frame, ever-treated pool only)
# ═══════════════════════════════════════════════════════════════════════════════

# Filter to ever-treated only: never-treated units (gname_cs == 0) are not valid
# controls under notyettreated and would simply be ignored by att_gt, so we drop
# them upfront to make the comparison pool explicit.
dat_notyet <- frames$direct[gname_cs > 0L]

cat("\n=== C. Robustness: not-yet-treated (direct frame) ===\n")
r_nyt_dr  <- run_spec(dat_notyet, OUTCOME_BEV, XFORMLA_CS, "notyettreated",
                      "dr",  "notyet/cond_dr")
r_nyt_unc <- run_spec(dat_notyet, OUTCOME_BEV, NULL,       "notyettreated",
                      "reg", "notyet/uncond")

# ═══════════════════════════════════════════════════════════════════════════════
# Write outputs
# ═══════════════════════════════════════════════════════════════════════════════

# Combined A+B: four-column table
# (direct cond | direct uncond | broad cond | broad uncond)
cat("Writing A+B combined table\n")
.ab <- function(r_list, field) lapply(r_list, `[[`, field)
r_ab <- list(dir_cond_dr = r_dir_dr, dir_uncond = r_dir_unc,
             brd_cond_dr = r_brd_dr, brd_uncond = r_brd_unc)
write_es_longtblr(
  es_list     = .ab(r_ab, "es_agg"),
  att_list    = .ab(r_ab, "att_agg"),
  pre_tests   = .ab(r_ab, "pre_tst"),
  n_treated   = .ab(r_ab, "n_treated"),
  n_control   = .ab(r_ab, "n_control"),
  est_methods = c(dir_cond_dr = "dr", dir_uncond = "reg",
                  brd_cond_dr = "dr", brd_uncond = "reg"),
  resize      = TRUE,
  stem        = file.path(out_dir, "tab_es_main_combined"),
  caption     = paste0("DiD Coefficients: BEV New Registrations per 100k",
                       " --- Direct and Broad Treatment"),
  label       = "tab:es_main_ab_combined",
  note        = paste0(
    "\\textcite{callaway2021difference} doubly robust estimates for direct and broad ",
    "treatment with simultaneous 95\\% confidence bands (multiplier bootstrap, $B=",
    CS_BITERS, "$, clustered at AGS5). Conditional specifications control for baseline ",
    "tax capacity, population density, and state Green vote share. Direct treatment ",
    "estimation has one post-period less due to few post treatment observations. ",
    "$p$-value of a joint pre-test for $H_0:\\text{ATT}(g,t)=0$ for all $g,t$ in the ",
    "last line. $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$."
  )
)

cat("Writing A+B combined graph (2x2: direct top, broad bottom)\n")
# Wide A4-friendly landscape rectangle (2 cols x 2 rows): override the square
# per-panel default so the figure reads horizontally across the page.
es_graph(
  es_list = list(dir_cond_dr = r_dir_dr$es_agg, dir_uncond  = r_dir_unc$es_agg,
                 brd_cond_dr = r_brd_dr$es_agg, brd_uncond  = r_brd_unc$es_agg),
  file    = file.path(out_dir, "es_main_combined.png"),
  ncol    = 2L,
  width   = 11.5,
  height  = 6.0
)
write_fig_tex(
  img_file = file.path(out_dir, "es_main_combined.png"),
  caption  = paste0("Event Study: BEV New Registrations per 100k --- ",
                    "Direct (top) and Broad (bottom) Treatment"),
  label    = "fig:es_main_combined",
  note     = paste0(
    "\\textcite{callaway2021difference} estimates for direct (top) and broad (bottom) ",
    "treatment with simultaneous 95\\% confidence bands (multiplier bootstrap, $B=",
    CS_BITERS, "$, clustered at AGS5). Baseline controls for tax capacity and population ",
    "density included. The dashed vertical line marks the last pre-treatment period. ",
    "Direct treatment estimation has one post-period less due to few post treatment observations."
  )
)

cat("Writing C: robustness not-yet-treated\n")
emit_section(
  results = list(nyt_dr = r_nyt_dr, nyt_uncond = r_nyt_unc),
  ems     = c(nyt_dr = "dr", nyt_uncond = "reg"),
  graph_file  = file.path(out_dir, "es_robust_notyet.png"),
  graph_title = "Robustness: Not-yet-treated Control (Direct Frame)",
  main_stem    = file.path(out_dir, "es_robust_notyet"),
  main_caption = paste0("Robustness: Not-yet-treated Control Group",
                        " (Direct Frame)"),
  main_label   = "tab:es_robust_notyet",
  main_note    = NOTE_NOTYET,
  main_display = c("nyt_dr", "nyt_uncond")
)

# ── D. ATT summary CSV (one row per estimated cell) ────────────────────────────

att_summary <- rbindlist(list(
  .att_row("dir_cond_dr",  "direct", "XFORMLA_CS", "dr",  "nevertreated",
           OUTCOME_BEV,  r_dir_dr),
  .att_row("dir_uncond",   "direct", "none",       "reg", "nevertreated",
           OUTCOME_BEV,  r_dir_unc),
  .att_row("brd_cond_dr",  "broad",  "XFORMLA_CS", "dr",  "nevertreated",
           OUTCOME_BEV,  r_brd_dr),
  .att_row("brd_uncond",   "broad",  "none",       "reg", "nevertreated",
           OUTCOME_BEV,  r_brd_unc),
  .att_row("nyt_cond_dr",  "direct", "XFORMLA_CS", "dr",  "notyettreated",
           OUTCOME_BEV,  r_nyt_dr),
  .att_row("nyt_uncond",   "direct", "none",       "reg", "notyettreated",
           OUTCOME_BEV,  r_nyt_unc)
), use.names = TRUE, fill = TRUE)

fwrite(att_summary, file.path(out_dir, "est_att_main.csv"))
cat(sprintf("\nAll DiD main outputs -> %s\n", out_dir))
