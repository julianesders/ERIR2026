# ─────────────────────────────────────────────────────────────────────────────
# 03_did_main.R   — Staggered DiD estimates (Callaway-Sant'Anna)
#
# Estimator : did::att_gt; est_method threaded per spec ("reg" / "dr")
# Bootstrap : multiplier, B = 2000, clustered on AGS5 (county)
# ES bands  : simultaneous (cband = TRUE in aggte), pointwise fallback
# Outcome   : bev_neuzulassungen_p100k (level; winsorised at 99th pct)
# Covariates (conditional spec, XFORMLA_CS from _dict.R):
#   sk_base_z + state_green_base_z + dens_base_z
#   (Tax capacity | State Green share | Log pop. density; baseline 2014-16)
#
# The CS estimation + output machinery (run_spec, es_graph, write_es_longtblr,
# emit_section, .att_row, constants, plot theme) lives in `_did_helpers.R` and
# is SHARED with 06_spillovers.R, so every CS section — main or spillover —
# produces byte-identical figures/tables/CSVs; only the SAMPLE differs.
#
# Conditional specs are run with BOTH est_method = "reg" (outcome regression)
# and est_method = "dr" (doubly robust). For the conditional sections the MAIN
# table shows Unconditional + Conditional (reg); the dr twin goes to an appendix
# table (`*_dr.{tex,csv}`). The unconditional column is estimator-invariant
# (reg == dr when xformla = NULL) and is reported once.
#
# Cohorts with < COHORT_MIN treated units are dropped upstream in
# 00_prep_analysis.R; the frames loaded here are already post-drop.
#
# Output sections
# ────────────────
# A. Direct frame  — es_main_direct.{png,tex,csv} + es_main_direct_dr.{tex,csv}
# B. Broad frame   — es_main_broad.{png,tex,csv}  + es_main_broad_dr.{tex,csv}
# C. Robustness: not-yet-treated (direct, ever-treated pool only)
#                  — es_robust_notyet.{png,tex,csv} + es_robust_notyet_dr.{tex,csv}
# D. Corporate BEV — es_corp.{png,tex,csv}
# E. Private BEV   — es_priv.{png,tex,csv}
# F. ICE placebo   — es_ice.{png,tex,csv}
# G. ATT summary   — est_att_main.csv  (one row per estimated cell)
#
# All outputs -> 04_results/03_did_main/
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
  stop("Cannot determine script path. Run as: Rscript 03_did_main.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "04_results", "03_did_main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))
source(file.path(code_dir, "03_analysis", "_did_helpers.R"))

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
# Dropping the ~6 small-cohort direct units leaves the never-treated pool
# unchanged, so this should pass; the guard exists to catch upstream breakage.

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
  "fallback when the simultaneous critical value is unavailable ",
  "(the conditional dr column may degrade to pointwise). ",
  "Table: pointwise SEs; ",
  "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$ (two-sided). ",
  "Conditional spec covariates (z-scored, earliest year 2014--16): ",
  "tax capacity, state Green vote share, log pop.\\ density. ",
  "Pre-treatment test: the package's joint Wald $\\chi^2$ of $H_0$ that all ",
  "pre-treatment $\\text{ATT}(g,t)=0$ (AGS5-clustered multiplier bootstrap; ",
  "$df=$ number of pre-treatment $(g,t)$ cells), shown as ``--'' if unavailable."
)

NOTE_MAIN <- paste0(
  "Callaway--Sant'Anna (2021) CS estimator. Main columns: unconditional and ",
  "conditional outcome-regression (reg); the conditional doubly-robust (dr) ",
  "twin is in the appendix. ", .NOTE_TAIL
)

NOTE_DR <- paste0(
  "Callaway--Sant'Anna (2021) CS estimator, conditional doubly-robust (dr) ",
  "twin to the main event-study table. ", .NOTE_TAIL
)

NOTE_NOTYET <- paste0(
  NOTE_MAIN,
  " Control group: not-yet-treated cohorts (ever-treated units only)."
)

NOTE_NOTYET_DR <- paste0(
  NOTE_DR,
  " Control group: not-yet-treated cohorts (ever-treated units only)."
)

note_sec <- function(outcome_lbl)
  paste0(
    "Specification as per the unconditional direct-frame spec ",
    "(nevertreated control, no covariates), but outcome: ",
    outcome_lbl, "."
  )

# ═══════════════════════════════════════════════════════════════════════════════
# A + B  Main specs — direct and broad frames (uncond | cond reg | cond dr)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n=== A. Direct frame ===\n")
r_dir_unc <- run_spec(frames$direct, OUTCOME_BEV, NULL,       "nevertreated",
                      "reg", "direct/uncond")
r_dir_reg <- run_spec(frames$direct, OUTCOME_BEV, XFORMLA_CS, "nevertreated",
                      "reg", "direct/cond_reg")
r_dir_dr  <- run_spec(frames$direct, OUTCOME_BEV, XFORMLA_CS, "nevertreated",
                      "dr",  "direct/cond_dr")

cat("\n=== B. Broad frame ===\n")
r_brd_unc <- run_spec(frames$broad, OUTCOME_BEV, NULL,       "nevertreated",
                      "reg", "broad/uncond")
r_brd_reg <- run_spec(frames$broad, OUTCOME_BEV, XFORMLA_CS, "nevertreated",
                      "reg", "broad/cond_reg")
r_brd_dr  <- run_spec(frames$broad, OUTCOME_BEV, XFORMLA_CS, "nevertreated",
                      "dr",  "broad/cond_dr")

# ═══════════════════════════════════════════════════════════════════════════════
# C  Robustness — not-yet-treated (direct frame, ever-treated pool only)
# ═══════════════════════════════════════════════════════════════════════════════

# Filter to ever-treated only: never-treated units (gname_cs == 0) are not valid
# controls under notyettreated and would simply be ignored by att_gt, so we drop
# them upfront to make the comparison pool explicit.
dat_notyet <- frames$direct[gname_cs > 0L]

cat("\n=== C. Robustness: not-yet-treated (direct frame) ===\n")
r_nyt_unc <- run_spec(dat_notyet, OUTCOME_BEV, NULL,       "notyettreated",
                      "reg", "notyet/uncond")
r_nyt_reg <- run_spec(dat_notyet, OUTCOME_BEV, XFORMLA_CS, "notyettreated",
                      "reg", "notyet/reg")
r_nyt_dr  <- run_spec(dat_notyet, OUTCOME_BEV, XFORMLA_CS, "notyettreated",
                      "dr",  "notyet/dr")

# ═══════════════════════════════════════════════════════════════════════════════
# D – F  Secondary outcomes (direct, nevertreated, unconditional)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n=== D-F. Secondary outcomes ===\n")
r_corp <- run_spec(frames$direct, OUTCOME_CORP, NULL, "nevertreated",
                   "reg", "direct/corp")
r_priv <- run_spec(frames$direct, OUTCOME_PRIV, NULL, "nevertreated",
                   "reg", "direct/priv")
r_ice  <- run_spec(frames$direct, OUTCOME_ICE,  NULL, "nevertreated",
                   "reg", "direct/ice")

# ═══════════════════════════════════════════════════════════════════════════════
# Write outputs
# ═══════════════════════════════════════════════════════════════════════════════

# ── A. Direct frame ────────────────────────────────────────────────────────────

cat("\nWriting A: direct frame\n")
emit_section(
  results = list(uncond = r_dir_unc, cond_reg = r_dir_reg, cond_dr = r_dir_dr),
  ems     = c(uncond = "reg", cond_reg = "reg", cond_dr = "dr"),
  graph_file  = file.path(out_dir, "es_main_direct.png"),
  graph_title = "BEV New Registrations — Direct Treatment Frame",
  main_stem    = file.path(out_dir, "es_main_direct"),
  main_caption = paste0("Event Study: BEV New Registrations per 100k",
                        " --- Direct Treatment"),
  main_label   = "tab:es_main_direct",
  main_note    = NOTE_MAIN,
  main_display = c("uncond", "cond_reg"),
  dr_key       = "cond_dr",
  dr_stem      = file.path(out_dir, "es_main_direct_dr"),
  dr_caption   = paste0("Event Study (Appendix, doubly-robust): BEV New",
                        " Registrations per 100k --- Direct Treatment"),
  dr_label     = "tab:es_main_direct_dr",
  dr_note      = NOTE_DR
)

# ── B. Broad frame ─────────────────────────────────────────────────────────────

cat("Writing B: broad frame\n")
emit_section(
  results = list(uncond = r_brd_unc, cond_reg = r_brd_reg, cond_dr = r_brd_dr),
  ems     = c(uncond = "reg", cond_reg = "reg", cond_dr = "dr"),
  graph_file  = file.path(out_dir, "es_main_broad.png"),
  graph_title = "BEV New Registrations --- Broad Treatment Frame",
  main_stem    = file.path(out_dir, "es_main_broad"),
  main_caption = paste0("Event Study: BEV New Registrations per 100k",
                        " --- Broad Treatment"),
  main_label   = "tab:es_main_broad",
  main_note    = NOTE_MAIN,
  main_display = c("uncond", "cond_reg"),
  dr_key       = "cond_dr",
  dr_stem      = file.path(out_dir, "es_main_broad_dr"),
  dr_caption   = paste0("Event Study (Appendix, doubly-robust): BEV New",
                        " Registrations per 100k --- Broad Treatment"),
  dr_label     = "tab:es_main_broad_dr",
  dr_note      = NOTE_DR
)

# ── C. Robustness: not-yet-treated ─────────────────────────────────────────────

cat("Writing C: robustness not-yet-treated\n")
emit_section(
  results = list(nyt_uncond = r_nyt_unc, nyt_reg = r_nyt_reg,
                 nyt_dr = r_nyt_dr),
  ems     = c(nyt_uncond = "reg", nyt_reg = "reg", nyt_dr = "dr"),
  graph_file  = file.path(out_dir, "es_robust_notyet.png"),
  graph_title = "Robustness: Not-yet-treated Control (Direct Frame)",
  main_stem    = file.path(out_dir, "es_robust_notyet"),
  main_caption = paste0("Robustness: Not-yet-treated Control Group",
                        " (Direct Frame)"),
  main_label   = "tab:es_robust_notyet",
  main_note    = NOTE_NOTYET,
  main_display = c("nyt_uncond", "nyt_reg"),
  dr_key       = "nyt_dr",
  dr_stem      = file.path(out_dir, "es_robust_notyet_dr"),
  dr_caption   = paste0("Robustness (Appendix, doubly-robust): Not-yet-treated",
                        " Control Group (Direct Frame)"),
  dr_label     = "tab:es_robust_notyet_dr",
  dr_note      = NOTE_NOTYET_DR
)

# ── D-F. Secondary outcomes ────────────────────────────────────────────────────

sec_specs <- list(
  list(res = r_corp, key = "corp", file = "es_corp",
       cap = "Event Study: Corporate BEV Registrations per 100k",
       lbl = "tab:es_corp",
       ylab = "ATT (corporate BEV registrations per 100k)",
       olbl = "corporate BEV new registrations per 100k"),
  list(res = r_priv, key = "priv", file = "es_priv",
       cap = "Event Study: Private BEV Registrations per 100k",
       lbl = "tab:es_priv",
       ylab = "ATT (private BEV registrations per 100k)",
       olbl = "private BEV new registrations per 100k"),
  list(res = r_ice, key = "ice", file = "es_ice",
       cap = "Placebo: ICE New Registrations per 100k (Direct Frame)",
       lbl = "tab:es_ice",
       ylab = "ATT (ICE new registrations per 100k)",
       olbl = "ICE new registrations per 100k (placebo)")
)

for (sp in sec_specs) {
  key <- sp$key; res <- sp$res
  if (is.null(res)) {
    cat(sprintf("  skip %s (estimation returned NULL)\n", key)); next
  }
  cat(sprintf("Writing %s\n", sp$file))
  emit_section(
    results = setNames(list(res), key),
    ems     = setNames("reg", key),
    graph_file  = file.path(out_dir, paste0(sp$file, ".png")),
    graph_title = sp$cap,
    main_stem    = file.path(out_dir, sp$file),
    main_caption = sp$cap,
    main_label   = sp$lbl,
    main_note    = note_sec(sp$olbl),
    ylab         = sp$ylab
  )
}

# ── G. ATT summary CSV (one row per estimated cell) ────────────────────────────

att_summary <- rbindlist(list(
  .att_row("dir_uncond",   "direct", "none",       "reg", "nevertreated",
           OUTCOME_BEV,  r_dir_unc),
  .att_row("dir_cond_reg", "direct", "XFORMLA_CS", "reg", "nevertreated",
           OUTCOME_BEV,  r_dir_reg),
  .att_row("dir_cond_dr",  "direct", "XFORMLA_CS", "dr",  "nevertreated",
           OUTCOME_BEV,  r_dir_dr),
  .att_row("brd_uncond",   "broad",  "none",       "reg", "nevertreated",
           OUTCOME_BEV,  r_brd_unc),
  .att_row("brd_cond_reg", "broad",  "XFORMLA_CS", "reg", "nevertreated",
           OUTCOME_BEV,  r_brd_reg),
  .att_row("brd_cond_dr",  "broad",  "XFORMLA_CS", "dr",  "nevertreated",
           OUTCOME_BEV,  r_brd_dr),
  .att_row("nyt_uncond",   "direct", "none",       "reg", "notyettreated",
           OUTCOME_BEV,  r_nyt_unc),
  .att_row("nyt_reg",      "direct", "XFORMLA_CS", "reg", "notyettreated",
           OUTCOME_BEV,  r_nyt_reg),
  .att_row("nyt_dr",       "direct", "XFORMLA_CS", "dr",  "notyettreated",
           OUTCOME_BEV,  r_nyt_dr),
  .att_row("corp",         "direct", "none",       "reg", "nevertreated",
           OUTCOME_CORP, r_corp),
  .att_row("priv",         "direct", "none",       "reg", "nevertreated",
           OUTCOME_PRIV, r_priv),
  .att_row("ice",          "direct", "none",       "reg", "nevertreated",
           OUTCOME_ICE,  r_ice)
), use.names = TRUE, fill = TRUE)

fwrite(att_summary, file.path(out_dir, "est_att_main.csv"))
cat(sprintf("\nAll DiD main outputs -> %s\n", out_dir))
