# ─────────────────────────────────────────────────────────────────────────────
# 06_spillovers.R   — Spatial spillovers (donut + descriptive event study)
#
# Both probes use the DIRECT frame (sharp Gemeinde-level treatment), where
# adjacency-based spillover is mechanically coherent: the spillover channel is a
# directly-treated 1st-order Gemeinde neighbour. The broad frame is NOT used for
# spillovers — broad treatment is Kreis-level coverage and does not map onto
# Gemeinde adjacency.
#
# Both probes go through the SAME CS machinery as 03_did_main.R
# (run_spec / emit_section from _did_helpers.R), so the outputs are byte-for-byte
# the same format as the main DiD specification — only the SAMPLE differs:
#
#   - Donut robustness: the main direct DiD spec (direct frame / never-treated /
#     unconditional CS — reproduces 03_did_main.R section A) re-estimated after
#     DROPPING never-treated units whose 1st-order Gemeinde neighbour was
#     directly treated (direct_treated_any_nbrs_gem_1 == 1).
#   - Descriptive spillover event study: among never-direct-treated units,
#     pseudo-treatment = first year direct_treated_any_nbrs_gem_1 == 1;
#     CS event study (descriptive only — selection on geography).
#
# Outputs (04_results/06_spillovers/), identical format to the main spec:
#   es_donut.{png,tex,csv}
#   es_spillover.{png,tex,csv}
#   est_att_spillover.csv   (one row per probe, same schema as est_att_main.csv)
# ─────────────────────────────────────────────────────────────────────────────

library(data.table)
library(did)
library(ggplot2)

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
  stop("Cannot determine script path. Run as: Rscript 06_spillovers.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "04_results", "06_spillovers")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))
source(file.path(code_dir, "03_analysis", "_did_helpers.R"))

frame_direct <- .read_frame(file.path(data_final, "frame_did_direct.csv"))

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Donut robustness — reproduces 03 section A on the contamination-trimmed pool
# ═══════════════════════════════════════════════════════════════════════════════

# Drop never-treated controls whose 1st-order Gemeinde neighbour was directly
# treated. This is the standard donut: it tests whether spillover into spatially
# adjacent controls drives the headline direct-treatment ATT.
contaminated <- frame_direct[
  is.na(first_treat_direct) & direct_treated_any_nbrs_gem_1 == 1L,
  unique(AGS8)
]
cat(sprintf("Donut: dropping %d never-treated units with treated neighbours\n",
            length(contaminated)))
dat_donut <- frame_direct[!(AGS8 %in% contaminated)]

cat("\n=== 1. Donut (direct frame, contaminated controls dropped) ===\n")
r_donut <- run_spec(dat_donut, OUTCOME_BEV, NULL, "nevertreated", "reg", "donut")

NOTE_DONUT <- paste0(
  "Callaway--Sant'Anna (2021) CS estimator, unconditional. Donut robustness: ",
  "the main direct-frame spec (reproduces \\texttt{03\\_did\\_main.R} section A) ",
  "re-estimated after dropping never-treated controls with a directly-treated ",
  "1st-order Gemeinde neighbour ",
  "(\\texttt{direct\\_treated\\_any\\_nbrs\\_gem\\_1}=1). ",
  "Outcome: BEV new registrations per 100,000 (level; winsorised at 99th pct). ",
  "Bootstrap SEs clustered at county (AGS5), $B=", CS_BITERS, "$. ",
  "Graph CI bands: simultaneous (95\\%). Table: pointwise SEs; ",
  "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$ (two-sided). ",
  "Pre-treatment test: the package's joint Wald $\\chi^2$ of $H_0$ that all ",
  "pre-treatment $\\text{ATT}(g,t)=0$ (AGS5-clustered), ``--'' if unavailable."
)

if (!is.null(r_donut))
  emit_section(
    results = list(donut = r_donut),
    ems     = c(donut = "reg"),
    graph_file  = file.path(out_dir, "es_donut.png"),
    graph_title = "Donut Robustness: BEV New Registrations (Direct Frame)",
    main_stem    = file.path(out_dir, "es_donut"),
    main_caption = paste0("Donut Robustness: BEV New Registrations per 100k",
                          " (Direct Frame, Contaminated Controls Dropped)"),
    main_label   = "tab:es_donut",
    main_note    = NOTE_DONUT
  )

# ═══════════════════════════════════════════════════════════════════════════════
# 2. Descriptive spillover event study — pseudo-treatment from a treated neighbour
# ═══════════════════════════════════════════════════════════════════════════════

# Among never-direct-treated units, pseudo-treatment = first year a 1st-order
# Gemeinde neighbour appears as directly treated. Never-pseudo-treated units (no
# neighbour ever directly treated) form the control group. Same CS machinery as
# the main spec, driven by the pseudo-cohort column `gname_pseudo`.
never <- frame_direct[is.na(first_treat_direct)]
first_nbr <- never[
  direct_treated_any_nbrs_gem_1 == 1L,
  .(pseudo_treat = min(year)), by = AGS8
]
sp <- merge(never, first_nbr, by = "AGS8", all.x = TRUE)
sp[, gname_pseudo := fifelse(is.na(pseudo_treat), 0L, as.integer(pseudo_treat))]
sp[, ags8_id := .GRP, by = AGS8]

cat("\n=== 2. Descriptive spillover ES (direct-neighbour pseudo-treatment) ===\n")
r_spill <- run_spec(sp, OUTCOME_BEV, NULL, "nevertreated", "reg", "spillover",
                    gname = "gname_pseudo")

NOTE_SPILL <- paste0(
  "Callaway--Sant'Anna (2021) CS estimator, unconditional. \\textbf{Descriptive ",
  "only.} Among never-direct-treated Gemeinden, pseudo-treatment = first year a ",
  "1st-order Gemeinde neighbour is directly treated ",
  "(\\texttt{direct\\_treated\\_any\\_nbrs\\_gem\\_1}); never-pseudo-treated ",
  "units are the control group. The pseudo-cohort is endogenous to spatial ",
  "geography, so this is NOT causal (selection on geography). ",
  "Outcome: BEV new registrations per 100,000 (level; winsorised at 99th pct). ",
  "Bootstrap SEs clustered at county (AGS5), $B=", CS_BITERS, "$. ",
  "Graph CI bands: simultaneous (95\\%). Table: pointwise SEs; ",
  "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$ (two-sided). ",
  "Pre-treatment test: the package's joint Wald $\\chi^2$, ``--'' if unavailable."
)

if (!is.null(r_spill))
  emit_section(
    results = list(spillover = r_spill),
    ems     = c(spillover = "reg"),
    graph_file  = file.path(out_dir, "es_spillover.png"),
    graph_title = "Descriptive Spillover ES (Direct-neighbour Pseudo-treatment)",
    main_stem    = file.path(out_dir, "es_spillover"),
    main_caption = paste0("Descriptive Spillover Event Study: BEV New ",
                          "Registrations per 100k (Never-treated, ",
                          "Direct-neighbour Pseudo-treatment)"),
    main_label   = "tab:es_spillover",
    main_note    = NOTE_SPILL,
    ylab         = "Pseudo-ATT (BEV new registrations per 100k)"
  )

# ═══════════════════════════════════════════════════════════════════════════════
# 3. ATT summary (same schema as est_att_main.csv)
# ═══════════════════════════════════════════════════════════════════════════════

att_summary <- rbindlist(list(
  .att_row("donut",     "direct",            "none", "reg", "nevertreated",
           OUTCOME_BEV, r_donut),
  .att_row("spillover", "direct_pseudo_nbr", "none", "reg", "nevertreated",
           OUTCOME_BEV, r_spill)
), use.names = TRUE, fill = TRUE)
fwrite(att_summary, file.path(out_dir, "est_att_spillover.csv"))

cat(sprintf("\nSpillover outputs -> %s\n", out_dir))
