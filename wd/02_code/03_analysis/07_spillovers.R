# ─────────────────────────────────────────────────────────────────────────────
# 07_spillovers.R   — Spatial spillover robustness (donut)
#
# Uses the DIRECT frame (sharp Gemeinde-level treatment), where adjacency-based
# spillover is mechanically coherent: the spillover channel is a directly-treated
# 1st-order Gemeinde neighbour. The broad frame is NOT used — broad treatment is
# Kreis-level coverage and does not map onto Gemeinde adjacency.
#
# Donut robustness: the main direct DiD spec re-estimated after DROPPING
# never-treated units whose 1st-order Gemeinde neighbour was directly treated
# (direct_treated_any_nbrs_gem_1 == 1). Tests whether spillover into spatially
# adjacent controls drives the headline ATT.
#
# Outputs (03_output/07_spillovers/):
#   es_donut.{png,tex,csv}
#   est_att_spillover.csv   (same schema as est_att_main.csv)
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
  stop("Cannot determine script path. Run as: Rscript 07_spillovers.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "03_output", "07_spillovers")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))
source(file.path(code_dir, "03_analysis", "_did_helpers.R"))
set.seed(202L)

frame_direct <- .read_frame(file.path(data_final, "frame_did_direct.csv"))

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Donut robustness — direct frame, contamination-trimmed control pool
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
r_donut <- run_spec(dat_donut, OUTCOME_BEV, XFORMLA_CS, "nevertreated",
                    "dr", "donut")

NOTE_DONUT <- paste0(
  "\\textcite{callaway2021difference} conditional doubly robust estimates dropping ",
  "never-treated municipalities with direct-treated 1st-degree neighbors from the ",
  "control group. The dependent variable is BEV registration per 100k inhabitants. ",
  "Controls include state Green vote share, tax capacity p.c., and population density. ",
  "Control group is never-treated municipalities. SEs are AGS5-clustered derived from ",
  "a multiplier boostrap ($B=", CS_BITERS, "$). $p$-value of a joint pre-test for ",
  "$H_0:\\text{ATT}(g,t)=0$ for all $g,t$ in the last line."
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
# 2. ATT summary (same schema as est_att_main.csv)
# ═══════════════════════════════════════════════════════════════════════════════

att_summary <- rbindlist(list(
  .att_row("donut", "direct", "XFORMLA_CS", "dr", "nevertreated",
           OUTCOME_BEV, r_donut)
), use.names = TRUE, fill = TRUE)
fwrite(att_summary, file.path(out_dir, "est_att_spillover.csv"))

cat(sprintf("\nSpillover outputs -> %s\n", out_dir))
