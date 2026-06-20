# ───────────────────────────────────────────────────────────────────────────────
# 06_spillovers.R   — Spatial spillovers (donut + descriptive event study)
#
#   - Donut robustness: rerun the main DiD spec (broad / never / CS-dr) but
#     DROP never-treated units flagged by either direct_treated_any_nbrs_gem_1
#     or broad_treated_any_nbrs_kreis (more conservative donut).
#     Writes est_att_donut.csv; 07_assemble.R picks it up by pattern.
#   - Descriptive spillover event study: among never-treated units, set
#     pseudo-treatment to the first year that emk_absorbing_any_nbrs_1 == 1;
#     estimate a BJS event study.
#
# Outputs (04_results/06_spillovers/):
#   est_att_donut.csv
#   es_spillover.pdf
# ───────────────────────────────────────────────────────────────────────────────

library(data.table)
library(did)
library(didimputation)
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

OUTCOME    <- "bev_neuzulassungen_p100k"
CS_BITERS  <- 1999L

# Match the headline 03 spec: unconditional CS (xformla = NULL).
xformla_cov <- NULL

.read_frame <- function(p) fread(p,
  colClasses = list(character = c("AGS8", "AGS5", "AGS2", "treat_type")))
frame_broad <- .read_frame(file.path(data_final, "frame_did_broad.csv"))

# -- 1. Donut robustness ------------------------------------------------------
# Drop never-treated controls with any 1st-order neighbour that was treated
# (broad). Uses `emk_absorbing_any_nbrs_1` from spatial_neighbors_ags8.csv —
# the granular Gemeinde / Kreis split indicators in the python script revision
# are not yet in the merged frame.

contaminated <- frame_broad[
  is.na(first_treat_broad) & emk_absorbing_any_nbrs_1 == 1L,
  unique(AGS8)
]
cat(sprintf("Donut: dropping %d never-treated units with treated neighbours\n",
            length(contaminated)))

dat_clean <- frame_broad[!(AGS8 %in% contaminated)]

cs <- tryCatch(att_gt(
  yname = OUTCOME, gname = "gname_cs", idname = "ags8_id", tname = "year",
  data = as.data.frame(dat_clean), control_group = "nevertreated",
  anticipation = 0L, xformla = xformla_cov, est_method = "dr",
  clustervars = "AGS5", bstrap = TRUE, biters = CS_BITERS,
  allow_unbalanced_panel = TRUE,
  print_details = FALSE
), error = function(e) {
  cat("  CS err:", conditionMessage(e), "\n"); NULL })

donut_rows <- list()
if (!is.null(cs)) {
  s <- tryCatch(aggte(cs, type = "simple", na.rm = TRUE),
                error = function(e) NULL)
  if (!is.null(s)) {
    donut_rows[[length(donut_rows) + 1L]] <- data.table(
      spec = "donut_drop_treated_neighbours",
      estimator = "CS-dr",
      outcome = OUTCOME,
      estimate = s$overall.att,
      se       = s$overall.se,
      ci_lo    = s$overall.att - 1.96 * s$overall.se,
      ci_hi    = s$overall.att + 1.96 * s$overall.se,
      n_treated = uniqueN(dat_clean[gname_cs > 0L, AGS8]),
      n_control = uniqueN(dat_clean[gname_cs == 0L, AGS8])
    )
  }
}
donut_dt <- rbindlist(donut_rows, fill = TRUE)
fwrite(donut_dt, file.path(out_dir, "est_att_donut.csv"))

# -- 2. Descriptive spillover event study -------------------------------------
# Among never-treated units, pseudo-treatment = first year a 1st-order
# neighbour appears as broad-treated. BJS event study on the outcome.

never <- frame_broad[is.na(first_treat_broad)]
first_nbr <- never[
  emk_absorbing_any_nbrs_1 == 1L,
  .(pseudo_treat = min(year)), by = AGS8
]
sp <- merge(never, first_nbr, by = "AGS8", all.x = TRUE)
# didimputation v0.5.1: never-treated coded as 0 (not Inf).
sp[, gname_pseudo := fifelse(is.na(pseudo_treat), 0L, as.integer(pseudo_treat))]
sp[, ags8_id := .GRP, by = AGS8]

n_pseudo <- uniqueN(sp[gname_pseudo > 0L, AGS8])
cat(sprintf("Spillover ES: %d never-treated AGS8 with treated neighbour\n",
            n_pseudo))

bjs <- tryCatch(did_imputation(
  data        = as.data.frame(sp),
  yname       = OUTCOME,
  gname       = "gname_pseudo",
  tname       = "year",
  idname      = "ags8_id",
  horizon     = TRUE,
  pretrends   = TRUE,
  cluster_var = "AGS5"
), error = function(e) {
  cat("  BJS err:", conditionMessage(e), "\n"); NULL })

if (!is.null(bjs)) {
  dt <- as.data.table(bjs)
  dt[, e := NA_integer_]
  dt[grepl("^tau[0-9]+$", term), e := as.integer(sub("tau", "", term))]
  dt[grepl("^pre[0-9]+$", term), e := -as.integer(sub("pre", "", term))]
  dt[is.na(e) & grepl("^-?[0-9]+$", term), e := as.integer(term)]
  dt <- dt[!is.na(e) & e >= ES_MIN & e <= ES_MAX]
  if (nrow(dt) > 0L) {
    dt[, ci_lo := estimate - 1.96 * std.error]
    dt[, ci_hi := estimate + 1.96 * std.error]
    p <- ggplot(dt, aes(x = e, y = estimate)) +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
      geom_vline(xintercept = -0.5, linetype = "dashed",
                 color = "grey50", linewidth = 0.4) +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                  alpha = 0.18, fill = "#7b3294") +
      geom_line(linewidth = 0.7, color = "#7b3294") +
      geom_point(size = 1.8, color = "#7b3294") +
      scale_x_continuous(breaks = ES_MIN:ES_MAX) +
      labs(
        x = "Years rel. to neighbour's first treatment",
        y = "Pseudo-ATT (BEV new registrations p100k)",
        caption = paste0("Descriptive only — selection-on-geography ",
                         "caveat applies (never-treated sample).")
      ) +
      theme_minimal(base_size = 11)
    ggsave(file.path(out_dir, "es_spillover.pdf"), p,
           width = 7, height = 5)
    fwrite(dt, file.path(out_dir, "es_spillover.csv"))
  }
}

cat(sprintf("\nSpillover outputs -> %s\n", out_dir))
