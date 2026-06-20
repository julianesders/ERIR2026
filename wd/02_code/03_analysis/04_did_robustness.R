# ───────────────────────────────────────────────────────────────────────────────
# 04_did_robustness.R   — Robustness for the main DiD result
#
# Headline outcome: bev_neuzulassungen_p100k (level + log1p), broad frame.
# Each variant runs BJS (static) + CS (simple ATT). PPML/count specifications
# have been dropped: the bare N_elektro_* columns are not in the panel, and
# the per-100k rate already gives the substantive answer.
#
# Variants:
#   baseline           anticip 0, no xformla              (matches 03 headline)
#   conditional        anticip 0, xformla=XFORMLA_CS      (covariate-adjusted)
#   anticipation_1     anticip 1, no xformla
#   drop_2016          anticip 0, exclude 2016 cohort
#   drop_covid         anticip 0, exclude 2020-21         (short-run only)
#   complete_case      anticip 0, drop rows where the outcome's _imp == TRUE
#   notyet_evertreated anticip 0, ever-treated only, control=notyettreated
#                                                         (selection probe)
#
# Outputs (04_results/04_did_robustness/):
#   est_att_robust.csv          one row per variant × scale × estimator
#   fig_robust_grid.pdf         forest: ATT × variant, by scale
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
  stop("Cannot determine script path. Run as: Rscript 04_did_robustness.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "04_results", "04_did_robustness")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))

.read_frame <- function(p) fread(p,
  colClasses = list(character = c("AGS8", "AGS5", "AGS2", "treat_type")))
frame_broad <- .read_frame(file.path(data_final, "frame_did_broad.csv"))

OUTCOMES <- c(
  level = "bev_neuzulassungen_p100k",
  log1p = "log1p_bev_neuzulassungen_p100k"
)
IMP_FLAG <- "N_elektro_overall_imp"

# -- Variant definitions -------------------------------------------------------

variants <- list(
  baseline = list(
    filter = function(d) d, anticip = 0L, xformla = NULL,
    control = "nevertreated", scope = "full"
  ),
  conditional = list(
    filter = function(d) d, anticip = 0L, xformla = XFORMLA_CS,
    control = "nevertreated", scope = "full"
  ),
  anticipation_1 = list(
    filter = function(d) d, anticip = 1L, xformla = NULL,
    control = "nevertreated", scope = "full"
  ),
  drop_2016 = list(
    filter = function(d) d[is.na(first_treat_broad) | first_treat_broad != 2016L],
    anticip = 0L, xformla = NULL,
    control = "nevertreated", scope = "full"
  ),
  drop_covid = list(
    filter = function(d) d[!(year %in% c(2020L, 2021L))],
    anticip = 0L, xformla = NULL,
    control = "nevertreated", scope = "short_run_only"
  ),
  complete_case = list(
    filter = function(d) {
      if (!(IMP_FLAG %in% names(d))) return(d)
      d[get(IMP_FLAG) == FALSE]
    },
    anticip = 0L, xformla = NULL,
    control = "nevertreated", scope = "full"
  ),
  notyet_evertreated = list(
    filter = function(d) d[!is.na(first_treat_broad)],
    anticip = 0L, xformla = NULL,
    control = "notyettreated", scope = "selection_probe"
  )
)

# -- Estimators ---------------------------------------------------------------

bjs_simple <- function(dat, yname) {
  r <- tryCatch(as.data.table(did_imputation(
    data        = as.data.frame(dat),
    yname       = yname,
    gname       = "gname_bjs",
    tname       = "year",
    idname      = "ags8_id",
    horizon     = NULL,
    cluster_var = "AGS5"
  )), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  row <- r[term == "treat"][1L]
  if (is.null(row) || nrow(row) == 0L) return(NULL)
  list(estimate = row$estimate, se = row$std.error)
}

cs_simple <- function(dat, yname, anticip, xformla, control_group) {
  cs <- tryCatch(att_gt(
    yname = yname, gname = "gname_cs", idname = "ags8_id", tname = "year",
    data = as.data.frame(dat), control_group = control_group,
    anticipation = anticip, xformla = xformla,
    est_method = "dr",
    clustervars = "AGS5", bstrap = TRUE, biters = 999L,
    allow_unbalanced_panel = TRUE,
    print_details = FALSE
  ), error = function(e) {
    cat("    CS err:", conditionMessage(e), "\n"); NULL })
  if (is.null(cs)) return(NULL)
  s <- tryCatch(aggte(cs, type = "simple", na.rm = TRUE),
                error = function(e) NULL)
  if (is.null(s)) return(NULL)
  list(estimate = s$overall.att, se = s$overall.se)
}

# -- Run loop -----------------------------------------------------------------

rows <- list()

for (vname in names(variants)) {
  v <- variants[[vname]]
  dat <- v$filter(frame_broad)
  n_treated <- uniqueN(dat[gname_cs > 0L, AGS8])
  cat(sprintf("\n=== %s (scope=%s, n=%d, treated=%d) ===\n",
              vname, v$scope, nrow(dat), n_treated))

  for (sc in names(OUTCOMES)) {
    yname <- OUTCOMES[[sc]]
    if (!(yname %in% names(dat))) {
      cat(sprintf("  skip: %s missing\n", yname)); next
    }
    dat_y <- dat[!is.na(get(yname))]

    b <- bjs_simple(dat_y, yname)
    if (!is.null(b))
      rows[[length(rows) + 1L]] <- data.table(
        variant = vname, scope = v$scope, scale = sc,
        outcome = yname, estimator = "BJS",
        estimate = b$estimate, se = b$se
      )

    c <- cs_simple(dat_y, yname, v$anticip, v$xformla, v$control)
    if (!is.null(c))
      rows[[length(rows) + 1L]] <- data.table(
        variant = vname, scope = v$scope, scale = sc,
        outcome = yname, estimator = "CS",
        estimate = c$estimate, se = c$se
      )
  }
}

robust_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
robust_dt[, ci_lo := estimate - 1.96 * se]
robust_dt[, ci_hi := estimate + 1.96 * se]
fwrite(robust_dt, file.path(out_dir, "est_att_robust.csv"))

# -- Forest figure: 2 panels (level, log1p), variants on y-axis ---------------

if (nrow(robust_dt) > 0L) {
  pdt <- copy(robust_dt)
  pdt[, label := paste(variant, estimator, sep = " / ")]
  pdt[, label := factor(label, levels = unique(label))]
  p <- ggplot(pdt, aes(x = label, y = estimate, color = scope)) +
    geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi)) +
    facet_wrap(~ scale, scales = "free_x") +
    coord_flip() +
    labs(x = NULL, y = "ATT (BEV new registrations per 100k)",
         color = "Scope",
         caption = "drop_covid is short-run only; notyet_evertreated is a selection probe.") +
    theme_minimal(base_size = 11)
  ggsave(file.path(out_dir, "fig_robust_grid.pdf"), p,
         width = 9, height = 5.5)
}

cat(sprintf("\nDiD robustness outputs -> %s\n", out_dir))
