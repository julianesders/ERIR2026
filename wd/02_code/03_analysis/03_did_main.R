# ───────────────────────────────────────────────────────────────────────────────
# 03_did_main.R   — Part (ii): adoption response to EMK funding
#
# Replaces the old 03_spillover_did.R. Two co-primary staggered-DiD estimators:
#   - Borusyak-Jaravel-Spiess (BJS) via didimputation (efficiency, short T)
#   - Callaway-Sant'Anna (CS) via did, with `est_method="dr"` and a baseline
#     covariate vector (covariate-adjusted benchmark)
#
# Designs: {broad, direct} × {never-treated control, not-yet-treated control}.
# Main spec: broad / never.   Sharp spec: direct / never.
# Selection check: notyet with `anticipation = 1`.
#
# Outputs (04_results/03_did_main/):
#   tab_att_main.{tex,csv}                  ATT × outcome × frame × estimator
#   est_att_main.csv                        long form for the manifest
#   es_<outcome>_<frame>.pdf                event studies (BJS + CS overlay)
#   pretests.csv                            BJS joint pretrends + CS Wald
# ───────────────────────────────────────────────────────────────────────────────

library(data.table)
library(fixest)
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
  stop("Cannot determine script path. Run as: Rscript 03_did_main.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "04_results", "03_did_main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))

# -- Constants -----------------------------------------------------------------

CS_BITERS <- 1999L

OUTCOMES <- c(
  "bev_neuzulassungen_p100k",   # main
  "bev_corporate_p100k",
  "bev_private_p100k",
  "bev_stock_p100k",            # cumulative, secondary
  "ice_neuzulassungen_p100k"    # placebo (from counts)
)

# CS covariate-adjusted spec: baseline (z'd) covariates, time-invariant.
xformla_cov <- ~ kk_base_z + sk_base_z + dens_base_z + green_base_z +
                  bev_base_z + chg_base_z

# -- Load frames ---------------------------------------------------------------

frames <- list(
  broad  = readRDS(file.path(data_final, "frame_did_broad.rds")),
  direct = readRDS(file.path(data_final, "frame_did_direct.rds"))
)

# Sanity: ags8_id integer, gname_cs integer, gname_bjs numeric (already set
# in 00_prep_analysis.R, but assert)
for (nm in names(frames)) {
  d <- frames[[nm]]
  stopifnot(is.integer(d$ags8_id),
            is.integer(d$gname_cs),
            is.numeric(d$gname_bjs))
}

# -- Helpers -------------------------------------------------------------------

run_cs <- function(yname, dat, control_group, anticipation = 0L,
                   xformla = NULL) {
  att_gt(
    yname                  = yname,
    gname                  = "gname_cs",
    idname                 = "ags8_id",
    tname                  = "year",
    data                   = as.data.frame(dat),
    control_group          = control_group,
    anticipation           = anticipation,
    xformla                = xformla,
    # est_method "dr" is the default and works with or without xformla;
    # IPW with NULL xformla is degenerate, so we always pick dr.
    est_method             = "dr",
    clustervars            = "AGS5",
    bstrap                 = TRUE,
    biters                 = CS_BITERS,
    # Outcomes only span 2014–2023; panel includes 1995–2023. Without this,
    # att_gt drops every unit because of missing pre-2014 outcomes.
    allow_unbalanced_panel = TRUE,
    print_details          = FALSE
  )
}

# horizon = TRUE -> event study + pretrends.  horizon = NULL -> static effect.
run_bjs <- function(yname, dat, horizon = TRUE) {
  did_imputation(
    data        = as.data.frame(dat),
    yname       = yname,
    gname       = "gname_bjs",
    tname       = "year",
    idname      = "ags8_id",
    horizon     = horizon,
    pretrends   = horizon,
    cluster_var = "AGS5"
  )
}

cs_to_dt <- function(es, label) {
  crit <- if (!is.null(es$crit.val.egt)) es$crit.val.egt else 1.96
  data.table(
    e          = es$egt,
    estimate   = es$att.egt,
    se         = es$se.egt,
    ci_lo_pw   = es$att.egt - 1.96 * es$se.egt,
    ci_hi_pw   = es$att.egt + 1.96 * es$se.egt,
    ci_lo_sim  = es$att.egt - crit * es$se.egt,
    ci_hi_sim  = es$att.egt + crit * es$se.egt,
    estimator  = label
  )
}

# BJS term parser (reuses the verbatim logic from the old 03_spillover_did.R)
bjs_to_dt <- function(res, label) {
  dt <- as.data.table(res)
  dt[, e := NA_integer_]
  dt[grepl("^tau[0-9]+$", term), e := as.integer(sub("tau", "", term))]
  dt[grepl("^pre[0-9]+$", term), e := -as.integer(sub("pre", "", term))]
  dt[is.na(e) & grepl("^-?[0-9]+$", term), e := as.integer(term)]
  dt <- dt[!is.na(e) & e >= ES_MIN & e <= ES_MAX]
  if (nrow(dt) == 0L) return(NULL)
  dt[, .(
    e,
    estimate  = estimate,
    se        = std.error,
    ci_lo_pw  = estimate - 1.96 * std.error,
    ci_hi_pw  = estimate + 1.96 * std.error,
    ci_lo_sim = estimate - 1.96 * std.error,  # BJS no simultaneous bands
    ci_hi_sim = estimate + 1.96 * std.error,
    estimator = label
  )]
}

# -- Design grid ---------------------------------------------------------------
# Plan: main = broad / never; sharp = direct / never; selection = broad / notyet
# with anticipation = 1. Each cell estimates both BJS and CS.

design <- rbindlist(list(
  data.table(frame = "broad",  ctrl = "never",  anticip = 0L, role = "main"),
  data.table(frame = "direct", ctrl = "never",  anticip = 0L, role = "sharp"),
  data.table(frame = "broad",  ctrl = "notyet", anticip = 1L, role = "selection")
))

att_rows  <- list()
pretests  <- list()

for (i in seq_len(nrow(design))) {
  fr_nm <- design$frame[i]
  ctrl  <- design$ctrl[i]
  antic <- design$anticip[i]
  role  <- design$role[i]
  dat   <- frames[[fr_nm]]
  cg    <- if (ctrl == "never") "nevertreated" else "notyettreated"

  n_treated <- uniqueN(dat[gname_cs > 0L, AGS8])
  n_control <- uniqueN(dat[gname_cs == 0L, AGS8])
  cat(sprintf("\n=== %s frame | %s ctrl | anticip=%d (%s) ===\n",
              fr_nm, ctrl, antic, role))
  cat(sprintf("  %d treated AGS8, %d control AGS8\n", n_treated, n_control))

  for (yname in OUTCOMES) {
    # Restrict to years where this outcome is observed. KBA outcomes start in
    # 2014; charging stations slightly later. Keeping pre-2014 rows in `dat`
    # makes BJS and CS scan 312k rows per call (3.5× slower) for no gain.
    dat_y <- dat[!is.na(get(yname))]
    cat(sprintf("  outcome: %s — n=%d (non-NA), years %d–%d\n",
                yname, nrow(dat_y),
                dat_y[, min(year)], dat_y[, max(year)]))

    # --- BJS ---
    bjs_es  <- tryCatch(run_bjs(yname, dat_y, horizon = TRUE),
                        error = function(e) {
                          cat("    BJS ES error:", conditionMessage(e), "\n"); NULL })
    # Static effect: horizon = NULL (per didimputation v0.5.1 docs).
    bjs_att <- tryCatch(run_bjs(yname, dat_y, horizon = NULL),
                        error = function(e) {
                          cat("    BJS ATT error:", conditionMessage(e), "\n"); NULL })

    if (!is.null(bjs_att)) {
      att_dt <- as.data.table(bjs_att)
      # Static-effect term label per the v0.5.1 docs is "treat".
      att_row <- att_dt[term == "treat"]
      if (nrow(att_row) >= 1L) {
        r <- att_row[1L]
        att_rows[[length(att_rows) + 1L]] <- data.table(
          frame = fr_nm, ctrl = ctrl, role = role, outcome = yname,
          estimator = "BJS",
          estimate  = r$estimate, se = r$std.error,
          ci_lo = r$estimate - 1.96 * r$std.error,
          ci_hi = r$estimate + 1.96 * r$std.error,
          n_treated = n_treated, n_control = n_control
        )
      }
    }

    # BJS joint pretrends p-value is not exposed by didimputation; pre-trend
    # evidence comes from the per-term BJS event-study leads in the plot and
    # from the CS Wald test below. No BJS pretest row written.

    # --- CS (covariate-adjusted) ---
    cs_obj <- tryCatch(run_cs(yname, dat_y, cg, anticip, xformla_cov),
                       error = function(e) {
                         cat("    CS error:", conditionMessage(e), "\n"); NULL })
    cs_simple <- if (!is.null(cs_obj))
      tryCatch(aggte(cs_obj, type = "simple", na.rm = TRUE),
               error = function(e) NULL) else NULL
    cs_es <- if (!is.null(cs_obj))
      tryCatch(aggte(cs_obj, type = "dynamic",
                     min_e = ES_MIN, max_e = ES_MAX,
                     # balance_e is numeric-or-NULL; NULL = no balancing.
                     balance_e = NULL, na.rm = TRUE, cband = TRUE),
               error = function(e) NULL) else NULL

    if (!is.null(cs_simple)) {
      att_rows[[length(att_rows) + 1L]] <- data.table(
        frame = fr_nm, ctrl = ctrl, role = role, outcome = yname,
        estimator = "CS-dr",
        estimate  = cs_simple$overall.att,
        se        = cs_simple$overall.se,
        ci_lo     = cs_simple$overall.att - 1.96 * cs_simple$overall.se,
        ci_hi     = cs_simple$overall.att + 1.96 * cs_simple$overall.se,
        n_treated = n_treated, n_control = n_control
      )
    }

    # CS Wald pre-test (W and Wpval, if present in att_gt object)
    if (!is.null(cs_obj)) {
      pretests[[length(pretests) + 1L]] <- data.table(
        frame = fr_nm, ctrl = ctrl, role = role, outcome = yname,
        estimator = "CS-dr",
        wald_W     = if (!is.null(cs_obj$Wpval)) cs_obj$W      else NA_real_,
        wald_pval  = if (!is.null(cs_obj$Wpval)) cs_obj$Wpval  else NA_real_
      )
    }

    # --- Event-study plot ---
    parts <- Filter(Negate(is.null), list(
      if (!is.null(cs_es))  cs_to_dt(cs_es,    "CS-dr") else NULL,
      if (!is.null(bjs_es)) bjs_to_dt(bjs_es,  "BJS")   else NULL
    ))
    if (length(parts) > 0L) {
      plot_dt <- rbindlist(parts, use.names = TRUE, fill = TRUE)
      out_lbl <- if (yname %in% names(OUTCOME_LABELS)) OUTCOME_LABELS[[yname]] else yname
      ttl     <- sprintf("%s — %s frame, %s control", out_lbl, fr_nm, ctrl)
      p_es <- ggplot(plot_dt, aes(x = e, color = estimator, fill = estimator)) +
        geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
        geom_vline(xintercept = -0.5, linetype = "dashed",
                   color = "grey50", linewidth = 0.4) +
        geom_ribbon(aes(ymin = ci_lo_pw, ymax = ci_hi_pw),
                    alpha = 0.12, color = NA) +
        geom_line(aes(y = estimate), linewidth = 0.7) +
        geom_point(aes(y = estimate), size = 1.8) +
        scale_x_continuous(breaks = ES_MIN:ES_MAX) +
        labs(title = ttl, x = "Years relative to first funding receipt",
             y = "ATT", color = NULL, fill = NULL,
             caption = "Shaded band: pointwise 95% CI.") +
        theme_minimal(base_size = 11) +
        theme(legend.position = "bottom")
      fname <- sprintf("es_%s_%s_%s.pdf", yname, fr_nm, ctrl)
      ggsave(file.path(out_dir, fname), p_es, width = 7, height = 5)
    }
  }
}

# -- Write summary tables -----------------------------------------------------

att_dt  <- rbindlist(att_rows,  use.names = TRUE, fill = TRUE)
fwrite(att_dt,   file.path(out_dir, "est_att_main.csv"))

pre_dt  <- rbindlist(pretests,  use.names = TRUE, fill = TRUE)
fwrite(pre_dt,   file.path(out_dir, "pretests.csv"))

# Compact LaTeX summary (one row per outcome × role; BJS | CS-dr columns)
tab <- dcast(att_dt[role %in% c("main", "sharp")],
             outcome + frame + ctrl + role ~ estimator,
             value.var = c("estimate", "se"))
fwrite(tab, file.path(out_dir, "tab_att_main.csv"))

tex_lines <- c(
  "\\begin{tabular}{lllrrrr}",
  "\\hline",
  "Outcome & Frame & Role & BJS ATT (se) & CS-dr ATT (se) \\\\",
  "\\hline",
  tab[, sprintf("%s & %s & %s & %.3f (%.3f) & %.3f (%.3f) \\\\",
                outcome, frame, role,
                estimate_BJS, se_BJS,
                `estimate_CS-dr`, `se_CS-dr`)],
  "\\hline",
  "\\end{tabular}"
)
writeLines(tex_lines, file.path(out_dir, "tab_att_main.tex"))

cat(sprintf("\nMain DiD outputs -> %s\n", out_dir))
