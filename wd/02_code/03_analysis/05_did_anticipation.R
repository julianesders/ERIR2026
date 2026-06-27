# ─────────────────────────────────────────────────────────────────────────────
# 05_did_anticipation.R  — CS anticipation robustness (k = 0, 1, 2)
#
# Reruns the headline direct-frame conditional DR specification from
# 04_did_main.R three times with 0, 1, and 2 periods of allowed anticipation
# (the `anticipation` argument in did::att_gt).
#
# Under anticipation = k, the estimator treats event times
# e ∈ {−k, …, −1} as potentially anticipatory and excludes them from the
# pre-trend Wald test. Only e ≤ −(k+1) count as "clean" pre-treatment.
# Comparing the three columns shows:
#   (a) whether pre-treatment movement is concentrated in the window
#       immediately before onset (true anticipation vs. pre-trend violation);
#   (b) whether post-treatment ATTs shift materially when the comparison
#       base period is moved back by one or two years.
#
# Estimation: conditional doubly-robust (est_method = "dr"), direct frame,
# never-treated control, XFORMLA_CS covariates — identical to section A of
# 04_did_main.R except for the `anticipation` argument.
#
# Outputs -> 03_output/05_did_anticipation/
#   es_anticipation.{png,tex,csv}   3-column event-study (k=0,1,2)
#   est_att_anticipation.csv        overall ATT per spec
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
  stop("Cannot determine script path. Run as: Rscript 05_did_anticipation.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "03_output", "05_did_anticipation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))
source(file.path(code_dir, "03_analysis", "_did_helpers.R"))
set.seed(202L)

# ── Extend COL_LABELS with anticipation-specific keys ─────────────────────────
# These are appended to the global COL_LABELS so es_graph and
# write_es_longtblr can resolve the panel/column headings from COL_LABELS[nm].

COL_LABELS <- c(COL_LABELS,
  antcp0 = "$\\delta=0$ (baseline)",
  antcp1 = "$\\delta=1$",
  antcp2 = "$\\delta=2$"
)

# ── Data ──────────────────────────────────────────────────────────────────────

dat <- .read_frame(file.path(data_final, "frame_did_direct.csv"))
stopifnot(is.integer(dat$ags8_id), is.integer(dat$gname_cs))
if (!"state_green_base_z" %in% names(dat))
  stop("frame missing state_green_base_z — re-run 00_prep_analysis.R")

cat(sprintf("Loaded direct frame: %d obs | %d AGS8 | %d ever-treated\n",
            nrow(dat), uniqueN(dat$AGS8),
            uniqueN(dat[gname_cs > 0L, AGS8])))

# ── Runner ────────────────────────────────────────────────────────────────────
# Mirrors run_spec() from _did_helpers.R but exposes `anticipation` as a
# parameter. Always uses conditional DR / never-treated; only the anticipation
# window varies across calls.

run_antcp_spec <- function(dat, yname, xformla, anticipation, label) {
  dat_y <- dat[!is.na(get(yname))]

  if (!is.null(xformla)) {
    cov_nms <- all.vars(xformla)
    missing_cols <- setdiff(cov_nms, names(dat_y))
    if (length(missing_cols)) {
      cat(sprintf("  [%s] MISSING covariate columns: %s — skip\n",
                  label, paste(missing_cols, collapse = ", ")))
      return(NULL)
    }
    complete_flag <- dat_y[, {
      ok <- TRUE
      for (cc in cov_nms) ok <- ok & any(!is.na(get(cc)))
      .(ok = ok)
    }, by = AGS8]
    ok_ags8 <- complete_flag[ok == TRUE, AGS8]
    n_drop  <- uniqueN(dat_y$AGS8) - length(ok_ags8)
    if (n_drop > 0L)
      cat(sprintf("  [%s] dropping %d AGS8 with NA in covariate(s)\n",
                  label, n_drop))
    dat_y <- dat_y[AGS8 %in% ok_ags8]
  }

  n_tr <- uniqueN(dat_y[gname_cs > 0L, AGS8])
  n_co <- uniqueN(dat_y[gname_cs == 0L, AGS8])

  cn <- dat_y[gname_cs > 0L, .(n = uniqueN(AGS8)),
              by = gname_cs][order(gname_cs)]
  cat(sprintf("  [%s] anticipation=%d  n_tr=%d  n_co=%d  es_max(tbd)\n",
              label, anticipation, n_tr, n_co))
  cat(sprintf("  [%s] cohort sizes: %s\n", label,
              paste(sprintf("%d(%d)", cn$gname_cs, cn$n), collapse = " ")))
  thin <- cn[n < 10L]
  if (nrow(thin))
    cat(sprintf("  [%s] WARNING: %d cohort(s) <10 treated: %s\n",
                label, nrow(thin), paste(thin$gname_cs, collapse = ", ")))

  es_max <- es_max_data_driven(dat_y, yname)
  cat(sprintf("  [%s] es_max=%d\n", label, es_max))

  cs_obj <- tryCatch(
    withCallingHandlers(
      att_gt(
        yname                  = yname,
        gname                  = "gname_cs",
        idname                 = "ags8_id",
        tname                  = "year",
        data                   = as.data.frame(dat_y),
        control_group          = "nevertreated",
        anticipation           = anticipation,
        xformla                = xformla,
        est_method             = "dr",
        clustervars            = "AGS5",
        bstrap                 = TRUE,
        biters                 = CS_BITERS,
        allow_unbalanced_panel = TRUE,
        print_details          = FALSE
      ),
      warning = function(w) {
        cat("  CS warn:", conditionMessage(w), "\n")
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      cat("  CS err:", conditionMessage(e), "\n"); NULL
    }
  )

  if (!is.null(cs_obj)) {
    n_na <- sum(is.na(cs_obj$att))
    if (n_na > 0L)
      cat(sprintf("  [%s] %d / %d ATT(g,t) are NA\n",
                  label, n_na, length(cs_obj$att)))
  }

  list(
    es_agg     = cs_es_agg(cs_obj, es_max),
    att_agg    = cs_att_agg(cs_obj),
    pre_tst    = cs_pre_test(cs_obj),
    n_treated  = n_tr,
    n_control  = n_co,
    est_method = "dr"
  )
}

# ── Estimate ───────────────────────────────────────────────────────────────────

cat("\n=== Anticipation robustness (direct frame, cond. dr) ===\n")

results <- list(
  antcp0 = run_antcp_spec(dat, OUTCOME_BEV, XFORMLA_CS, 0L, "antcp0"),
  antcp1 = run_antcp_spec(dat, OUTCOME_BEV, XFORMLA_CS, 1L, "antcp1"),
  antcp2 = run_antcp_spec(dat, OUTCOME_BEV, XFORMLA_CS, 2L, "antcp2")
)
ems <- c(antcp0 = "dr", antcp1 = "dr", antcp2 = "dr")

# ── Outputs ────────────────────────────────────────────────────────────────────

cat("\nWriting outputs...\n")

NOTE <- paste0(
  "\\textcite{callaway2021difference} doubly robust estimates for direct and broad ",
  "treatment accounting for $\\delta\\in \\{0,1,2\\}$ periods of anticipation. ",
  "Simultaneous 95\\% confidence bands (multiplier bootstrap, $B=", CS_BITERS, "$, ",
  "clustered at AGS5). Conditional specifications control for baseline tax capacity, ",
  "population density, and state Green vote share. ",
  "$p$-value of a joint pre-test for $H_0:\\text{ATT}(g,t)=0$ for all $g,t$ ",
  "in the last line. $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$."
)

emit_section(
  results      = results,
  ems          = ems,
  graph_file   = file.path(out_dir, "es_anticipation.png"),
  graph_title  = "Anticipation Robustness ($\\delta = 0, 1, 2$)",
  main_stem    = file.path(out_dir, "es_anticipation"),
  main_caption = "Anticipation Robustness: CS, Direct Frame ($\\delta = 0, 1, 2$)",
  main_label   = "tab:es_anticipation",
  main_note    = NOTE,
  ncol         = 3L,
  fig_caption  = "Anticipation Robustness: CS, Direct Frame ($\\delta = 0, 1, 2$)",
  fig_label    = "fig:es_anticipation",
  fig_note     = paste0(
    "CS conditional DR event study, direct frame, never-treated control. ",
    "Each panel allows $\\delta$ periods of pre-treatment anticipation ",
    "($\\delta=0$: baseline; $\\delta=1$: $e=-1$ anticipatory; ",
    "$\\delta=2$: $e\\in\\{-1,-2\\}$ anticipatory). ",
    "Simultaneous 95\\% CI bands (multiplier bootstrap, ",
    "$B=", CS_BITERS, "$, clustered at AGS5). ",
    "The dashed vertical line marks the last pre-treatment period."
  ),
  fig_width  = 12.0,
  fig_height =  4.5
)

# ATT summary CSV
att_rows <- rbindlist(lapply(names(results), function(nm) {
  r <- results[[nm]]
  if (is.null(r) || is.null(r$att_agg)) return(NULL)
  a <- r$att_agg
  k <- as.integer(sub("antcp", "", nm))
  pt <- r$pre_tst
  data.table(
    spec         = nm,
    anticipation = k,
    frame        = "direct",
    xformla      = "XFORMLA_CS",
    est_method   = "dr",
    control      = "nevertreated",
    outcome      = OUTCOME_BEV,
    att          = a$overall.att,
    se           = a$overall.se,
    ci_lo        = a$overall.att - 1.96 * a$overall.se,
    ci_hi        = a$overall.att + 1.96 * a$overall.se,
    pre_test_pval   = pt$pval,
    pre_test_method = pt$method,
    n_treated    = r$n_treated,
    n_control    = r$n_control
  )
}), use.names = TRUE, fill = TRUE)
fwrite(att_rows, file.path(out_dir, "est_att_anticipation.csv"))

cat(sprintf("\nAnticipation robustness outputs -> %s\n", out_dir))
