# ─────────────────────────────────────────────────────────────────────────────
# 03_did_main.R   — Staggered DiD estimates (Callaway-Sant'Anna DR)
#
# Estimator : did::att_gt, est_method = "dr"
# Bootstrap : multiplier, B = 2000, clustered on AGS5 (county)
# ES bands  : simultaneous (cband = TRUE in aggte)
# Outcome   : bev_neuzulassungen_p100k (level; winsorised at 99th pct)
# Covariates (conditional spec, XFORMLA_CS from _dict.R):
#   sk_base_z + state_green_base_z + dens_base_z
#   (Tax capacity | State Green share | Log pop. density; baseline 2014-16)
#
# Output sections
# ────────────────
# A. Direct frame  — es_main_direct.{png,tex,csv}   (uncond | cond)
# B. Broad frame   — es_main_broad.{png,tex,csv}    (uncond | cond)
# C. Robustness: not-yet-treated (direct, ever-treated pool only)
#                  — es_robust_notyet.{png,tex,csv}  (uncond | cond)
# D. Corporate BEV — es_corp.{png,tex,csv}
# E. Private BEV   — es_priv.{png,tex,csv}
# F. ICE placebo   — es_ice.{png,tex,csv}
# G. ATT summary   — est_att_main.csv
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

# ── Constants ─────────────────────────────────────────────────────────────────

CS_BITERS    <- 2000L
OUTCOME_BEV  <- "bev_neuzulassungen_p100k"
OUTCOME_CORP <- "bev_corporate_p100k"
OUTCOME_PRIV <- "bev_private_p100k"
OUTCOME_ICE  <- "ice_neuzulassungen_p100k"

# Column keys and their display labels (graph facet strips + table headers)
COL_LABELS <- c(
  uncond  = "Unconditional",
  cond    = "Conditional",
  nyt_unc = "Not-yet-treated (uncond.)",
  nyt_con = "Not-yet-treated (cond.)",
  corp    = "Corporate BEV",
  priv    = "Private BEV",
  ice     = "ICE (placebo)"
)

# ── Plot theme (matches dual-map colour scheme) ────────────────────────────────

PLOT_BLUE <- "#004CFF"
PLOT_FILL <- "#9ecae1"

.font_family <- "sans"
.sf_pro <- "/Library/Fonts/SF-Pro-Display-Medium.otf"
if (requireNamespace("systemfonts", quietly = TRUE) &&
    file.exists(.sf_pro)) {
  systemfonts::register_font(name = "SFProDisplayMedium",
                             plain = .sf_pro)
  .font_family <- "SFProDisplayMedium"
}

theme_es <- theme_minimal(base_size = 11,
                          base_family = .font_family) +
  theme(
    strip.text       = element_text(face = "bold", size = 10),
    axis.line        = element_line(colour = "grey70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    plot.caption     = element_text(size = 7, colour = "grey50",
                                    hjust = 0),
    plot.title       = element_text(size = 11, face = "bold")
  )

# ── Load frames ───────────────────────────────────────────────────────────────

.read_frame <- function(p) fread(p, colClasses = list(
  character = c("AGS8", "AGS5", "AGS2", "treat_type")))

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

for (nm in names(frames)) {
  d            <- frames[[nm]]
  never_u      <- uniqueN(d$ags8_id[d$gname_cs == 0L])
  untreated_n  <- sum(d$gname_cs == 0L)
  cat(sprintf("CS [%s]: never_units=%d; untreated_obs=%d (%.0f%%)\n",
              nm, never_u, untreated_n,
              100 * untreated_n / nrow(d)))
  stopifnot(untreated_n > 0.5 * nrow(d), never_u > 1000L)
}

# ── Core CS estimator ─────────────────────────────────────────────────────────

run_cs <- function(yname, dat,
                   xformla      = NULL,
                   control_group = "nevertreated") {
  att_gt(
    yname                  = yname,
    gname                  = "gname_cs",
    idname                 = "ags8_id",
    tname                  = "year",
    data                   = as.data.frame(dat),
    control_group          = control_group,
    anticipation           = 0L,
    xformla                = xformla,
    est_method             = "dr",
    clustervars            = "AGS5",
    bstrap                 = TRUE,
    biters                 = CS_BITERS,
    allow_unbalanced_panel = TRUE,
    print_details          = FALSE
  )
}

# ── Aggregation helpers ────────────────────────────────────────────────────────

cs_es_agg <- function(cs_obj, es_max) {
  if (is.null(cs_obj)) return(NULL)
  res <- tryCatch(
    aggte(cs_obj, type = "dynamic",
          min_e = ES_MIN, max_e = es_max,
          balance_e = NULL, na.rm = TRUE, cband = TRUE),
    error = function(e) {
      cat("  aggte(dynamic, cband=TRUE) err:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(res)) {
    cat("  aggte: retrying with cband=FALSE\n")
    res <- tryCatch(
      aggte(cs_obj, type = "dynamic",
            min_e = ES_MIN, max_e = es_max,
            balance_e = NULL, na.rm = TRUE, cband = FALSE),
      error = function(e) {
        cat("  aggte(dynamic, cband=FALSE) err:", conditionMessage(e), "\n")
        NULL
      }
    )
    # When cband=FALSE, crit.val.egt is NULL; replace with 1.96 so CI formulas
    # in es_graph/write_es_longtblr degrade gracefully to pointwise 95% bands.
    if (!is.null(res) && is.null(res$crit.val.egt))
      res$crit.val.egt <- 1.96
  }
  res
}

cs_att_agg <- function(cs_obj) {
  if (is.null(cs_obj)) return(NULL)
  tryCatch(aggte(cs_obj, type = "simple", na.rm = TRUE),
           error = function(e) NULL)
}

# Pre-treatment Wald test: joint H0 that all pre-period ATTs == 0.
# Uses the influence-function covariance matrix from the ES aggte object.
# Falls back to a diagonal chi-squared if the matrix is unavailable.
cs_pre_test <- function(es_agg) {
  if (is.null(es_agg))
    return(list(stat = NA_real_, df = NA_integer_, pval = NA_real_))
  pre_idx <- which(es_agg$egt < 0)
  if (length(pre_idx) == 0)
    return(list(stat = NA_real_, df = 0L, pval = NA_real_))

  att_pre <- es_agg$att.egt[pre_idx]
  inf     <- es_agg$inf.func.egt

  W <- if (!is.null(inf) && ncol(inf) >= max(pre_idx)) {
    inf_pre <- inf[, pre_idx, drop = FALSE]
    n       <- nrow(inf_pre)
    V       <- (t(inf_pre) %*% inf_pre) / n^2
    tryCatch(
      as.numeric(t(att_pre) %*% solve(V) %*% att_pre),
      error = function(e)
        sum((att_pre / es_agg$se.egt[pre_idx])^2, na.rm = TRUE)
    )
  } else {
    sum((att_pre / es_agg$se.egt[pre_idx])^2, na.rm = TRUE)
  }
  df <- length(pre_idx)
  list(stat = W, df = df, pval = pchisq(W, df = df, lower.tail = FALSE))
}

# ── Single-cell runner ────────────────────────────────────────────────────────

run_spec <- function(dat, yname, xformla, control_group, label) {
  dat_y <- dat[!is.na(get(yname))]

  # For conditional specs: drop AGS8s with NA in any xformla covariate.
  # Baselines are time-invariant, so NA in one year = NA in all years.
  if (!is.null(xformla)) {
    cov_nms <- all.vars(xformla)
    missing_cols <- setdiff(cov_nms, names(dat_y))
    if (length(missing_cols)) {
      cat(sprintf("  [%s] MISSING covariate columns: %s — skip\n",
                  label, paste(missing_cols, collapse = ", ")))
      return(NULL)
    }
    # One non-NA check per AGS8 (baselines are constant within unit)
    cc <- dat_y[, lapply(.SD, function(x) any(!is.na(x))),
                by = AGS8, .SDcols = cov_nms]
    # .SDcols in the outer [ makes .SD refer to the covariate columns only
    ok_ags8 <- cc[rowSums(as.matrix(.SD)) == length(cov_nms), AGS8,
                  .SDcols = cov_nms]
    n_drop  <- uniqueN(dat_y$AGS8) - length(ok_ags8)
    if (n_drop > 0L) {
      cat(sprintf("  [%s] dropping %d AGS8 with NA in xformla covariate(s)\n",
                  label, n_drop))
      dat_y <- dat_y[AGS8 %in% ok_ags8]
    }
  }

  n_tr  <- uniqueN(dat_y[gname_cs > 0L, AGS8])
  n_co  <- uniqueN(dat_y[gname_cs == 0L, AGS8])
  if (n_tr < 5L) {
    cat(sprintf("  [%s] n_treated=%d; skip\n", label, n_tr))
    return(NULL)
  }
  es_max <- es_max_data_driven(dat_y, yname)
  cat(sprintf("  [%s] n_tr=%d  n_co=%d  es_max=%d\n",
              label, n_tr, n_co, es_max))

  cs_obj <- tryCatch(
    withCallingHandlers(
      run_cs(yname, dat_y, xformla, control_group),
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
    n_na  <- sum(is.na(cs_obj$att))
    n_tot <- length(cs_obj$att)
    if (n_na > 0L)
      cat(sprintf("  [%s] %d / %d (g,t) ATTs are NA\n", label, n_na, n_tot))
  }

  es_agg  <- cs_es_agg(cs_obj, es_max)
  att_agg <- cs_att_agg(cs_obj)
  pre_tst <- cs_pre_test(es_agg)

  list(es_agg    = es_agg,
       att_agg   = att_agg,
       pre_tst   = pre_tst,
       n_treated = n_tr,
       n_control = n_co)
}

# ── ES graph ──────────────────────────────────────────────────────────────────
# Uses simultaneous CI bands (crit.val.egt from aggte with cband = TRUE).

es_graph <- function(es_list, file,
                     ylab  = "ATT (BEV new registrations per 100k)",
                     title = NULL) {
  dt_list <- lapply(names(es_list), function(nm) {
    ea <- es_list[[nm]]
    if (is.null(ea)) return(NULL)
    cv <- ea$crit.val.egt
    data.table(
      spec     = factor(COL_LABELS[nm], levels = COL_LABELS),
      e        = ea$egt,
      estimate = ea$att.egt,
      ci_lo    = ea$att.egt - cv * ea$se.egt,
      ci_hi    = ea$att.egt + cv * ea$se.egt
    )
  })
  dt <- rbindlist(Filter(Negate(is.null), dt_list), use.names = TRUE)
  if (nrow(dt) == 0L) {
    cat("  es_graph: no data for", file, "\n"); return(invisible(NULL))
  }

  n_panels <- uniqueN(dt$spec)
  cap <- paste0(
    "Simultaneous 95% CI (multiplier bootstrap, B=", CS_BITERS,
    ", AGS5 cluster). CS-DR estimator."
  )

  p <- ggplot(dt, aes(x = e)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.35) +
    geom_vline(xintercept = -0.5, linetype = "dashed",
               colour = "grey60", linewidth = 0.35) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                fill = PLOT_FILL, alpha = 0.40, colour = NA) +
    geom_line(aes(y = estimate),
              colour = PLOT_BLUE, linewidth = 0.75) +
    geom_point(aes(y = estimate),
               colour = PLOT_BLUE, size = 1.8) +
    facet_wrap(~ spec, ncol = n_panels, scales = "free_y") +
    labs(x     = "Years relative to first funding receipt",
         y     = ylab,
         title = title,
         caption = cap) +
    theme_es +
    theme(plot.background = element_rect(fill = "white", colour = NA))
  ggsave(file, p,
         width  = 3.5 * n_panels + 1.5,
         height = 4.5,
         dpi    = 300,
         device = ragg::agg_png)
  invisible(NULL)
}

# ── ES table (longtblr) ───────────────────────────────────────────────────────

.st <- function(b, s) {
  p <- 2 * pnorm(abs(b / s), lower.tail = FALSE)
  ifelse(is.na(p), "",
    ifelse(p < 0.01, "***",
      ifelse(p < 0.05, "**",
        ifelse(p < 0.1, "*", ""))))
}
.f3 <- function(x) ifelse(is.na(x) | !is.finite(x), "--",
                           sprintf("%.3f", x))
.fn <- function(x) {
  if (is.null(x) || is.na(x) || x == 0L) "--"
  else format(as.integer(x), big.mark = ",", scientific = FALSE)
}

write_es_longtblr <- function(es_list, att_list, pre_tests,
                               n_treated, n_control,
                               stem, caption, label, note) {
  col_nms <- names(es_list)
  hdrs    <- COL_LABELS[col_nms]
  n_cols  <- length(col_nms)

  all_e <- sort(unique(unlist(
    lapply(es_list, function(ea) if (!is.null(ea)) ea$egt)
  )))

  # Header
  hdr_row <- paste(
    c("$e$", hdrs), collapse = " & "
  )
  hdr_row <- paste0(hdr_row, " \\\\")

  # Body rows
  body <- character(0)
  for (ev in all_e) {
    r_est <- sprintf("%+d", as.integer(ev))
    r_se  <- ""
    for (nm in col_nms) {
      ea <- es_list[[nm]]
      if (is.null(ea)) {
        r_est <- paste0(r_est, " & ")
        r_se  <- paste0(r_se,  " & ")
        next
      }
      idx <- which(ea$egt == ev)
      if (length(idx) == 0L) {
        r_est <- paste0(r_est, " & ")
        r_se  <- paste0(r_se,  " & ")
        next
      }
      b <- ea$att.egt[idx]; s <- ea$se.egt[idx]
      r_est <- paste0(r_est, " & ", sprintf("%.3f%s", b, .st(b, s)))
      r_se  <- paste0(r_se,  " & ", sprintf("(%.3f)", s))
    }
    body <- c(body,
              paste0(r_est, " \\\\"),
              paste0(r_se,  " \\\\"))
    if (ev == -1L) body <- c(body, "\\hline")
  }

  # Footer
  foot <- "\\hline"

  # Overall ATT
  att_e <- vapply(col_nms, function(nm) {
    a <- att_list[[nm]]
    if (is.null(a)) return("--")
    sprintf("%.3f%s", a$overall.att, .st(a$overall.att, a$overall.se))
  }, character(1))
  att_s <- vapply(col_nms, function(nm) {
    a <- att_list[[nm]]
    if (is.null(a)) return("") else sprintf("(%.3f)", a$overall.se)
  }, character(1))
  foot <- c(foot,
    paste(c("Overall ATT", att_e), collapse = " & "), " \\\\",
    paste(c("",            att_s), collapse = " & "), " \\\\",
    paste(c("$N_{\\text{treated}}$",
            vapply(col_nms, function(nm) .fn(n_treated[[nm]]),
                   character(1))),
          collapse = " & "), " \\\\",
    paste(c("$N_{\\text{control}}$",
            vapply(col_nms, function(nm) .fn(n_control[[nm]]),
                   character(1))),
          collapse = " & "), " \\\\",
    "\\hline",
    paste(c("Pre-test $\\chi^2$",
            vapply(col_nms, function(nm)
              .f3(pre_tests[[nm]]$stat), character(1))),
          collapse = " & "), " \\\\",
    paste(c("\\quad $df$",
            vapply(col_nms, function(nm) {
              d <- pre_tests[[nm]]$df
              if (is.null(d) || is.na(d)) "--" else as.character(d)
            }, character(1))),
          collapse = " & "), " \\\\",
    paste(c("\\quad $p$-value",
            vapply(col_nms, function(nm)
              .f3(pre_tests[[nm]]$pval), character(1))),
          collapse = " & "), " \\\\"
  )

  write_longtblr(
    stem        = stem,
    caption     = caption,
    label       = label,
    note        = note,
    colspec     = sprintf("r *{%d}{r}", n_cols),
    header_rows = hdr_row,
    body_rows   = body,
    footer_rows = foot
  )

  # CSV twin
  csv_rows <- lapply(col_nms, function(nm) {
    ea  <- es_list[[nm]]
    att <- att_list[[nm]]
    pt  <- pre_tests[[nm]]
    if (is.null(ea)) return(NULL)
    cv <- ea$crit.val.egt
    data.table(
      spec            = nm,
      e               = ea$egt,
      estimate        = ea$att.egt,
      se              = ea$se.egt,
      crit_val        = cv,
      ci_lo_simult    = ea$att.egt - cv * ea$se.egt,
      ci_hi_simult    = ea$att.egt + cv * ea$se.egt,
      ci_lo_pw        = ea$att.egt - 1.96 * ea$se.egt,
      ci_hi_pw        = ea$att.egt + 1.96 * ea$se.egt,
      att_overall     = if (!is.null(att)) att$overall.att else NA_real_,
      att_se          = if (!is.null(att)) att$overall.se  else NA_real_,
      pre_test_stat   = pt$stat,
      pre_test_df     = pt$df,
      pre_test_pval   = pt$pval,
      n_treated       = n_treated[[nm]],
      n_control       = n_control[[nm]]
    )
  })
  fwrite(
    rbindlist(Filter(Negate(is.null), csv_rows),
              use.names = TRUE, fill = TRUE),
    paste0(stem, ".csv")
  )
  invisible(NULL)
}

# ── Table note strings ─────────────────────────────────────────────────────────

NOTE_MAIN <- paste0(
  "Callaway--Sant'Anna (2021) CS-DR estimator. ",
  "Outcome: BEV new registrations per 100,000 population ",
  "(level; winsorised at 99th percentile). ",
  "Bootstrap SEs clustered at county (AGS5), $B=", CS_BITERS, "$. ",
  "Graph CI bands: simultaneous (95\\%) where available; ",
  "pointwise (95\\%) as fallback for conditional spec. ",
  "Table: pointwise SEs; ",
  "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$ (two-sided). ",
  "Conditional spec covariates (z-scored, earliest year 2014--16): ",
  "tax capacity, state Green vote share, log pop.\\ density. ",
  "Pre-treatment test: Wald $\\chi^2$ on joint pre-period ATTs ",
  "(covariance via bootstrap influence functions)."
)

NOTE_NOTYET <- paste0(
  NOTE_MAIN,
  " Control group: not-yet-treated cohorts (ever-treated units only)."
)

note_sec <- function(outcome_lbl)
  paste0(
    "Specification as per the unconditional direct-frame spec ",
    "(nevertreated control, no covariates), but outcome: ",
    outcome_lbl, "."
  )

# ── ATT summary helper ────────────────────────────────────────────────────────

.att_row <- function(nm, frame, spec_lbl, ctrl, outcome, res) {
  if (is.null(res) || is.null(res$att_agg)) return(NULL)
  a <- res$att_agg
  data.table(
    spec = nm, frame = frame, xformla = spec_lbl, control = ctrl,
    outcome  = outcome,
    att      = a$overall.att,
    se       = a$overall.se,
    ci_lo    = a$overall.att - 1.96 * a$overall.se,
    ci_hi    = a$overall.att + 1.96 * a$overall.se,
    n_treated = res$n_treated,
    n_control = res$n_control
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
# A + B  Main specs — direct and broad frames
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n=== A. Direct frame ===\n")
r_dir_unc <- run_spec(frames$direct, OUTCOME_BEV, NULL,
                      "nevertreated", "direct/uncond")
r_dir_con <- run_spec(frames$direct, OUTCOME_BEV, XFORMLA_CS,
                      "nevertreated", "direct/cond")

cat("\n=== B. Broad frame ===\n")
r_brd_unc <- run_spec(frames$broad, OUTCOME_BEV, NULL,
                      "nevertreated", "broad/uncond")
r_brd_con <- run_spec(frames$broad, OUTCOME_BEV, XFORMLA_CS,
                      "nevertreated", "broad/cond")

# ═══════════════════════════════════════════════════════════════════════════════
# C  Robustness — not-yet-treated (direct frame, ever-treated pool only)
# ═══════════════════════════════════════════════════════════════════════════════

# Filter to ever-treated only: never-treated units (gname_cs == 0) are not
# valid controls under notyettreated and would simply be ignored by att_gt,
# so we drop them upfront for clarity.
dat_notyet <- frames$direct[gname_cs > 0L]

cat("\n=== C. Robustness: not-yet-treated (direct frame) ===\n")
r_nyt_unc <- run_spec(dat_notyet, OUTCOME_BEV, NULL,
                      "notyettreated", "notyet/uncond")
r_nyt_con <- run_spec(dat_notyet, OUTCOME_BEV, XFORMLA_CS,
                      "notyettreated", "notyet/cond")

# ═══════════════════════════════════════════════════════════════════════════════
# D – F  Secondary outcomes (direct, nevertreated, unconditional)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n=== D-F. Secondary outcomes ===\n")
r_corp <- run_spec(frames$direct, OUTCOME_CORP, NULL,
                   "nevertreated", "direct/corp")
r_priv <- run_spec(frames$direct, OUTCOME_PRIV, NULL,
                   "nevertreated", "direct/priv")
r_ice  <- run_spec(frames$direct, OUTCOME_ICE,  NULL,
                   "nevertreated", "direct/ice")

# ═══════════════════════════════════════════════════════════════════════════════
# Write outputs
# ═══════════════════════════════════════════════════════════════════════════════

# Helper: build named lists for one or two-column groups
.L <- function(...) {
  args <- list(...)
  setNames(args, names(args))
}

# ── A. Direct frame ────────────────────────────────────────────────────────────

cat("\nWriting A: direct frame\n")
es_graph(
  es_list = list(uncond = r_dir_unc$es_agg, cond = r_dir_con$es_agg),
  file    = file.path(out_dir, "es_main_direct.png"),
  title   = "BEV New Registrations — Direct Treatment Frame"
)
write_es_longtblr(
  es_list   = list(uncond = r_dir_unc$es_agg,  cond = r_dir_con$es_agg),
  att_list  = list(uncond = r_dir_unc$att_agg, cond = r_dir_con$att_agg),
  pre_tests = list(uncond = r_dir_unc$pre_tst, cond = r_dir_con$pre_tst),
  n_treated = list(uncond = r_dir_unc$n_treated,
                   cond   = r_dir_con$n_treated),
  n_control = list(uncond = r_dir_unc$n_control,
                   cond   = r_dir_con$n_control),
  stem    = file.path(out_dir, "es_main_direct"),
  caption = paste0("Event Study: BEV New Registrations per 100k",
                   " --- Direct Treatment"),
  label   = "tab:es_main_direct",
  note    = NOTE_MAIN
)

# ── B. Broad frame ─────────────────────────────────────────────────────────────

cat("Writing B: broad frame\n")
es_graph(
  es_list = list(uncond = r_brd_unc$es_agg, cond = r_brd_con$es_agg),
  file    = file.path(out_dir, "es_main_broad.png"),
  title   = "BEV New Registrations --- Broad Treatment Frame"
)
write_es_longtblr(
  es_list   = list(uncond = r_brd_unc$es_agg,  cond = r_brd_con$es_agg),
  att_list  = list(uncond = r_brd_unc$att_agg, cond = r_brd_con$att_agg),
  pre_tests = list(uncond = r_brd_unc$pre_tst, cond = r_brd_con$pre_tst),
  n_treated = list(uncond = r_brd_unc$n_treated,
                   cond   = r_brd_con$n_treated),
  n_control = list(uncond = r_brd_unc$n_control,
                   cond   = r_brd_con$n_control),
  stem    = file.path(out_dir, "es_main_broad"),
  caption = paste0("Event Study: BEV New Registrations per 100k",
                   " --- Broad Treatment"),
  label   = "tab:es_main_broad",
  note    = NOTE_MAIN
)

# ── C. Robustness: not-yet-treated ─────────────────────────────────────────────

cat("Writing C: robustness not-yet-treated\n")
es_graph(
  es_list = list(nyt_unc = r_nyt_unc$es_agg,
                 nyt_con = r_nyt_con$es_agg),
  file    = file.path(out_dir, "es_robust_notyet.png"),
  title   = "Robustness: Not-yet-treated Control (Direct Frame)"
)
write_es_longtblr(
  es_list   = list(nyt_unc = r_nyt_unc$es_agg,
                   nyt_con = r_nyt_con$es_agg),
  att_list  = list(nyt_unc = r_nyt_unc$att_agg,
                   nyt_con = r_nyt_con$att_agg),
  pre_tests = list(nyt_unc = r_nyt_unc$pre_tst,
                   nyt_con = r_nyt_con$pre_tst),
  n_treated = list(nyt_unc = r_nyt_unc$n_treated,
                   nyt_con = r_nyt_con$n_treated),
  n_control = list(nyt_unc = r_nyt_unc$n_control,
                   nyt_con = r_nyt_con$n_control),
  stem    = file.path(out_dir, "es_robust_notyet"),
  caption = paste0("Robustness: Not-yet-treated Control Group",
                   " (Direct Frame)"),
  label   = "tab:es_robust_notyet",
  note    = NOTE_NOTYET
)

# ── D-F. Secondary outcomes ────────────────────────────────────────────────────

sec_specs <- list(
  list(
    res   = r_corp, key = "corp",
    file  = "es_corp",
    cap   = "Event Study: Corporate BEV Registrations per 100k",
    lbl   = "tab:es_corp",
    ylab  = "ATT (corporate BEV registrations per 100k)",
    olbl  = "corporate BEV new registrations per 100k"
  ),
  list(
    res   = r_priv, key = "priv",
    file  = "es_priv",
    cap   = "Event Study: Private BEV Registrations per 100k",
    lbl   = "tab:es_priv",
    ylab  = "ATT (private BEV registrations per 100k)",
    olbl  = "private BEV new registrations per 100k"
  ),
  list(
    res   = r_ice, key = "ice",
    file  = "es_ice",
    cap   = "Placebo: ICE New Registrations per 100k (Direct Frame)",
    lbl   = "tab:es_ice",
    ylab  = "ATT (ICE new registrations per 100k)",
    olbl  = "ICE new registrations per 100k (placebo)"
  )
)

for (sp in sec_specs) {
  key <- sp$key; res <- sp$res
  if (is.null(res)) {
    cat(sprintf("  skip %s (estimation returned NULL)\n", key)); next
  }
  cat(sprintf("Writing %s\n", sp$file))
  es_graph(
    es_list = setNames(list(res$es_agg), key),
    file    = file.path(out_dir, paste0(sp$file, ".png")),
    ylab    = sp$ylab,
    title   = sp$cap
  )
  write_es_longtblr(
    es_list   = setNames(list(res$es_agg),  key),
    att_list  = setNames(list(res$att_agg), key),
    pre_tests = setNames(list(res$pre_tst), key),
    n_treated = setNames(list(res$n_treated), key),
    n_control = setNames(list(res$n_control), key),
    stem    = file.path(out_dir, sp$file),
    caption = sp$cap,
    label   = sp$lbl,
    note    = note_sec(sp$olbl)
  )
}

# ── G. ATT summary CSV ────────────────────────────────────────────────────────

att_summary <- rbindlist(list(
  .att_row("dir_unc", "direct", "unconditional", "nevertreated",
           OUTCOME_BEV, r_dir_unc),
  .att_row("dir_con", "direct", "conditional",   "nevertreated",
           OUTCOME_BEV, r_dir_con),
  .att_row("brd_unc", "broad",  "unconditional", "nevertreated",
           OUTCOME_BEV, r_brd_unc),
  .att_row("brd_con", "broad",  "conditional",   "nevertreated",
           OUTCOME_BEV, r_brd_con),
  .att_row("nyt_unc", "direct", "unconditional", "notyettreated",
           OUTCOME_BEV, r_nyt_unc),
  .att_row("nyt_con", "direct", "conditional",   "notyettreated",
           OUTCOME_BEV, r_nyt_con),
  .att_row("corp",    "direct", "unconditional", "nevertreated",
           OUTCOME_CORP, r_corp),
  .att_row("priv",    "direct", "unconditional", "nevertreated",
           OUTCOME_PRIV, r_priv),
  .att_row("ice",     "direct", "unconditional", "nevertreated",
           OUTCOME_ICE,  r_ice)
), use.names = TRUE, fill = TRUE)

fwrite(att_summary, file.path(out_dir, "est_att_main.csv"))
cat(sprintf("\nAll DiD main outputs -> %s\n", out_dir))
