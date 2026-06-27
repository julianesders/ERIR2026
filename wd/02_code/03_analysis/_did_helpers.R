# ─────────────────────────────────────────────────────────────────────────────
# _did_helpers.R — shared Callaway-Sant'Anna estimation + output machinery.
#
# Sourced by 04_did_main.R AND 07_spillovers.R (after _dict.R and the
# data.table / did / ggplot2 libraries are loaded). Contains ONLY definitions —
# no top-level estimation — so every script that uses CS produces byte-identical
# event-study figures (es_graph), longtblr tables + CSV twins (write_es_longtblr)
# and ATT summaries (.att_row). The ONLY difference between callers is the
# estimation SAMPLE handed to run_spec().
#
# Depends on _dict.R: ES_MIN, es_max_data_driven(), XFORMLA_CS, write_longtblr().
# ─────────────────────────────────────────────────────────────────────────────

library(data.table)
library(did)
library(ggplot2)

# ── Constants ─────────────────────────────────────────────────────────────────

CS_BITERS   <- 2000L
OUTCOME_BEV <- "bev_neuzulassungen_p100k"

# Column keys and their display labels (graph facet strips + table headers).
# Conditional CS (doubly robust) is the PRIMARY spec; unconditional is the
# comparison column. Order here = left-to-right display order, so conditional
# prints first.
COL_LABELS <- c(
  cond_dr     = "Conditional",
  uncond      = "Unconditional",
  nyt_dr      = "Not-yet-treated (conditional)",
  nyt_uncond  = "Not-yet-treated (unconditional)",
  donut       = "Donut (direct frame)",
  spillover   = "Spillover",
  # Combined-table keys (direct vs broad, used in es_main_combined)
  dir_cond_dr = "Direct: Conditional",
  dir_uncond  = "Direct: Unconditional",
  brd_cond_dr = "Broad: Conditional",
  brd_uncond  = "Broad: Unconditional"
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

# ── Frame reader ──────────────────────────────────────────────────────────────

.read_frame <- function(p) fread(p, colClasses = list(
  character = c("AGS8", "AGS5", "AGS2", "treat_type")))

# ── Core CS estimator ─────────────────────────────────────────────────────────
# est_method threaded through (NOT boot_type — that argument does not exist in
# `did` and would crash). Unconditional specs (xformla = NULL) are estimator-
# invariant, so passing "reg" vs "dr" there is immaterial.

run_cs <- function(yname, dat,
                   xformla       = NULL,
                   control_group = "nevertreated",
                   est_method    = "dr",
                   gname         = "gname_cs") {
  att_gt(
    yname                  = yname,
    gname                  = gname,
    idname                 = "ags8_id",
    tname                  = "year",
    data                   = as.data.frame(dat),
    control_group          = control_group,
    anticipation           = 0L,
    xformla                = xformla,
    est_method             = est_method,
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

# Joint pre-trends Wald. `did::att_gt` computes a correctly AGS5-clustered Wald
# test of H0: all pre-treatment ATT(g,t) == 0 (multiplier bootstrap) and stores
# it on the MP object as `cs_obj$W` (statistic) and `cs_obj$Wpval` (p-value) —
# this is the "P-value for pre-test of parallel trends assumption" the package
# prints, with Wpval == pchisq(W, df, lower.tail = FALSE) and df = number of
# pre-treatment (g,t) cells.
cs_pre_test <- function(cs_obj) {
  na_out <- list(stat = NA_real_, df = NA_integer_, pval = NA_real_,
                 method = "unavailable")
  if (is.null(cs_obj)) return(na_out)
  n_pre <- sum(cs_obj$t < cs_obj$group)
  if (!is.finite(n_pre) || n_pre == 0L)
    return(list(stat = NA_real_, df = 0L, pval = NA_real_,
                method = "no_pre_periods"))
  W  <- cs_obj$W
  pv <- cs_obj$Wpval
  if (is.null(W) || is.null(pv) || !is.finite(W) || !is.finite(pv))
    return(na_out)
  list(stat = as.numeric(W), df = as.integer(n_pre),
       pval = as.numeric(pv), method = "joint_wald")
}

# ── Single-cell runner ────────────────────────────────────────────────────────
# `gname` defaults to "gname_cs"; pass a different cohort column (e.g. a
# pseudo-treatment) to reuse the identical machinery on a different design.

run_spec <- function(dat, yname, xformla, control_group, est_method, label,
                     gname = "gname_cs") {
  dat_y <- dat[!is.na(get(yname))]

  # Conditional covariate handling — CORRECTED IDIOM.
  # Baselines are time-invariant within AGS8, so NA in any year => NA in all.
  # Compute a per-AGS8 "all covariates non-NA" flag in ONE pass, then subset
  # AGS8 from THAT result. Do NOT fold the completeness reduction and the AGS8
  # selection into a single [ ] with .SDcols — that was the original bug.
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

  n_tr <- uniqueN(dat_y[get(gname) > 0L, AGS8])
  n_co <- uniqueN(dat_y[get(gname) == 0L, AGS8])
  if (n_tr < 5L) {
    cat(sprintf("  [%s] n_treated=%d; skip\n", label, n_tr))
    return(NULL)
  }

  # dr overlap diagnostic: print treated-N per surviving cohort; warn if <10
  # (the within-cohort propensity score may be unstable below that).
  if (est_method == "dr" && !is.null(xformla)) {
    cn <- dat_y[get(gname) > 0L, .(n = uniqueN(AGS8)),
                by = c(gname)][order(get(gname))]
    cat(sprintf("  [%s] dr cohort sizes: %s\n", label,
                paste(sprintf("%d(%d)", cn[[gname]], cn$n), collapse = " ")))
    thin <- cn[n < 10L]
    if (nrow(thin))
      cat(sprintf(paste0("  [%s] WARNING: %d cohort(s) <10 treated under dr ",
                         "(propensity may be unstable): %s\n"),
                  label, nrow(thin), paste(thin[[gname]], collapse = ", ")))
  }

  es_max <- es_max_data_driven(dat_y, yname, gname_col = gname)
  cat(sprintf("  [%s] n_tr=%d  n_co=%d  es_max=%d  est=%s\n",
              label, n_tr, n_co, es_max, est_method))

  cs_obj <- tryCatch(
    withCallingHandlers(
      run_cs(yname, dat_y, xformla, control_group, est_method, gname = gname),
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
  pre_tst <- cs_pre_test(cs_obj)

  list(es_agg     = es_agg,
       att_agg    = att_agg,
       pre_tst    = pre_tst,
       n_treated  = n_tr,
       n_control  = n_co,
       est_method = est_method)
}

# Default CS event-study figure note (CS-specific; write_fig_tex itself lives in
# _dict.R as a generic float writer). emit_section uses this unless overridden.
FIG_NOTE_DEFAULT <- paste0(
  "Callaway--Sant'Anna (2021) CS event-study estimates with simultaneous ",
  "95\\% confidence bands (multiplier bootstrap, $B=", CS_BITERS, "$, ",
  "clustered at county/AGS5); pointwise 95\\% fallback where the simultaneous ",
  "critical value is unavailable. The dashed vertical line marks the last ",
  "pre-treatment period; $e$ is years relative to first funding receipt."
)

# ── ES graph ──────────────────────────────────────────────────────────────────
# Shows ALL estimable columns for a section side-by-side (the main/appendix
# split applies to the *tables*, not the figure). Uses simultaneous CI bands
# (crit.val.egt from aggte with cband = TRUE; pointwise 1.96 fallback). The plot
# carries no title/caption — those go in the companion .tex (write_fig_tex).
# `width`/`height` (inches) override the per-panel defaults; pass them for the
# combined 2x2 main figure so it lays out as a wide A4-friendly rectangle.

es_graph <- function(es_list, file,
                     ylab   = "ATT (BEV new registrations per 100k)",
                     ncol   = NULL,
                     width  = NULL,
                     height = NULL) {
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
  # Re-enforce COL_LABELS ordering after rbindlist (which may coerce factors)
  present <- unique(as.character(dt$spec))
  ordered_specs <- as.character(COL_LABELS[COL_LABELS %in% present])
  dt[, spec := factor(as.character(spec), levels = ordered_specs)]

  n_panels  <- uniqueN(dt$spec)
  ncol_use  <- if (!is.null(ncol)) ncol else n_panels
  nrow_use  <- ceiling(n_panels / ncol_use)

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
    facet_wrap(~ spec, ncol = ncol_use, scales = "free_y") +
    labs(x = "Years relative to first funding receipt",
         y = ylab) +
    theme_es +
    theme(plot.background = element_rect(fill = "white", colour = NA))
  ggsave(file, p,
         width  = if (!is.null(width))  width  else 3.5 * ncol_use + 1.5,
         height = if (!is.null(height)) height else 4.5 * nrow_use,
         dpi    = 300,
         device = ragg::agg_png)
  invisible(NULL)
}

# ── ES table (longtblr) ───────────────────────────────────────────────────────

.st <- function(b, s) {
  p <- 2 * pnorm(abs(b / s), lower.tail = FALSE)
  ifelse(is.na(p), "",
    ifelse(p < 0.01, "$^{***}$",
      ifelse(p < 0.05, "$^{**}$",
        ifelse(p < 0.1, "$^{*}$", ""))))
}
.f3 <- function(x) ifelse(is.na(x) | !is.finite(x), "--",
                           sprintf("%.3f", x))
.fn <- function(x) {
  if (is.null(x) || is.na(x) || x == 0L) "--"
  else format(as.integer(x), big.mark = ",", scientific = FALSE)
}

# `display_cols` are the columns rendered in the TeX table; `est_methods` is a
# named char vector (column key -> "reg"/"dr"). The CSV twin records EVERY
# estimated column (the union), so no estimate is lost to a main/appendix split.
write_es_longtblr <- function(es_list, att_list, pre_tests,
                              n_treated, n_control, est_methods,
                              stem, caption, label, note,
                              display_cols = names(es_list),
                              resize = FALSE) {
  col_nms  <- names(es_list)          # all estimated columns -> CSV (union)
  disp_nms <- display_cols            # displayed columns      -> TeX
  hdrs     <- COL_LABELS[disp_nms]
  n_cols   <- length(disp_nms)

  all_e <- sort(unique(unlist(
    lapply(disp_nms, function(nm) {
      ea <- es_list[[nm]]; if (!is.null(ea)) ea$egt
    })
  )))

  # Header
  hdr_row <- paste0(paste(c("Relative period", hdrs), collapse = " & "), " \\\\")

  # Body rows
  body <- character(0)
  for (ev in all_e) {
    r_est <- sprintf("%+d", as.integer(ev))
    r_se  <- ""
    for (nm in disp_nms) {
      ea  <- es_list[[nm]]
      idx <- if (is.null(ea)) integer(0) else which(ea$egt == ev)
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

  # Footer: overall ATT, Ns, joint-Wald pre-test.
  att_e <- vapply(disp_nms, function(nm) {
    a <- att_list[[nm]]
    if (is.null(a)) return("--")
    sprintf("%.3f%s", a$overall.att, .st(a$overall.att, a$overall.se))
  }, character(1))
  att_s <- vapply(disp_nms, function(nm) {
    a <- att_list[[nm]]
    if (is.null(a)) return("") else sprintf("(%.3f)", a$overall.se)
  }, character(1))

  # Pre-test footer: joint Wald only; "--" for ALL fields when method is not
  # "joint_wald" (no diagonal fallback).
  .pt_field <- function(nm, field) {
    pt <- pre_tests[[nm]]
    if (is.null(pt) || is.null(pt$method) || pt$method != "joint_wald")
      return("--")
    if (field == "df") return(as.character(pt$df))
    .f3(pt[[field]])
  }

  foot <- c("\\hline",
    paste(c("Overall ATT", att_e), collapse = " & "), " \\\\",
    paste(c("",            att_s), collapse = " & "), " \\\\",
    paste(c("$N_{\\text{treated}}$",
            vapply(disp_nms, function(nm) .fn(n_treated[[nm]]),
                   character(1))),
          collapse = " & "), " \\\\",
    paste(c("$N_{\\text{control}}$",
            vapply(disp_nms, function(nm) .fn(n_control[[nm]]),
                   character(1))),
          collapse = " & "), " \\\\",
    "\\hline",
    paste(c("Pre-test $p$-value",
            vapply(disp_nms, .pt_field, character(1), field = "pval")),
          collapse = " & "), " \\\\"
  )

  write_longtblr(
    stem        = stem,
    caption     = caption,
    label       = label,
    note        = note,
    colspec     = sprintf("l *{%d}{l}", n_cols),
    header_rows = hdr_row,
    body_rows   = body,
    footer_rows = foot,
    resize      = resize
  )

  # CSV twin — union of ALL estimated columns (incl. dr).
  csv_rows <- lapply(col_nms, function(nm) {
    ea  <- es_list[[nm]]
    att <- att_list[[nm]]
    pt  <- pre_tests[[nm]]
    if (is.null(ea)) return(NULL)
    cv <- ea$crit.val.egt
    data.table(
      spec            = nm,
      est_method      = est_methods[[nm]],
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
      pre_test_method = pt$method,
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

# ── ATT summary helper ────────────────────────────────────────────────────────

.att_row <- function(nm, frame, spec_lbl, est_method, ctrl, outcome, res) {
  if (is.null(res) || is.null(res$att_agg)) return(NULL)
  a <- res$att_agg
  data.table(
    spec       = nm, frame = frame, xformla = spec_lbl,
    est_method = est_method, control = ctrl, outcome = outcome,
    att        = a$overall.att,
    se         = a$overall.se,
    ci_lo      = a$overall.att - 1.96 * a$overall.se,
    ci_hi      = a$overall.att + 1.96 * a$overall.se,
    n_treated  = res$n_treated,
    n_control  = res$n_control
  )
}

# ── Section emitter ────────────────────────────────────────────────────────────
# `results` is a named list of run_spec() outputs keyed by COL_LABELS keys;
# `ems` maps the same keys to "reg"/"dr". The figure shows all estimable
# columns; the main TeX shows `main_display`; the optional dr appendix twin
# shows `dr_key` only. The main CSV is always the union of all `results`.

emit_section <- function(results, ems,
                         graph_file, graph_title,
                         main_stem, main_caption, main_label, main_note,
                         main_display = names(results),
                         dr_key = NULL, dr_stem = NULL,
                         dr_caption = NULL, dr_label = NULL, dr_note = NULL,
                         ylab = "ATT (BEV new registrations per 100k)",
                         plot = TRUE,
                         ncol = NULL,
                         fig_caption = NULL, fig_label = NULL,
                         fig_note = FIG_NOTE_DEFAULT,
                         fig_width = NULL, fig_height = NULL) {
  es  <- lapply(results, `[[`, "es_agg")
  att <- lapply(results, `[[`, "att_agg")
  pre <- lapply(results, `[[`, "pre_tst")
  ntr <- lapply(results, `[[`, "n_treated")
  nco <- lapply(results, `[[`, "n_control")

  if (plot) {
    es_graph(es_list = es, file = graph_file, ylab = ylab, ncol = ncol,
             width = fig_width, height = fig_height)
    write_fig_tex(
      img_file = graph_file,
      caption  = if (!is.null(fig_caption)) fig_caption else graph_title,
      label    = if (!is.null(fig_label)) fig_label
                 else sub("^tab:", "fig:", main_label),
      note     = fig_note
    )
  }

  write_es_longtblr(
    es_list = es, att_list = att, pre_tests = pre,
    n_treated = ntr, n_control = nco, est_methods = ems,
    display_cols = main_display,
    stem = main_stem, caption = main_caption,
    label = main_label, note = main_note
  )

  if (!is.null(dr_key)) {
    write_es_longtblr(
      es_list = es[dr_key], att_list = att[dr_key], pre_tests = pre[dr_key],
      n_treated = ntr[dr_key], n_control = nco[dr_key],
      est_methods = ems[dr_key], display_cols = dr_key,
      stem = dr_stem, caption = dr_caption,
      label = dr_label, note = dr_note
    )
  }
}
