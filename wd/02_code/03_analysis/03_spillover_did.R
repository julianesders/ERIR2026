# Required packages (not in requirements.txt -- install once):
#   install.packages("did")
#   remotes::install_github("kylebutts/didimputation")
#   install.packages("ggplot2")

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
  stop("Cannot determine script path. Run as: Rscript 03_spillover_did.R")
}
root       <- dirname(dirname(dirname(self)))
data_int   <- file.path(root, "01_data", "02_intermediate")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "04_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -- Load panel ----------------------------------------------------------------

panel <- fread(
  file.path(data_final, "emk_inkar_panel_ags8.csv"),
  colClasses = list(character = c("AGS8", "AGS5", "AGS2"))
)
panel <- panel[year >= 2015]

# -- Outcome construction ------------------------------------------------------
# bev_neuzulassungen_p100k, bev_corporate_p100k, bev_private_p100k already
# computed in the Python merge step; no need to reconstruct from raw counts.

panel[, log_bev1 := log1p(bev_neuzulassungen_p100k)]

# ICE placebo: invert EV share applied to existing BEV neuzulassungen p100k.
# Verify scale before trusting these values.
share_med <- median(panel$N_ev_share_overall, na.rm = TRUE)
share_fac <- if (share_med < 0.5) 1 else 100   # 1 if fraction, 100 if pct
cat(sprintf("N_ev_share_overall: median=%.4f -> treating as %s\n",
  share_med, if (share_fac == 1) "fraction [0,1]" else "percentage [0,100]"))
panel[,
  ice_overall_p100k := fifelse(
    N_ev_share_overall / share_fac > 0.001,
    bev_neuzulassungen_p100k * (share_fac / N_ev_share_overall - 1),
    NA_real_
  )
]

# -- Treatment variables -------------------------------------------------------

panel[,
  first_treat_yr := if (any(emk_absorbing == 1L))
    min(year[emk_absorbing == 1L])
  else
    .Machine$integer.max,
  by = AGS8
]
panel[, ever_treated := first_treat_yr != .Machine$integer.max]

# did package convention: gname = 0 for never-treated
panel[, first_treat_cs  := fifelse(ever_treated, as.integer(first_treat_yr), 0L)]
# didimputation convention: gname = Inf for never-treated
panel[, first_treat_bor := fifelse(ever_treated, as.numeric(first_treat_yr), Inf)]

# did package requires an integer unit ID
panel[, ags8_id := as.integer(.GRP), by = AGS8]

cat(sprintf(
  "Panel: %d obs | %d AGS8 | %d ever-treated | %.1f%% onset rate\n",
  nrow(panel), panel[, uniqueN(AGS8)],
  panel[ever_treated == TRUE, uniqueN(AGS8)],
  100 * panel[, mean(emk_absorbing)]
))

# -- Broadcasting sensitivity setup -------------------------------------------
# Broadcast-treated units (AGS5 project broadcast to all Gemeinden in the Kreis)
# are reclassified as never-treated. Direct AGS8-assigned units are unchanged.
# This tests whether any estimated effect survives once the noise from mislabeled
# units is removed. Attenuation toward zero in the full sample is expected if
# broadcast units have no real outcome change.

emk <- fread(
  file.path(data_int, "emk", "emk_ags_matched.csv"),
  colClasses = list(character = c("AGS8", "AGS5")),
  na.strings  = c("", "NA")
)
direct_ags8 <- emk[!is.na(AGS8), unique(AGS8)]
n_direct_treated    <- panel[ever_treated == TRUE & AGS8 %in% direct_ags8,   uniqueN(AGS8)]
n_broadcast_treated <- panel[ever_treated == TRUE & !(AGS8 %in% direct_ags8), uniqueN(AGS8)]
cat(sprintf(
  "Ever-treated AGS8: %d direct, %d broadcast-only\n",
  n_direct_treated, n_broadcast_treated
))

panel_direct <- copy(panel)
panel_direct[
  ever_treated == TRUE & !(AGS8 %in% direct_ags8),
  `:=`(first_treat_cs  = 0L,
       first_treat_bor = Inf,
       ever_treated     = FALSE)
]

# -- Estimation constants and helpers -----------------------------------------

OUTCOMES <- c(
  "bev_neuzulassungen_p100k", "bev_corporate_p100k", "bev_private_p100k",
  "log_bev1", "ice_overall_p100k"
)
ES_MIN <- -4L
ES_MAX <-  4L

# Callaway-Sant'Anna: att_gt + aggte(type="dynamic")
# bstrap=FALSE for speed; switch to TRUE + cband=TRUE for publication.
run_cs <- function(yname, data_dt, control_group) {
  att <- att_gt(
    yname         = yname,
    gname         = "first_treat_cs",
    idname        = "ags8_id",
    tname         = "year",
    data          = as.data.frame(data_dt),
    control_group = control_group,
    anticipation  = 0L,
    clustervars   = "AGS5",
    bstrap        = FALSE,
    print_details = FALSE
  )
  aggte(att, type = "dynamic", min_e = ES_MIN, max_e = ES_MAX,
        balance_e = TRUE, na.rm = TRUE)
}

# Borusyak-Jaravel-Spiess: did_imputation
run_bor <- function(yname, data_dt) {
  did_imputation(
    data        = as.data.frame(data_dt),
    yname       = yname,
    gname       = "first_treat_bor",
    tname       = "year",
    idname      = "ags8_id",
    horizon     = TRUE,
    pretrends   = TRUE,
    cluster_var = "AGS5"
  )
}

# Convert aggte output to plotting data.table
cs_to_dt <- function(es, label) {
  crit <- if (!is.null(es$crit.val.egt)) es$crit.val.egt else 1.96
  data.table(
    e         = es$egt,
    att       = es$att.egt,
    lo_pw     = es$att.egt - 1.96 * es$se.egt,
    hi_pw     = es$att.egt + 1.96 * es$se.egt,
    lo_sim    = es$att.egt - crit * es$se.egt,
    hi_sim    = es$att.egt + crit * es$se.egt,
    estimator = label
  )
}

# Convert didimputation output to plotting data.table.
# Terms: "tau0", "tau1", ... for post; "pre1", "pre2", ... for pre (pre1 = t-1).
bor_to_dt <- function(res, label) {
  dt <- as.data.table(res)
  dt[, e := NA_integer_]
  # "tau0"/"tau1" format (older didimputation)
  dt[grepl("^tau[0-9]+$", term), e := as.integer(sub("tau", "", term))]
  # "pre1"/"pre2" format where pre1 = t-1 (older didimputation)
  dt[grepl("^pre[0-9]+$", term), e := -as.integer(sub("pre", "", term))]
  # plain signed integers "0", "1", "-1", ... (newer didimputation)
  dt[is.na(e) & grepl("^-?[0-9]+$", term), e := as.integer(term)]
  dt <- dt[!is.na(e) & e >= ES_MIN & e <= ES_MAX]
  if (nrow(dt) == 0L) {
    cat(sprintf("    [Borusyak] no event-study terms matched; first terms: %s\n",
                paste(head(res$term, 8), collapse = ", ")))
    return(NULL)
  }
  dt[, .(
    e,
    att    = estimate,
    lo_pw  = estimate - 1.96 * std.error,
    hi_pw  = estimate + 1.96 * std.error,
    lo_sim = estimate - 1.96 * std.error,  # no simultaneous bands in Borusyak
    hi_sim = estimate + 1.96 * std.error,
    estimator = label
  )]
}

es_plot <- function(plot_dt, title, ylab = "ATT (per 100k)") {
  ggplot(plot_dt, aes(x = e, color = estimator, fill = estimator)) +
    geom_hline(yintercept = 0, color = "gray50", linewidth = 0.4) +
    geom_vline(xintercept = -0.5, linetype = "dashed",
               color = "gray50", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lo_pw, ymax = hi_pw), alpha = 0.12, color = NA) +
    geom_line(aes(y = att), linewidth = 0.7) +
    geom_point(aes(y = att), size = 1.8) +
    scale_x_continuous(breaks = ES_MIN:ES_MAX) +
    labs(
      title = title,
      x     = "Years relative to first funding receipt",
      y     = ylab,
      color = NULL, fill = NULL,
      caption = paste0(
        "Shaded band: pointwise 95% CI. Leads left of dashed line adjudicate ",
        "parallel trends.\n2016 cohort has 1 clean pre-period; pre-trend ",
        "evidence is stronger for later cohorts."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.caption = element_text(size = 7))
}

# -- Label maps ----------------------------------------------------------------

OUTCOME_LABELS <- c(
  bev_neuzulassungen_p100k = "BEV new registrations (per 100k population)",
  bev_corporate_p100k      = "BEV new registrations, corporate (per 100k population)",
  bev_private_p100k        = "BEV new registrations, private households (per 100k population)",
  log_bev1                 = "log(1 + BEV new registrations per 100k population)",
  ice_overall_p100k        = "ICE new registrations (per 100k population) — placebo outcome"
)

CTRL_SHORT <- c(
  A_notyet = "Not-yet-treated control\n(sample: ever-treated only)",
  B_never  = "Never-treated control\n(full sample)",
  C_joint  = "Joint: not-yet- and never-treated control\n(full sample)"
)

# -- Control group specifications ---------------------------------------------
# (A) Not-yet-treated only: restrict data to ever-treated, use notyettreated.
#     Local comparison; units are themselves selected (will eventually be treated).
# (B) Never-treated only: nevertreated control group.
#     Larger pool but most negatively selected; pre-trend divergence expected.
# (C) Joint (CS default): notyettreated on full data (both not-yet and never).
#     Valid only if both (A) and (B) satisfy parallel trends.

ctrl_specs <- list(
  A_notyet = list(
    label    = "Control: not-yet-treated (sample restricted to ever-treated Gemeinden)",
    cg       = "notyettreated",
    ever_only = TRUE
  ),
  B_never  = list(
    label    = "Control: never-treated Gemeinden (full sample)",
    cg       = "nevertreated",
    ever_only = FALSE
  ),
  C_joint  = list(
    label    = "Control: not-yet- and never-treated Gemeinden (Callaway-Sant’Anna default)",
    cg       = "notyettreated",
    ever_only = FALSE
  )
)

# -- Estimation loop -----------------------------------------------------------

run_and_plot <- function(panel_use, outcomes, ctrl_specs, file_suffix = "") {
  for (yname in outcomes) {
    ylab <- switch(yname,
      log_bev1                  = "ATT (log BEV neuz. per 100k + 1)",
      ice_overall_p100k         = "ATT (ICE neuz. per 100k) [placebo]",
      bev_neuzulassungen_p100k  = "ATT (BEV neuz. overall per 100k)",
      bev_corporate_p100k       = "ATT (BEV neuz. corporate per 100k)",
      bev_private_p100k         = "ATT (BEV neuz. private per 100k)",
      "ATT (per 100k)"
    )

    meta_parts <- list()

    for (cid in names(ctrl_specs)) {
      spec     <- ctrl_specs[[cid]]
      data_use <- if (spec$ever_only)
        panel_use[ever_treated == TRUE] else panel_use

      cat(sprintf("  %-30s | %-10s | n=%d\n", yname, cid, nrow(data_use)))

      es_cs  <- tryCatch(run_cs(yname, data_use, spec$cg), error = function(e) {
        cat(sprintf("    CS error: %s\n", conditionMessage(e))); NULL
      })
      es_bor <- tryCatch(run_bor(yname, data_use), error = function(e) {
        cat(sprintf("    Borusyak error: %s\n", conditionMessage(e))); NULL
      })

      parts <- Filter(Negate(is.null), list(
        if (!is.null(es_cs))  cs_to_dt(es_cs,  "Callaway-Sant'Anna") else NULL,
        if (!is.null(es_bor)) bor_to_dt(es_bor, "Borusyak et al.")   else NULL
      ))
      if (length(parts) == 0L) next

      dt <- rbindlist(parts)
      dt[, ctrl_id := cid]
      meta_parts[[cid]] <- dt
    }

    if (length(meta_parts) == 0L) next

    meta_dt    <- rbindlist(meta_parts)
    ctrl_order <- intersect(names(ctrl_specs), unique(meta_dt$ctrl_id))
    meta_dt[, facet_lbl := factor(CTRL_SHORT[ctrl_id],
                                   levels = CTRL_SHORT[ctrl_order])]

    out_lbl    <- OUTCOME_LABELS[yname]
    if (is.na(out_lbl)) out_lbl <- yname
    direct_tag <- if (nchar(file_suffix)) "\nDirect project assignment only" else ""

    p <- ggplot(meta_dt, aes(x = e, color = estimator, fill = estimator)) +
      geom_hline(yintercept = 0, color = "gray50", linewidth = 0.4) +
      geom_vline(xintercept = -0.5, linetype = "dashed",
                 color = "gray50", linewidth = 0.4) +
      geom_ribbon(aes(ymin = lo_pw, ymax = hi_pw), alpha = 0.12, color = NA) +
      geom_line(aes(y = att), linewidth = 0.7) +
      geom_point(aes(y = att), size = 1.8) +
      scale_x_continuous(breaks = ES_MIN:ES_MAX) +
      facet_wrap(~facet_lbl, ncol = 3) +
      labs(
        title   = paste0(out_lbl, direct_tag),
        x       = "Years relative to first funding receipt",
        y       = ylab,
        color   = NULL, fill = NULL,
        caption = paste0(
          "Shaded band: pointwise 95% CI. Leads left of dashed line adjudicate ",
          "parallel trends.\n2016 cohort has 1 clean pre-period; pre-trend ",
          "evidence is stronger for later cohorts."
        )
      ) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom",
            plot.caption    = element_text(size = 7),
            strip.text      = element_text(size = 8, face = "bold"),
            panel.spacing   = unit(1, "lines"))

    fname <- sprintf("es_%s%s.pdf", yname, file_suffix)
    ggsave(file.path(out_dir, fname), p, width = 15, height = 5)
  }
}

cat("\n=== Primary estimation (all treated, broadcast included) ===\n")
run_and_plot(panel, OUTCOMES, ctrl_specs)

cat("\n=== Broadcasting sensitivity (direct AGS8-assigned only) ===\n")
run_and_plot(panel_direct, c("bev_neuzulassungen_p100k"), ctrl_specs, "_direct")

cat(sprintf("\nAll plots written to: %s\n", out_dir))
