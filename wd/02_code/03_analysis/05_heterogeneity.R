# ───────────────────────────────────────────────────────────────────────────────
# 05_heterogeneity.R   — Inequality payoff
#
# CS (dr, never-ctrl, broad frame) estimated within Steuerkraft (tax capacity)
# terciles, plus treat_type heterogeneity (direct vs broadcast-only treated
# units, broad frame). Steuerkraft is the only wealth ranking used; Kaufkraft
# is deliberately excluded from the DiD heterogeneity.
#
# Outputs (04_results/05_heterogeneity/):
#   did_heterogeneity_terciles.{tex,csv}    ATT × Steuerkraft tercile
#   fig_att_by_tercile.pdf                  event-studies (facet) + dot-whisker
#   did_heterogeneity_treat_type.{tex,csv}  direct vs broadcast-only ATT
# ───────────────────────────────────────────────────────────────────────────────

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
  stop("Cannot determine script path. Run as: Rscript 05_heterogeneity.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "04_results", "05_heterogeneity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))

CS_BITERS <- 1999L
OUTCOME   <- "bev_neuzulassungen_p100k"

.read_frame <- function(p) fread(p,
  colClasses = list(character = c("AGS8", "AGS5", "AGS2", "treat_type")))
frame_broad  <- .read_frame(file.path(data_final, "frame_did_broad.csv"))
frame_direct <- .read_frame(file.path(data_final, "frame_did_direct.csv"))

# Headline 03 spec is unconditional CS (xformla = NULL); the heterogeneity
# stratification keeps the same identifying assumption per tercile.
xformla_cov <- NULL

# -- Tercile-stratified CS ----------------------------------------------------
# Treated AND control units split by the same baseline-tercile.

run_cs_terc <- function(dat, terc_col) {
  out_rows <- list()
  es_rows  <- list()
  for (t in sort(unique(na.omit(dat[[terc_col]])))) {
    sub <- dat[get(terc_col) == t]
    n_tr <- uniqueN(sub[gname_cs > 0L, AGS8])
    n_co <- uniqueN(sub[gname_cs == 0L, AGS8])
    if (n_tr < 5L || n_co < 5L) {
      cat(sprintf("  tercile %d: skipped (n_tr=%d, n_co=%d)\n", t, n_tr, n_co))
      next
    }
    cs <- tryCatch(att_gt(
      yname = OUTCOME, gname = "gname_cs", idname = "ags8_id", tname = "year",
      data = as.data.frame(sub), control_group = "nevertreated",
      anticipation = 0L, xformla = xformla_cov, est_method = "dr",
      clustervars = "AGS5", bstrap = TRUE, biters = CS_BITERS,
      allow_unbalanced_panel = TRUE,
      print_details = FALSE
    ), error = function(e) {
      cat("  CS err tercile", t, ":", conditionMessage(e), "\n"); NULL })
    if (is.null(cs)) next
    s  <- tryCatch(aggte(cs, type = "simple",  na.rm = TRUE),
                   error = function(e) NULL)
    dy <- tryCatch(aggte(cs, type = "dynamic",
                         min_e = ES_MIN, max_e = ES_MAX,
                         balance_e = NULL, na.rm = TRUE, cband = TRUE),
                   error = function(e) NULL)
    if (!is.null(s))
      out_rows[[length(out_rows) + 1L]] <- data.table(
        tercile = t, n_treated = n_tr, n_control = n_co,
        att = s$overall.att, se = s$overall.se
      )
    if (!is.null(dy))
      es_rows[[length(es_rows) + 1L]] <- data.table(
        tercile = t,
        e         = dy$egt,
        att       = dy$att.egt,
        se        = dy$se.egt
      )
  }
  list(att = rbindlist(out_rows, fill = TRUE),
       es  = rbindlist(es_rows, fill = TRUE))
}

res_sk <- run_cs_terc(frame_broad, "sk_base_terc")

terc_dt <- copy(res_sk$att)[, rank := "sk_base"]
terc_dt[, ci_lo := att - 1.96 * se]
terc_dt[, ci_hi := att + 1.96 * se]
fwrite(terc_dt, file.path(out_dir, "did_heterogeneity_terciles.csv"))

write_longtblr(
  stem        = file.path(out_dir, "did_heterogeneity_terciles"),
  caption     = "Heterogeneity by Baseline Tax-Capacity Tercile: CS Overall ATT",
  label       = "tab:did_heterogeneity_terciles",
  note        = paste0(
    "Callaway--Sant'Anna (2021) overall ATT (\\texttt{bev\\_neuzulassungen\\_p100k}) ",
    "estimated separately within each tercile of baseline tax capacity ",
    "(\\texttt{sk\\_base}, Steuerkraft). ",
    "Terciles defined on treated + never-treated units; never-treated control group. ",
    "Broad frame, broad treatment definition. ",
    "95\\% CI: pointwise 1.96-SE (cluster bootstrap on AGS5, $B=",
    CS_BITERS, "$)."
  ),
  colspec     = "l r r r r r r",
  header_rows = "Rank & Tercile & $N_{\\text{treated}}$ & $N_{\\text{control}}$ & ATT & SE & 95\\% CI \\\\",
  body_rows   = terc_dt[, sprintf(
    "%s & %d & %d & %d & %.3f & %.3f & [%.3f,\\,%.3f] \\\\",
    rank, tercile, n_treated, n_control, att, se, ci_lo, ci_hi)],
  footer_rows = character(0)
)

# Difference top vs bottom tercile (sk): SE-based normal approx (independence
# across split-sample bootstraps); full split-sample bootstrap optional.
diff_sk <- if (nrow(res_sk$att) >= 2L) {
  hi <- res_sk$att[tercile == max(tercile)]
  lo <- res_sk$att[tercile == min(tercile)]
  data.table(
    rank = "sk_base",
    diff = hi$att - lo$att,
    se   = sqrt(hi$se^2 + lo$se^2)
  )[, ci_lo := diff - 1.96 * se][, ci_hi := diff + 1.96 * se][]
} else NULL
if (!is.null(diff_sk))
  fwrite(diff_sk, file.path(out_dir, "tab_diff_top_bottom_sk.csv"))

# -- Figure: event-studies (facet by Steuerkraft tercile) ---------------------

if (nrow(res_sk$es) > 0L) {
  es_dt <- copy(res_sk$es)
  es_dt[, ci_lo := att - 1.96 * se]
  es_dt[, ci_hi := att + 1.96 * se]
  es_dt[, facet := sprintf("sk tercile %d", tercile)]
  p_es <- ggplot(es_dt, aes(x = e, y = att)) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = -0.5, linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.18,
                fill = "#1b6ca8") +
    geom_line(linewidth = 0.7, color = "#1b6ca8") +
    geom_point(size = 1.8, color = "#1b6ca8") +
    facet_wrap(~ facet, ncol = 3) +
    scale_x_continuous(breaks = ES_MIN:ES_MAX) +
    labs(x = "Years relative to first funding receipt",
         y = "ATT (BEV new registrations p100k)",
         caption = "CS-dr, never-treated control, broad frame.") +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(out_dir, "fig_att_by_tercile.pdf"), p_es,
         width = 15, height = 5)
}

# -- Treat-type heterogeneity: broad-frame ATT for direct vs broadcast_only ---
# Run CS on the broad frame twice: once dropping broadcast_only from the
# treated set, once dropping direct. Compare the resulting ATTs.

dat_dir_only  <- frame_broad[treat_type %in% c("direct", "never")]
dat_brod_only <- frame_broad[treat_type %in% c("broadcast_only", "never")]

run_simple_att <- function(d, lbl) {
  if (uniqueN(d[gname_cs > 0L, AGS8]) < 5L) return(NULL)
  cs <- tryCatch(att_gt(
    yname = OUTCOME, gname = "gname_cs", idname = "ags8_id", tname = "year",
    data = as.data.frame(d), control_group = "nevertreated",
    anticipation = 0L, xformla = xformla_cov, est_method = "dr",
    clustervars = "AGS5", bstrap = TRUE, biters = CS_BITERS,
    allow_unbalanced_panel = TRUE,
    print_details = FALSE
  ), error = function(e) NULL)
  if (is.null(cs)) return(NULL)
  s <- tryCatch(aggte(cs, type = "simple", na.rm = TRUE),
                error = function(e) NULL)
  if (is.null(s)) return(NULL)
  data.table(
    treated_subset = lbl,
    att = s$overall.att, se = s$overall.se,
    n_treated = uniqueN(d[gname_cs > 0L, AGS8]),
    n_control = uniqueN(d[gname_cs == 0L, AGS8])
  )
}

tt_dt <- rbindlist(list(
  run_simple_att(dat_dir_only,  "direct only"),
  run_simple_att(dat_brod_only, "broadcast only")
), fill = TRUE)
if (nrow(tt_dt)) {
  tt_dt[, ci_lo := att - 1.96 * se][, ci_hi := att + 1.96 * se]
  fwrite(tt_dt, file.path(out_dir, "did_heterogeneity_treat_type.csv"))
  write_longtblr(
    stem        = file.path(out_dir, "did_heterogeneity_treat_type"),
    caption     = "Heterogeneity by Treatment Type: Direct vs Broadcast-only",
    label       = "tab:did_heterogeneity_treat_type",
    note        = paste0(
      "Callaway--Sant'Anna (2021) overall ATT (\\texttt{bev\\_neuzulassungen\\_p100k}) ",
      "estimated on the broad frame restricted by treated subset: ",
      "``direct only'' keeps \\texttt{treat\\_type\\,=\\,direct} as treated ",
      "(municipalities with their own EMK project); ",
      "``broadcast only'' keeps \\texttt{treat\\_type\\,=\\,broadcast\\_only} ",
      "(covered by a county-level project but no own project). ",
      "In both cases the never-treated pool is unchanged. ",
      "Unconditional parallel trends, never-treated control group. ",
      "95\\% CI: pointwise 1.96-SE (cluster bootstrap on AGS5, $B=",
      CS_BITERS, "$)."
    ),
    colspec     = "l r r r r r",
    header_rows = "Treated subset & $N_{\\text{treated}}$ & $N_{\\text{control}}$ & ATT & SE & 95\\% CI \\\\",
    body_rows   = tt_dt[, sprintf(
      "%s & %d & %d & %.3f & %.3f & [%.3f,\\,%.3f] \\\\",
      treated_subset, n_treated, n_control, att, se, ci_lo, ci_hi)],
    footer_rows = character(0)
  )
}

cat(sprintf("\nHeterogeneity outputs -> %s\n", out_dir))
