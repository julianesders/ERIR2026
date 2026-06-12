# ───────────────────────────────────────────────────────────────────────────────
# 05_heterogeneity.R   — Inequality payoff
#
# CS (dr, never-ctrl, broad frame) estimated within Kaufkraft terciles and
# Steuerkraft terciles, plus treat_type heterogeneity (direct vs broadcast-
# only treated units, broad frame).
#
# Outputs (04_results/05_heterogeneity/):
#   tab_att_terciles.{tex,csv}     ATT × tercile × rank
#   fig_att_by_tercile.pdf         event-studies (facet) + dot-whisker
#   tab_att_treat_type.{tex,csv}   direct vs broadcast-only ATT
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

frame_broad  <- readRDS(file.path(data_final, "frame_did_broad.rds"))
frame_direct <- readRDS(file.path(data_final, "frame_did_direct.rds"))

xformla_cov <- ~ kk_base_z + sk_base_z + dens_base_z + green_base_z +
                  bev_base_z + chg_base_z

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

res_kk <- run_cs_terc(frame_broad, "kk_base_terc")
res_sk <- run_cs_terc(frame_broad, "sk_base_terc")

terc_dt <- rbindlist(list(
  res_kk$att[, rank := "kk_base"],
  res_sk$att[, rank := "sk_base"]
), fill = TRUE)
terc_dt[, ci_lo := att - 1.96 * se]
terc_dt[, ci_hi := att + 1.96 * se]
fwrite(terc_dt, file.path(out_dir, "tab_att_terciles.csv"))

tex_lines <- c(
  "\\begin{tabular}{llrrrrr}",
  "\\hline",
  "Rank & Tercile & N treated & N control & ATT & SE & 95\\% CI \\\\",
  "\\hline",
  terc_dt[, sprintf("%s & %d & %d & %d & %.3f & %.3f & [%.3f, %.3f] \\\\",
                    rank, tercile, n_treated, n_control,
                    att, se, ci_lo, ci_hi)],
  "\\hline",
  "\\end{tabular}"
)
writeLines(tex_lines, file.path(out_dir, "tab_att_terciles.tex"))

# Bootstrap difference top vs bottom tercile (kk): rerun with cluster bootstrap
# Heuristic — use the SE-based normal approx for now; full bootstrap optional.
diff_kk <- if (nrow(res_kk$att) >= 2L) {
  hi <- res_kk$att[tercile == max(tercile)]
  lo <- res_kk$att[tercile == min(tercile)]
  data.table(
    rank = "kk_base",
    diff = hi$att - lo$att,
    se   = sqrt(hi$se^2 + lo$se^2)
  )[, ci_lo := diff - 1.96 * se][, ci_hi := diff + 1.96 * se][]
} else NULL
if (!is.null(diff_kk))
  fwrite(diff_kk, file.path(out_dir, "tab_diff_top_bottom_kk.csv"))

# -- Figure: event-studies (facet by tercile) + dot-whisker (kk) --------------

if (nrow(res_kk$es) > 0L) {
  es_dt <- copy(res_kk$es)
  es_dt[, ci_lo := att - 1.96 * se]
  es_dt[, ci_hi := att + 1.96 * se]
  es_dt[, facet := sprintf("kk tercile %d", tercile)]
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
  fwrite(tt_dt, file.path(out_dir, "tab_att_treat_type.csv"))
  tex_lines <- c(
    "\\begin{tabular}{lrrrrrr}",
    "\\hline",
    "Treated subset & N treated & N control & ATT & SE & 95\\% CI \\\\",
    "\\hline",
    tt_dt[, sprintf("%s & %d & %d & %.3f & %.3f & [%.3f, %.3f] \\\\",
                    treated_subset, n_treated, n_control,
                    att, se, ci_lo, ci_hi)],
    "\\hline",
    "\\end{tabular}"
  )
  writeLines(tex_lines, file.path(out_dir, "tab_att_treat_type.tex"))
}

cat(sprintf("\nHeterogeneity outputs -> %s\n", out_dir))
