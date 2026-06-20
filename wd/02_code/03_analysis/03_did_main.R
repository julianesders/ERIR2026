# ───────────────────────────────────────────────────────────────────────────────
# 03_did_main.R   — Headline staggered-DiD estimates
#
# Two estimators (co-primary):
#   - BJS (didimputation): conditional adjustment via unit + year FE
#   - CS  (did, est_method="dr" without xformla): unconditional parallel
#     trends. The covariate-adjusted CS variant lives in 04_did_robustness.R.
#
# Two outcome scales for every estimate:
#   - level rate  : bev_*_p100k         (additional registrations / 100k)
#   - log1p rate  : log1p_bev_*_p100k   (≈ percent change in the rate)
#
# Two frames: broad (main), direct (sharp).
# Three outcome tiers: headline (overall flow), secondaries (corp/priv/stock),
# placebo (ICE flow). Selection-control / anticipation variants live in 04.
#
# Outputs (04_results/03_did_main/):
#   est_att_main.csv                 long: outcome × frame × estimator × scale
#   tab_att_main.{tex,csv}           headline outcome rows
#   es_headline_<frame>.pdf          2-panel: raw & log1p (BJS+CS overlay)
#   es_secondary_<frame>_<scale>.pdf 3-panel: corp / priv / stock
#   es_placebo_<frame>.pdf           2-panel: raw & log1p ICE
#   composition_<frame>.csv          per-horizon cohort/treated counts
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
  stop("Cannot determine script path. Run as: Rscript 03_did_main.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "04_results", "03_did_main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))

# -- Constants -----------------------------------------------------------------

CS_BITERS <- 2000L

HEADLINE   <- "bev_neuzulassungen_p100k"
SECONDARY  <- c("bev_corporate_p100k", "bev_private_p100k", "bev_stock_p100k")
PLACEBO    <- "ice_neuzulassungen_p100k"
ALL_LEVEL  <- c(HEADLINE, SECONDARY, PLACEBO)

# -- Load frames ---------------------------------------------------------------

.read_frame <- function(p) fread(p,
  colClasses = list(character = c("AGS8", "AGS5", "AGS2", "treat_type")))
frames <- list(
  broad  = .read_frame(file.path(data_final, "frame_did_broad.csv")),
  direct = .read_frame(file.path(data_final, "frame_did_direct.csv"))
)
for (nm in names(frames)) {
  d <- frames[[nm]]
  stopifnot(is.integer(d$ags8_id),
            is.integer(d$gname_cs),
            is.numeric(d$gname_bjs))
}

# -- A1 guard ------------------------------------------------------------------

for (nm in names(frames)) {
  d  <- frames[[nm]]
  g  <- d$gname_bjs; y <- d$year
  untreated_obs <- sum(is.na(g) | g == BJS_NEVER | y < g)
  never_units   <- length(unique(d$ags8_id[is.na(g) | g == BJS_NEVER]))
  cat(sprintf("BJS [%s]: untreated_obs=%d (%.0f%%); never_units=%d\n",
              nm, untreated_obs, 100 * untreated_obs / nrow(d), never_units))
  stopifnot(untreated_obs > 0.5 * nrow(d), never_units > 1000L)
}

# -- Estimators (single cell wrappers) ----------------------------------------

run_cs <- function(yname, dat) {
  att_gt(
    yname                  = yname,
    gname                  = "gname_cs",
    idname                 = "ags8_id",
    tname                  = "year",
    data                   = as.data.frame(dat),
    control_group          = "nevertreated",
    anticipation           = 0L,
    xformla                = NULL,
    est_method             = "dr",
    clustervars            = "AGS5",
    bstrap                 = TRUE,
    biters                 = CS_BITERS,
    allow_unbalanced_panel = TRUE,
    print_details          = FALSE
  )
}


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

# -- ES extractors -------------------------------------------------------------

cs_es_dt <- function(cs_obj, es_max_eff) {
  es <- tryCatch(aggte(cs_obj, type = "dynamic",
                       min_e = ES_MIN, max_e = es_max_eff,
                       balance_e = NULL, na.rm = TRUE, cband = TRUE),
                 error = function(e) NULL)
  if (is.null(es)) return(NULL)
  data.table(
    e         = es$egt,
    estimate  = es$att.egt,
    se        = es$se.egt,
    ci_lo     = es$att.egt - 1.96 * es$se.egt,
    ci_hi     = es$att.egt + 1.96 * es$se.egt,
    estimator = "CS"
  )[e >= ES_MIN & e <= es_max_eff]
}

bjs_es_dt <- function(bjs_res, es_max_eff) {
  if (is.null(bjs_res)) return(NULL)
  dt <- as.data.table(bjs_res)
  dt[, e := NA_integer_]
  dt[grepl("^tau[0-9]+$", term), e := as.integer(sub("tau", "", term))]
  dt[grepl("^pre[0-9]+$", term), e := -as.integer(sub("pre", "", term))]
  dt[is.na(e) & grepl("^-?[0-9]+$", term), e := as.integer(term)]
  dt <- dt[!is.na(e) & e >= ES_MIN & e <= es_max_eff]
  if (nrow(dt) == 0L) return(NULL)
  dt[, .(e, estimate, se = std.error,
         ci_lo = estimate - 1.96 * std.error,
         ci_hi = estimate + 1.96 * std.error,
         estimator = "BJS")]
}

# Per-horizon composition from att_gt object (B7)
horizon_composition <- function(cs_obj, dat) {
  if (is.null(cs_obj)) return(NULL)
  gt <- data.table(g = cs_obj$group, t = cs_obj$t)
  gt[, e := t - g]
  cohort_sizes <- dat[gname_cs > 0L, .(n_units = uniqueN(AGS8)), by = gname_cs]
  setnames(cohort_sizes, "gname_cs", "g")
  gt <- merge(gt, cohort_sizes, by = "g", all.x = TRUE)
  gt[, .(n_cohorts = uniqueN(g), n_treated_units = sum(n_units)), by = e]
}

# -- Single-cell estimate (returns ATT row + ES dt) ---------------------------

estimate_cell <- function(dat, yname, frame_nm, scale) {
  dat_y <- dat[!is.na(get(yname))]
  if (nrow(dat_y) == 0L || uniqueN(dat_y[gname_cs > 0L, AGS8]) < 5L)
    return(list(att = NULL, es = NULL, comp = NULL))
  es_max_eff <- es_max_data_driven(dat_y, yname, "gname_cs")
  n_treated  <- uniqueN(dat_y[gname_cs > 0L, AGS8])
  n_control  <- uniqueN(dat_y[gname_cs == 0L, AGS8])
  cat(sprintf("  [%s/%s/%s] n=%d, es_max=%d\n",
              frame_nm, scale, yname, nrow(dat_y), es_max_eff))

  # CS
  cs_obj <- tryCatch(run_cs(yname, dat_y), error = function(e) {
    cat("    CS error:", conditionMessage(e), "\n"); NULL })
  cs_simple <- if (!is.null(cs_obj))
    tryCatch(aggte(cs_obj, type = "simple", na.rm = TRUE),
             error = function(e) NULL) else NULL
  cs_es <- if (!is.null(cs_obj)) cs_es_dt(cs_obj, es_max_eff) else NULL

  # BJS
  bjs_es_raw  <- tryCatch(run_bjs(yname, dat_y, horizon = TRUE),
                          error = function(e) {
                            cat("    BJS ES err:", conditionMessage(e), "\n")
                            NULL })
  bjs_att_raw <- tryCatch(run_bjs(yname, dat_y, horizon = NULL),
                          error = function(e) NULL)
  bjs_es <- bjs_es_dt(bjs_es_raw, max(es_max_eff, ES_MAX))

  rows <- list()
  if (!is.null(cs_simple))
    rows[[length(rows) + 1L]] <- data.table(
      frame = frame_nm, outcome = yname, scale = scale, estimator = "CS",
      estimate = cs_simple$overall.att, se = cs_simple$overall.se,
      ci_lo = cs_simple$overall.att - 1.96 * cs_simple$overall.se,
      ci_hi = cs_simple$overall.att + 1.96 * cs_simple$overall.se,
      n_treated = n_treated, n_control = n_control
    )
  if (!is.null(bjs_att_raw)) {
    r <- as.data.table(bjs_att_raw)[term == "treat"][1L]
    if (!is.null(r) && nrow(r) > 0L)
      rows[[length(rows) + 1L]] <- data.table(
        frame = frame_nm, outcome = yname, scale = scale, estimator = "BJS",
        estimate = r$estimate, se = r$std.error,
        ci_lo = r$estimate - 1.96 * r$std.error,
        ci_hi = r$estimate + 1.96 * r$std.error,
        n_treated = n_treated, n_control = n_control
      )
  }
  att_dt <- if (length(rows)) rbindlist(rows) else NULL

  es_parts <- Filter(Negate(is.null), list(cs_es, bjs_es))
  es_dt <- if (length(es_parts))
    rbindlist(es_parts, use.names = TRUE, fill = TRUE)[
      , `:=`(frame = frame_nm, outcome = yname, scale = scale)] else NULL

  list(att = att_dt, es = es_dt,
       comp = horizon_composition(cs_obj, dat_y))
}

# -- Outcome × scale grid ------------------------------------------------------

outcome_grid <- rbindlist(list(
  data.table(outcome = ALL_LEVEL,                      scale = "level"),
  data.table(outcome = paste0("log1p_", ALL_LEVEL),    scale = "log1p")
))

# -- Run loop -----------------------------------------------------------------

att_all <- list()
es_all  <- list()
comp_rows <- list()

for (fr_nm in names(frames)) {
  dat <- frames[[fr_nm]]
  cat(sprintf("\n=== %s frame ===\n", fr_nm))
  for (i in seq_len(nrow(outcome_grid))) {
    yname <- outcome_grid$outcome[i]
    scale <- outcome_grid$scale[i]
    if (!(yname %in% names(dat))) {
      cat(sprintf("  skip: %s not in %s frame\n", yname, fr_nm)); next
    }
    res <- estimate_cell(dat, yname, fr_nm, scale)
    if (!is.null(res$att)) att_all[[length(att_all) + 1L]] <- res$att
    if (!is.null(res$es))  es_all[[length(es_all) + 1L]]   <- res$es
    if (!is.null(res$comp))
      comp_rows[[length(comp_rows) + 1L]] <-
        res$comp[, `:=`(frame = fr_nm, outcome = yname, scale = scale)]
  }
}

att_dt <- rbindlist(att_all, use.names = TRUE, fill = TRUE)
es_dt  <- rbindlist(es_all,  use.names = TRUE, fill = TRUE)
comp_dt <- if (length(comp_rows))
  rbindlist(comp_rows, use.names = TRUE, fill = TRUE) else data.table()

fwrite(att_dt, file.path(out_dir, "est_att_main.csv"))
if (nrow(comp_dt))
  fwrite(comp_dt, file.path(out_dir, "composition.csv"))

# -- Headline table (level + log1p, both frames, both estimators) ------------

tab_long <- att_dt[outcome %in% c(HEADLINE, paste0("log1p_", HEADLINE))]
tab <- dcast(tab_long, frame + scale ~ estimator,
             value.var = c("estimate", "se"))
fwrite(tab, file.path(out_dir, "tab_att_main.csv"))

.fmt <- function(x, d = 3L) ifelse(is.na(x), "--",
                                   sprintf(paste0("%.", d, "f"), x))
tex_lines <- c(
  "\\begin{tabular}{llrrrr}",
  "\\hline",
  "Frame & Scale & BJS ATT & BJS SE & CS ATT & CS SE \\\\",
  "\\hline",
  tab[, sprintf("%s & %s & %s & %s & %s & %s \\\\",
                frame, scale,
                .fmt(estimate_BJS), .fmt(se_BJS),
                .fmt(estimate_CS),  .fmt(se_CS))],
  "\\hline",
  "\\end{tabular}"
)
writeLines(tex_lines, file.path(out_dir, "tab_att_main.tex"))

# -- Plots --------------------------------------------------------------------
# Layout rules: ≤3 panels per image.
#   Headline (1 outcome × 2 scales) → 1 figure per frame, 2 panels (raw, log1p)
#   Secondary (3 outcomes × 2 scales) → 1 figure per frame per scale, 3 panels
#   Placebo (1 outcome × 2 scales) → 1 figure per frame, 2 panels

es_plot <- function(dt, file, ylab, ttl, facet_var = "outcome") {
  if (is.null(dt) || nrow(dt) == 0L) return(invisible(NULL))
  dt <- copy(dt)
  if (facet_var == "scale_outcome") {
    dt[, panel := sprintf("%s\n(%s)",
                          ifelse(outcome %in% names(OUTCOME_LABELS),
                                 OUTCOME_LABELS[outcome], outcome),
                          scale)]
  } else if (facet_var == "outcome") {
    dt[, panel := ifelse(outcome %in% names(OUTCOME_LABELS),
                         OUTCOME_LABELS[outcome], outcome)]
  } else {
    dt[, panel := get(facet_var)]
  }
  n_facets <- uniqueN(dt$panel)
  p <- ggplot(dt, aes(x = e, color = estimator, fill = estimator)) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = -0.5, linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15,
                color = NA) +
    geom_line(aes(y = estimate), linewidth = 0.7) +
    geom_point(aes(y = estimate), size = 1.6) +
    facet_wrap(~ panel, ncol = n_facets, scales = "free_y") +
    labs(title = ttl, x = "Years relative to first funding receipt",
         y = ylab, color = NULL, fill = NULL,
         caption = "Shaded band: pointwise 95% CI. CS unconditional; BJS imputation.") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))
  ggsave(file, p, width = 3.5 * n_facets + 1, height = 4.2)
}

for (fr_nm in names(frames)) {
  # Headline: 2 panels (level + log1p) for HEADLINE outcome
  head_dt <- es_dt[frame == fr_nm &
                   outcome %in% c(HEADLINE, paste0("log1p_", HEADLINE))]
  es_plot(head_dt,
          file.path(out_dir, sprintf("es_headline_%s.pdf", fr_nm)),
          ylab = "ATT",
          ttl  = sprintf("Headline — %s frame", fr_nm),
          facet_var = "scale_outcome")

  # Placebo: 2 panels (level + log1p) for PLACEBO outcome
  plac_dt <- es_dt[frame == fr_nm &
                   outcome %in% c(PLACEBO, paste0("log1p_", PLACEBO))]
  es_plot(plac_dt,
          file.path(out_dir, sprintf("es_placebo_%s.pdf", fr_nm)),
          ylab = "ATT",
          ttl  = sprintf("Placebo (ICE) — %s frame", fr_nm),
          facet_var = "scale_outcome")

  # Secondary: per scale, 3 panels (corp / priv / stock)
  for (sc in c("level", "log1p")) {
    sec_cols <- if (sc == "level") SECONDARY else paste0("log1p_", SECONDARY)
    sec_dt <- es_dt[frame == fr_nm & scale == sc & outcome %in% sec_cols]
    es_plot(sec_dt,
            file.path(out_dir,
                      sprintf("es_secondary_%s_%s.pdf", fr_nm, sc)),
            ylab = "ATT",
            ttl  = sprintf("Secondary outcomes — %s frame (%s)", fr_nm, sc),
            facet_var = "outcome")
  }
}

cat(sprintf("\nMain DiD outputs -> %s\n", out_dir))
