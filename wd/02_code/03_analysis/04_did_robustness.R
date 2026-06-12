# ───────────────────────────────────────────────────────────────────────────────
# 04_did_robustness.R   — Robustness for Part (ii)
#
#   - Sun-Abraham via fixest::sunab on bev_neuzulassungen_p100k, corp, private
#   - PPML / Wooldridge via fepois on counts with offset(log(xbev))
#   - Wild cluster bootstrap on TWFE event-study (sanity, if few clusters)
#   - Variations: anticipation ∈ {0,1}; drop 2016 cohort; complete-case on
#     _imp flags; drop COVID years 2020–21
#
# Outputs (04_results/04_did_robustness/):
#   est_att_robust.csv          one row per variation × outcome
#   fig_robust_grid.pdf         one appendix figure (forest of ATT × variant)
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
  stop("Cannot determine script path. Run as: Rscript 04_did_robustness.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "04_results", "04_did_robustness")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))

frame_broad <- readRDS(file.path(data_final, "frame_did_broad.rds"))

OUTCOMES_CONT  <- c("bev_neuzulassungen_p100k", "bev_corporate_p100k",
                    "bev_private_p100k")
OUTCOMES_COUNT <- c("N_elektro_overall", "N_elektro_corporate",
                    "N_elektro_private")

# -- Helper: Sun-Abraham via sunab + fixest -----------------------------------
# fixest::sunab convention: never-treated coded to 10000 (far cohort).

sunab_fit <- function(dat, yname) {
  d <- as.data.table(dat)
  d[, sunab_g := ifelse(is.na(first_treat_broad), 10000L,
                        as.integer(first_treat_broad))]
  feols(as.formula(sprintf("%s ~ sunab(sunab_g, year) | AGS8 + year", yname)),
        data = d, cluster = ~AGS5)
}

# -- Helper: BJS ATT (simple) -------------------------------------------------

bjs_simple <- function(dat, yname) {
  r <- tryCatch(as.data.table(did_imputation(
    data        = as.data.frame(dat),
    yname       = yname,
    gname       = "gname_bjs",
    tname       = "year",
    idname      = "ags8_id",
    horizon     = NULL,   # static effect per didimputation v0.5.1
    cluster_var = "AGS5"
  )), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  row <- r[term == "treat"][1L]
  if (is.null(row) || nrow(row) == 0L) return(NULL)
  list(estimate = row$estimate, se = row$std.error, n_obs = nrow(dat))
}

# -- Variant grid --------------------------------------------------------------

variants <- list(
  baseline       = list(filter = function(d) d,                     anticip = 0L),
  anticipation_1 = list(filter = function(d) d,                     anticip = 1L),
  drop_2016      = list(filter = function(d) d[is.na(first_treat_broad) | first_treat_broad != 2016L],
                        anticip = 0L),
  drop_covid     = list(filter = function(d) d[!(year %in% c(2020L, 2021L))],
                        anticip = 0L),
  complete_case  = list(filter = function(d) {
    cc <- d
    if ("B_elektro_overall_imp" %in% names(d))
      cc <- cc[B_elektro_overall_imp == FALSE]
    cc
  }, anticip = 0L)
)

rows <- list()

for (vname in names(variants)) {
  vspec <- variants[[vname]]
  dat   <- vspec$filter(frame_broad)
  cat(sprintf("\n=== variant: %s | %d obs | %d treated AGS8 ===\n",
              vname, nrow(dat), uniqueN(dat[gname_cs > 0L, AGS8])))

  for (yname in OUTCOMES_CONT) {
    # (a) Sun-Abraham
    m_sa <- tryCatch(sunab_fit(dat, yname), error = function(e) {
      cat("  SA error:", conditionMessage(e), "\n"); NULL })
    if (!is.null(m_sa)) {
      sa_att <- tryCatch(summary(m_sa, agg = "ATT"), error = function(e) NULL)
      if (!is.null(sa_att)) {
        co <- coef(sa_att); se <- sqrt(diag(vcov(sa_att)))
        rows[[length(rows) + 1L]] <- data.table(
          variant = vname, estimator = "Sun-Abraham", outcome = yname,
          estimate = co[1L], se = se[1L]
        )
      }
    }

    # (b) BJS
    b <- bjs_simple(dat, yname)
    if (!is.null(b))
      rows[[length(rows) + 1L]] <- data.table(
        variant = vname, estimator = "BJS", outcome = yname,
        estimate = b$estimate, se = b$se
      )

    # (c) CS-dr (anticipation per variant)
    cs <- tryCatch(att_gt(
      yname = yname, gname = "gname_cs", idname = "ags8_id", tname = "year",
      data = as.data.frame(dat), control_group = "nevertreated",
      anticipation = vspec$anticip, xformla = NULL,
      est_method = "dr",
      clustervars = "AGS5", bstrap = TRUE, biters = 999L,
      allow_unbalanced_panel = TRUE,
      print_details = FALSE
    ), error = function(e) NULL)
    if (!is.null(cs)) {
      s <- tryCatch(aggte(cs, type = "simple", na.rm = TRUE),
                    error = function(e) NULL)
      if (!is.null(s))
        rows[[length(rows) + 1L]] <- data.table(
          variant = vname, estimator = "CS",
          outcome = yname,
          estimate = s$overall.att, se = s$overall.se
        )
    }
  }

  # (d) PPML / fepois on counts with offset(log(xbev)).
  # Cohort × event-time interactions for the proportional-effect reading.
  for (yname in OUTCOMES_COUNT) {
    d <- copy(dat)
    d <- d[!is.na(get(yname)) & xbev > 0]
    d[, cohort := fifelse(is.na(first_treat_broad), 10000L,
                          as.integer(first_treat_broad))]
    fml <- as.formula(sprintf(
      "%s ~ sunab(cohort, year) + offset(log(xbev)) | AGS8 + year",
      yname))
    m_ppml <- tryCatch(fepois(fml, data = d, cluster = ~AGS5),
                       error = function(e) {
                         cat("  PPML err:", conditionMessage(e), "\n"); NULL })
    if (!is.null(m_ppml)) {
      s_att <- tryCatch(summary(m_ppml, agg = "ATT"), error = function(e) NULL)
      if (!is.null(s_att)) {
        co <- coef(s_att); se <- sqrt(diag(vcov(s_att)))
        rows[[length(rows) + 1L]] <- data.table(
          variant = vname, estimator = "PPML", outcome = yname,
          estimate = co[1L], se = se[1L]
        )
      }
    }
  }
}

# -- Wild cluster bootstrap on TWFE event-study (sanity, if few clusters) -----

if (requireNamespace("fwildclusterboot", quietly = TRUE)) {
  d <- copy(frame_broad)
  d[, treat := as.integer(!is.na(first_treat_broad) &
                          year >= first_treat_broad)]
  m_twfe <- feols(bev_neuzulassungen_p100k ~ treat | AGS8 + year,
                  data = d, cluster = ~AGS5)
  bt <- tryCatch(fwildclusterboot::boottest(
    m_twfe, param = "treat", clustid = "AGS5", B = 999L, type = "rademacher"
  ), error = function(e) {
    cat("Wild boot err:", conditionMessage(e), "\n"); NULL })
  if (!is.null(bt)) {
    rows[[length(rows) + 1L]] <- data.table(
      variant = "wildboot", estimator = "TWFE_wildboot",
      outcome = "bev_neuzulassungen_p100k",
      estimate = coef(m_twfe)["treat"],
      se       = NA_real_,
      ci_lo    = bt$conf_int[1L],
      ci_hi    = bt$conf_int[2L],
      p_value  = bt$p_val
    )
  }
} else {
  cat("fwildclusterboot not installed; skipping wild-cluster bootstrap.\n")
}

robust_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
robust_dt[, ci_lo := fifelse(is.na(ci_lo) & !is.na(se), estimate - 1.96 * se, ci_lo)]
robust_dt[, ci_hi := fifelse(is.na(ci_hi) & !is.na(se), estimate + 1.96 * se, ci_hi)]
fwrite(robust_dt, file.path(out_dir, "est_att_robust.csv"))

# -- Forest figure: ATT × variant (main outcome only) -------------------------

fig_dt <- robust_dt[outcome == "bev_neuzulassungen_p100k"]
if (nrow(fig_dt) > 0L) {
  p <- ggplot(fig_dt,
              aes(x = paste(variant, estimator, sep = "/"),
                  y = estimate)) +
    geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi)) +
    coord_flip() +
    labs(x = NULL, y = "ATT (BEV new registrations p100k)",
         caption = "Robustness grid — broad frame, main outcome.") +
    theme_minimal(base_size = 11)
  ggsave(file.path(out_dir, "fig_robust_grid.pdf"), p,
         width = 7, height = 5.5)
}

cat(sprintf("\nDiD robustness outputs -> %s\n", out_dir))
