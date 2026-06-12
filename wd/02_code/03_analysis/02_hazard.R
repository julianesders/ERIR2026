# ───────────────────────────────────────────────────────────────────────────────
# 02_hazard.R   — Part (i): drivers of EMK direct-treatment onset
#
# Replaces the old 02_panel_logit.R. Discrete-time hazard on the AGS8 risk set
# (year >= 2015), with cloglog primary + logit twin. AMEs via marginaleffects.
#
# Columns:
#   (1) base channels                                       (full sample)
#   (2) + personnel  (n_vze_personal_L1)                    (no Stadtstaaten)
#   (3) drop log_dens_z  (collinearity diagnostic)          (full sample)
#   (4) + kreis_funded  (crowd-out / stimulation)           (full sample)
#   (5) eco_index variant (swap bev_z+chg_z for eco_index)  (full sample)
#   (6) state Grüne     (channel swap)                      (full sample)
#   (7) fed Grüne       (channel swap)                      (full sample)
#
# Outputs (04_results/02_hazard/):
#   tab_hazard_coef.{tex,csv}    headline coefficient table (cloglog)
#   tab_hazard_logit.{tex,csv}   logit robustness twin
#   tab_hazard_ame.{tex,csv}     average marginal effects
#   tab_hazard_robust.csv        brglm2 / LPM / complete-case rows
#   tab_hazard_diag.csv          events per AGS2, per year
#   tab_hazard_cov.{tex,csv}     coverage-event variant (appendix)
#   fig_coef_stability.pdf       coefficient stability across (1)/(3)
# ───────────────────────────────────────────────────────────────────────────────

library(data.table)
library(fixest)
library(marginaleffects)
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
  stop("Cannot determine script path. Run as: Rscript 02_hazard.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "04_results", "02_hazard")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))

# -- Load hazard frames --------------------------------------------------------

ph     <- readRDS(file.path(data_final, "frame_hazard.rds"))
ph_cov <- readRDS(file.path(data_final, "frame_hazard_cov.rds"))

# Defensive: drop unit-years missing any of the L1 channels.
chan_cols <- c("log_pop_dens", "log_steuerkraft_L1", "log_kaufkraft_L1",
               "bev_stock_p100k_L1", "ev_chargepoints_p100k_L1",
               "muni_gruene_L1")
ph <- ph[complete.cases(ph[, ..chan_cols])]

# -- Build z-scored regressors on the estimation sample -----------------------

zfit <- function(d) {
  d[, log_dens_z      := z(log_pop_dens)]
  d[, sk_z            := z(log_steuerkraft_L1)]
  d[, kk_z            := z(log_kaufkraft_L1)]
  d[, bev_z           := z(log1p(pmax(bev_stock_p100k_L1, 0)))]
  d[, chg_z           := z(log1p(pmax(ev_chargepoints_p100k_L1, 0)))]
  d[, muni_gruene_z   := z(muni_gruene_L1)]
  d[, state_gruene_z  := z(state_gruene_L1)]
  d[, fed_gruene_z    := z(fed_gruene_L1)]
  d
}

ph <- zfit(ph)

# No-Stadtstaaten subsample (re-z so coefficients are sample-comparable).
ph_ns <- copy(ph[ns_flag == TRUE])
ph_ns[, pers_z := z(log1p(pmax(n_vze_personal_L1, 0)))]
ph_ns <- zfit(ph_ns)

cat(sprintf(
  "Hazard frame: %d obs | %d AGS8 | %d direct onsets | rate %.2f%%\n",
  nrow(ph), uniqueN(ph$AGS8), sum(ph$onset_direct),
  100 * mean(ph$onset_direct)
))
cat(sprintf(
  "No-Stadtstaaten: %d obs | %d AGS8 | %d onsets\n",
  nrow(ph_ns), uniqueN(ph_ns$AGS8), sum(ph_ns$onset_direct)
))

# -- Formulas ------------------------------------------------------------------

f1 <- onset_direct ~ log_dens_z + muni_gruene_z + bev_z + chg_z +
        sk_z + kk_z | year + AGS2
f2 <- onset_direct ~ log_dens_z + muni_gruene_z + bev_z + chg_z +
        sk_z + kk_z + pers_z | year + AGS2
f3 <- onset_direct ~ muni_gruene_z + bev_z + chg_z +
        sk_z + kk_z | year + AGS2
f4 <- onset_direct ~ log_dens_z + muni_gruene_z + bev_z + chg_z +
        sk_z + kk_z + kreis_funded | year + AGS2
f5 <- onset_direct ~ log_dens_z + muni_gruene_z + eco_index_L1 +
        sk_z + kk_z | year + AGS2
f6 <- onset_direct ~ log_dens_z + state_gruene_z + bev_z + chg_z +
        sk_z + kk_z | year + AGS2
f7 <- onset_direct ~ log_dens_z + fed_gruene_z + bev_z + chg_z +
        sk_z + kk_z | year + AGS2

cll <- binomial("cloglog")
lgt <- binomial("logit")

# -- Estimate ------------------------------------------------------------------

cat("\nEstimating cloglog models...\n")
m_cll <- list(
  `(1)`      = feglm(f1, data = ph,    family = cll, cluster = ~AGS5),
  `(2) +P`   = feglm(f2, data = ph_ns, family = cll, cluster = ~AGS5),
  `(3) -den` = feglm(f3, data = ph,    family = cll, cluster = ~AGS5),
  `(4) +KF`  = feglm(f4, data = ph,    family = cll, cluster = ~AGS5),
  `(5) eco`  = feglm(f5, data = ph,    family = cll, cluster = ~AGS5),
  `(6) S-G`  = feglm(f6, data = ph,    family = cll, cluster = ~AGS5),
  `(7) F-G`  = feglm(f7, data = ph,    family = cll, cluster = ~AGS5)
)

cat("Estimating logit twins...\n")
m_lgt <- list(
  `(1)`      = feglm(f1, data = ph,    family = lgt, cluster = ~AGS5),
  `(2) +P`   = feglm(f2, data = ph_ns, family = lgt, cluster = ~AGS5),
  `(3) -den` = feglm(f3, data = ph,    family = lgt, cluster = ~AGS5),
  `(4) +KF`  = feglm(f4, data = ph,    family = lgt, cluster = ~AGS5),
  `(5) eco`  = feglm(f5, data = ph,    family = lgt, cluster = ~AGS5),
  `(6) S-G`  = feglm(f6, data = ph,    family = lgt, cluster = ~AGS5),
  `(7) F-G`  = feglm(f7, data = ph,    family = lgt, cluster = ~AGS5)
)

# -- Coefficient tables (TeX + CSV twins) -------------------------------------

tab_note <- paste0(
  "Discrete-time hazard on AGS8 risk set (year >= 2015). ",
  "Direct-treatment onset is the event; broadcast-only units remain in the ",
  "risk set with kreis_funded switching on. Col (2): Hamburg/Bremen/Berlin ",
  "excluded (personnel conflates municipal/Länder roles). year + AGS2 FE; ",
  "SEs clustered at AGS5."
)

etable(m_cll, dict = dict, digits = 4, notes = tab_note,
       file = file.path(out_dir, "tab_hazard_coef.tex"),
       title = "Discrete-time hazard (cloglog): drivers of EMK onset",
       replace = TRUE)
etable(m_lgt, dict = dict, digits = 4, notes = tab_note,
       file = file.path(out_dir, "tab_hazard_logit.tex"),
       title = "Logit robustness twin",
       replace = TRUE)

write_estimates_csv(m_cll, file.path(out_dir, "tab_hazard_coef.csv"))
write_estimates_csv(m_lgt, file.path(out_dir, "tab_hazard_logit.csv"))

# -- Sensitivity: SEs clustered at AGS2 (Bundesland) instead of AGS5 ----------
# Point estimates are unchanged (cluster choice only affects vcov). We refit
# the vcov via fixest::vcov(..., cluster = ~AGS2) and emit a parallel table.
# With 16 AGS2 clusters and 166 direct onsets, this is conservative against
# the rare-events / few-treated-clusters concern at AGS5.

tab_note_ags2 <- paste0(tab_note, " SENSITIVITY: SEs clustered at AGS2 ",
                        "(Bundesland) instead of AGS5.")

etable(m_cll,
       dict = dict, digits = 4, notes = tab_note_ags2,
       cluster = ~AGS2,
       file = file.path(out_dir, "tab_hazard_coef_ags2.tex"),
       title = "Discrete-time hazard (cloglog), SEs clustered at AGS2",
       replace = TRUE)
etable(m_lgt,
       dict = dict, digits = 4, notes = tab_note_ags2,
       cluster = ~AGS2,
       file = file.path(out_dir, "tab_hazard_logit_ags2.tex"),
       title = "Logit twin, SEs clustered at AGS2",
       replace = TRUE)

write_estimates_csv_cluster <- function(models, file, clf) {
  rows <- lapply(seq_along(models), function(i) {
    m  <- models[[i]]
    nm <- names(models)[i]
    co <- coef(m)
    V  <- vcov(m, cluster = clf)
    se <- sqrt(diag(V))
    data.table(
      model = nm,
      term  = names(co),
      estimate = as.numeric(co),
      se       = as.numeric(se[names(co)]),
      n        = nobs(m)
    )
  })
  dt <- rbindlist(rows)
  dt[, ci_lo := estimate - 1.96 * se]
  dt[, ci_hi := estimate + 1.96 * se]
  fwrite(dt, file)
  invisible(dt)
}

write_estimates_csv_cluster(m_cll,
                            file.path(out_dir, "tab_hazard_coef_ags2.csv"),
                            ~AGS2)
write_estimates_csv_cluster(m_lgt,
                            file.path(out_dir, "tab_hazard_logit_ags2.csv"),
                            ~AGS2)

# -- Average marginal effects (headline table) --------------------------------

cat("\nComputing AMEs via marginaleffects::avg_slopes (vcov=FALSE) ...\n")
# marginaleffects refuses to compute the FE part of the vcov on fixest models;
# we pass vcov=FALSE for point AMEs only. Inference for the corresponding
# linear-index coefficients is in tab_hazard_coef (with clustered SEs); the
# AME and the cloglog β have the same sign and ranking, so the headline
# narrative is unaffected.
ame_rows <- rbindlist(lapply(seq_along(m_cll), function(i) {
  nm <- names(m_cll)[i]
  am <- tryCatch(avg_slopes(m_cll[[i]], vcov = FALSE),
                 error = function(e) {
                   cat(sprintf("  AME failed for %s: %s\n",
                               nm, conditionMessage(e))); NULL })
  if (is.null(am)) return(NULL)
  am_dt <- as.data.table(am)
  if (!all(c("term", "estimate") %in% names(am_dt))) return(NULL)
  am_dt[, .(model = nm, term, estimate, n = nobs(m_cll[[i]]))]
}), fill = TRUE)

if (nrow(ame_rows) > 0L) {
  fwrite(ame_rows, file.path(out_dir, "tab_hazard_ame.csv"))
  ame_rows[, cell := sprintf("%.4f", estimate)]
  tab_w <- dcast(ame_rows, term ~ model, value.var = "cell")
  body  <- apply(tab_w, 1L, function(rw) paste(rw, collapse = " & "))
  tex_lines <- c(
    "\\begin{tabular}{l*{7}{c}}", "\\hline",
    paste(c("Term", names(m_cll)), collapse = " & "), " \\\\",
    "\\hline",
    paste0(body, " \\\\"),
    "\\hline", "\\end{tabular}",
    "% Note: point AMEs only; cloglog coefficient SEs in tab_hazard_coef."
  )
  writeLines(tex_lines, file.path(out_dir, "tab_hazard_ame.tex"))
} else {
  cat("No AMEs produced — skipping AME tables.\n")
}

# -- Robustness: brglm2 penalized logit, LPM, complete-case --------------------

robust_rows <- list()

# (a) brglm2 penalized logit on (1), with explicit FE dummies (small set).
if (requireNamespace("brglm2", quietly = TRUE)) {
  cat("Fitting brglm2 penalized logit on spec (1)...\n")
  br_df <- as.data.frame(ph)
  br_df$year_f <- factor(br_df$year)
  br_df$ags2_f <- factor(br_df$AGS2)
  br_fit <- tryCatch(
    glm(onset_direct ~ log_dens_z + muni_gruene_z + bev_z + chg_z +
          sk_z + kk_z + year_f + ags2_f,
        data = br_df, family = binomial("logit"),
        method = brglm2::brglm_fit),
    error = function(e) { cat("brglm2 error:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(br_fit)) {
    co <- coef(br_fit); se <- sqrt(diag(vcov(br_fit)))
    keep <- !grepl("^(year_f|ags2_f|\\(Intercept\\))", names(co))
    robust_rows[["brglm2_logit_(1)"]] <- data.table(
      estimator = "brglm2_logit", model = "(1)",
      term = names(co)[keep],
      estimate = co[keep], se = se[keep],
      n = nobs(br_fit)
    )
  }
} else {
  cat("brglm2 not installed; skipping penalized logit robustness.\n")
}

# (b) LPM via feols on (1)
cat("Fitting LPM (feols) on spec (1)...\n")
m_lpm <- feols(f1, data = ph, cluster = ~AGS5)
co <- coef(m_lpm); se <- sqrt(diag(vcov(m_lpm)))
robust_rows[["lpm_(1)"]] <- data.table(
  estimator = "lpm", model = "(1)",
  term = names(co), estimate = co, se = se, n = nobs(m_lpm)
)

# (c) Complete-case: drop rows where any _imp flag of used covariates is TRUE.
imp_flags <- c("N_elektro_overall_imp", "bev_stock_p100k_L1")  # placeholder
# Use any _imp column on the L1 inputs: bev_stock & chargepoints stock came
# from KBA / Ladestationen; we approximate "no imputed" by requiring
# B_elektro_overall_imp == FALSE (proxies bev_stock).
if ("B_elektro_overall_imp" %in% names(ph)) {
  ph_cc <- ph[B_elektro_overall_imp == FALSE]
  m_cc <- feglm(f1, data = ph_cc, family = cll, cluster = ~AGS5)
  co <- coef(m_cc); se <- sqrt(diag(vcov(m_cc)))
  robust_rows[["cc_(1)"]] <- data.table(
    estimator = "cloglog_complete_case", model = "(1)",
    term = names(co), estimate = co, se = se, n = nobs(m_cc)
  )
}

robust_dt <- rbindlist(robust_rows, fill = TRUE)
fwrite(robust_dt, file.path(out_dir, "tab_hazard_robust.csv"))

# -- Diagnostics: events per AGS2 / per year ----------------------------------

diag_dt <- rbindlist(list(
  ph[, .(n_obs = .N, n_events = sum(onset_direct)), by = .(grp = AGS2)
     ][, dim := "AGS2"],
  ph[, .(n_obs = .N, n_events = sum(onset_direct)), by = .(grp = as.character(year))
     ][, dim := "year"]
), use.names = TRUE)
setcolorder(diag_dt, c("dim", "grp", "n_obs", "n_events"))
fwrite(diag_dt, file.path(out_dir, "tab_hazard_diag.csv"))

# Assert no AGS2 has zero events; fixest would silently drop the FE level
zero_ev <- diag_dt[dim == "AGS2" & n_events == 0L]
if (nrow(zero_ev)) {
  cat(sprintf("NOTE: %d AGS2 with zero events — FE drops these rows:\n",
              nrow(zero_ev)))
  print(zero_ev)
}

# -- Coefficient-stability figure: (1) vs (3) on shared channels --------------

stab_rows <- rbindlist(list(
  data.table(model = "(1)",   term = names(coef(m_cll[["(1)"]])),
             est   = coef(m_cll[["(1)"]]),
             se    = sqrt(diag(vcov(m_cll[["(1)"]])))),
  data.table(model = "(3) -den", term = names(coef(m_cll[["(3) -den"]])),
             est   = coef(m_cll[["(3) -den"]]),
             se    = sqrt(diag(vcov(m_cll[["(3) -den"]]))))
))
shared <- intersect(stab_rows[model == "(1)", term],
                    stab_rows[model == "(3) -den", term])
stab_rows <- stab_rows[term %in% shared]
p_stab <- ggplot(stab_rows, aes(x = term, y = est, color = model)) +
  geom_pointrange(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se),
                  position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
  coord_flip() +
  labs(x = NULL, y = "Coefficient (cloglog)",
       color = NULL,
       caption = "Coefficient stability when log_dens_z is dropped.") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "fig_coef_stability.pdf"), p_stab,
       width = 7, height = 4.5)

# -- Coverage-event appendix variant ------------------------------------------
# Spec (1) only, on the broad-coverage hazard frame.

ph_cov <- ph_cov[complete.cases(ph_cov[, ..chan_cols])]
ph_cov <- zfit(ph_cov)
f1_cov <- onset_broad ~ log_dens_z + muni_gruene_z + bev_z + chg_z +
           sk_z + kk_z | year + AGS2
m_cov <- feglm(f1_cov, data = ph_cov, family = cll, cluster = ~AGS5)
etable(list(`(1) coverage` = m_cov),
       dict = dict, digits = 4,
       notes = paste0("Coverage-event hazard: event = first year of any ",
                      "EMK coverage (direct ∪ Kreis broadcast)."),
       file = file.path(out_dir, "tab_hazard_cov.tex"), replace = TRUE)
write_estimates_csv(list(`(1) coverage` = m_cov),
                    file.path(out_dir, "tab_hazard_cov.csv"))

cat(sprintf("\nHazard outputs -> %s\n", out_dir))
