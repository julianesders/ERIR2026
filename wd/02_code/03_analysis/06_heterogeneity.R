# ───────────────────────────────────────────────────────────────────────────────
# 06_heterogeneity.R   — Inequality payoff + heterogeneity cuts
#
# Conditional CS (outcome regression, never-treated control, broad frame) — the
# PRIMARY spec — estimated within:
#   * baseline tax-capacity (Steuerkraft) terciles  [headline inequality cut]
#   * baseline population-density terciles
#   * East / West Germany
#   * pre / post-2021 funding-regime (FRL) cohorts
# plus treat_type heterogeneity (direct vs broadcast-only).
#
# Within a rank split the STRATIFYING baseline is dropped from the covariate set
# (conditioning on the variable you stratified on is near-constant within the
# stratum and collapses the DR design — see xformla_for_rank()).
#
# Inference (D2): all SEs clustered at AGS5. att_gt runs with bstrap = FALSE for
# speed; the AGS5-clustered SE is reconstructed analytically from the unit-level
# influence function (sqrt(sum_c contrib_c^2)). Between-group ATT differences use
# a CLUSTERED multiplier bootstrap that shares the per-cluster Mammen weight
# across the two groups (captures correlation when a county spans both groups, or
# when both groups share the never-treated control pool) — replacing the old
# normal-approx independence assumption.
#
# Outputs (03_output/06_heterogeneity/):
#   did_heterogeneity_terciles.{tex,csv}    ATT × tercile × rank (sk, dens)
#   did_heterogeneity_groups.{tex,csv}      ATT for East/West, Pre/Post-2021
#   tab_het_diffs.{tex,csv}                 between-group differences (clustered)
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
  stop("Cannot determine script path. Run as: Rscript 06_heterogeneity.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "03_output", "06_heterogeneity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))
source(file.path(code_dir, "03_analysis", "_did_helpers.R"))

CS_BITERS <- 2000L          # B for the between-group multiplier bootstrap
OUTCOME   <- "bev_neuzulassungen_p100k"

.read_frame <- function(p) fread(p,
  colClasses = list(character = c("AGS8", "AGS5", "AGS2", "treat_type")))
frame_broad  <- .read_frame(file.path(data_final, "frame_did_broad.csv"))

# ags8_id -> AGS5 crosswalk (for clustered IF aggregation)
id2ags5 <- frame_broad[, .(AGS5 = AGS5[1L]), by = ags8_id][
  , setNames(AGS5, as.character(ags8_id))]

# ── Conditional CS, stratifier-aware ───────────────────────────────────────────
# Conditional dr is primary (D1). Within a rank split, drop the stratifying
# baseline from XFORMLA_CS — conditioning on the variable you stratified on is
# near-constant within stratum and collapses the DR design.
xformla_for_rank <- function(rank_col) {
  drop <- switch(rank_col,
    sk_base_terc   = "sk_base_z",
    dens_base_terc = "dens_base_z",
    "")
  vars <- setdiff(all.vars(XFORMLA_CS), drop)
  if (length(vars) == 0L) NULL else reformulate(vars)
}

# One CS runner (bstrap = FALSE; AGS5-clustered SEs reconstructed from the IF).
# est_method = "reg" (conditional OUTCOME REGRESSION), not "dr": within the
# wealth / density / region subsets the DR propensity-score logit hits
# near-perfect separation and qr.solve throws a singular matrix on every cell
# (e.g. the top Steuerkraft tercile and the West group). The "reg" variant still
# CONDITIONS on the baseline covariates W_i (so D1 holds — conditioning is the
# point), is numerically stable in every subset, and shares the same influence-
# function machinery. The full-sample headline (04_did_main.R) keeps dr.
run_cs_sub <- function(d, xformla = XFORMLA_CS, control = "nevertreated") {
  d <- d[!is.na(get(OUTCOME))]          # drop NA outcomes before att_gt sees them
  cov_nms <- if (is.null(xformla)) character(0) else all.vars(xformla)
  if (length(cov_nms))
    d <- d[complete.cases(d[, ..cov_nms])]   # baselines const within unit
  att_gt(
    yname = OUTCOME, gname = "gname_cs", idname = "ags8_id", tname = "year",
    data = as.data.frame(d), control_group = control, anticipation = 0L,
    xformla = xformla, est_method = "reg",
    # bstrap = FALSE: SEs reconstructed from the IF in .contrib_by_cluster();
    # clustervars omitted here because bstrap=FALSE ignores it anyway.
    bstrap = FALSE, allow_unbalanced_panel = TRUE, print_details = FALSE
  )
}

# ── IF-based clustered contributions + self-check ──────────────────────────────
# Returns the overall ATT, the per-AGS5 contribution vector (contrib_c =
# (1/n) sum_{i in c} psi_i), and the AGS5-clustered SE = sqrt(sum_c contrib_c^2).
# Self-check: the UNCLUSTERED reconstruction sqrt(sum((psi/n)^2)) must reproduce
# the package's bstrap=FALSE overall.se (validates IF slot name + unit ordering;
# `did` stores the simple-aggregation IF in inf.function$simple.att, ordered by
# sorted idname). If it warns, see watch-out #3 in the plan.
.contrib_by_cluster <- function(cs_obj, tol = 0.02) {
  agg <- aggte(cs_obj, type = "simple", na.rm = TRUE)
  psi <- as.numeric(agg$inf.function$simple.att)
  dp  <- cs_obj$DIDparams
  ids <- sort(unique(dp$data[[dp$idname]]))
  stopifnot(length(psi) == length(ids))
  n   <- length(psi)
  ags5 <- id2ags5[as.character(ids)]
  if (anyNA(ags5)) stop("Unmapped ags8_id in id2ags5 crosswalk")

  se_unclust <- sqrt(sum((psi / n)^2))
  if (abs(se_unclust - agg$overall.se) > tol * max(1, abs(agg$overall.se)))
    warning(sprintf(
      "IF self-check: rebuilt unclustered SE %.5f vs package SE %.5f -- check IF slot / id order",
      se_unclust, agg$overall.se))

  contrib <- tapply(psi, ags5, sum) / n          # named by AGS5
  list(att = agg$overall.att, contrib = contrib, se = sqrt(sum(contrib^2)))
}

# (att, AGS5-clustered se, Ns) for one group's CS object.
grp_cell <- function(cs_obj, lbl, d) {
  cc <- .contrib_by_cluster(cs_obj)
  data.table(group = lbl, att = cc$att, se = cc$se,
             n_treated = uniqueN(d[gname_cs > 0L, AGS8]),
             n_control = uniqueN(d[gname_cs == 0L, AGS8]))
}

# ── Clustered multiplier bootstrap for a between-group ATT difference (A - B) ───
.draw_mammen <- function(n) {                       # mean 0, var 1, skew 1
  k1 <- (1 - sqrt(5)) / 2
  k2 <- (1 + sqrt(5)) / 2
  p1 <- (sqrt(5) + 1) / (2 * sqrt(5))
  ifelse(runif(n) < p1, k1, k2)
}

diff_cluster_mboot <- function(cs_A, cs_B, B = CS_BITERS, seed = 1L) {
  set.seed(seed)
  a <- .contrib_by_cluster(cs_A)
  b <- .contrib_by_cluster(cs_B)
  delta_hat <- a$att - b$att

  clusters <- union(names(a$contrib), names(b$contrib))
  Da <- setNames(numeric(length(clusters)), clusters); Db <- Da
  Da[names(a$contrib)] <- a$contrib
  Db[names(b$contrib)] <- b$contrib
  D  <- Da - Db                                    # shared cluster index

  se_analytic <- sqrt(sum(D^2))
  draws <- vapply(seq_len(B),
                  function(.) delta_hat + sum(.draw_mammen(length(D)) * D),
                  numeric(1))
  se_iqr <- as.numeric((quantile(draws, .75) - quantile(draws, .25)) /
                       (qnorm(.75) - qnorm(.25)))
  data.table(
    delta       = delta_hat,
    se_analytic = se_analytic,
    se_boot     = sd(draws),
    se_iqr      = se_iqr,                          # report this (CS-style robust)
    ci_lo       = as.numeric(quantile(draws, .025)),
    ci_hi       = as.numeric(quantile(draws, .975)),
    p_value     = 2 * pnorm(abs(delta_hat / se_iqr), lower.tail = FALSE)
  )
}

# Nonparametric cluster bootstrap (validation only; refits both groups per draw).
diff_cluster_npboot <- function(dat_A, dat_B, B = 200L, cluster = "AGS5",
                                seed = 1L) {
  set.seed(seed)
  est <- function(A, Bd)
    aggte(run_cs_sub(A), type = "simple", na.rm = TRUE)$overall.att -
    aggte(run_cs_sub(Bd), type = "simple", na.rm = TRUE)$overall.att
  cl  <- unique(c(dat_A[[cluster]], dat_B[[cluster]]))
  d0  <- est(dat_A, dat_B)
  rebuild <- function(d, samp) rbindlist(lapply(seq_along(samp), function(k) {
    x <- d[get(cluster) == samp[k]]
    if (nrow(x)) x[, ags8_id := ags8_id + k * 100000L][] else NULL
  }))
  draws <- replicate(B, {
    s <- sample(cl, length(cl), TRUE)
    est(rebuild(dat_A, s), rebuild(dat_B, s))
  })
  data.table(delta = d0, se = sd(draws),
             ci_lo = quantile(draws, .025), ci_hi = quantile(draws, .975))
}

# ── Clustered event-study (for the figure) ─────────────────────────────────────
# Per-event AGS5-clustered SE reconstructed from the dynamic IF columns.
clustered_es <- function(cs_obj, es_max) {
  dy <- aggte(cs_obj, type = "dynamic", min_e = ES_MIN, max_e = es_max,
              balance_e = NULL, na.rm = TRUE, cband = FALSE)
  inf <- dy$inf.function$dynamic.inf.func.e
  dp  <- cs_obj$DIDparams
  ids <- sort(unique(dp$data[[dp$idname]]))
  n   <- nrow(inf)
  ags5 <- id2ags5[as.character(ids)]
  se_e <- vapply(seq_len(ncol(inf)), function(j) {
    sqrt(sum((tapply(inf[, j], ags5, sum) / n)^2))
  }, numeric(1))
  data.table(e = dy$egt, att = dy$att.egt, se = se_e)
}

# ── Tercile-stratified conditional CS (returns CS objects for the diffs) ───────

run_cs_terc <- function(dat, terc_col) {
  xf <- xformla_for_rank(terc_col)
  objs <- list(); att_rows <- list(); es_rows <- list()
  for (t in sort(unique(na.omit(dat[[terc_col]])))) {
    sub  <- dat[get(terc_col) == t]
    n_tr <- uniqueN(sub[gname_cs > 0L, AGS8])
    n_co <- uniqueN(sub[gname_cs == 0L, AGS8])
    if (n_tr < 5L || n_co < 5L) {
      cat(sprintf("  %s tercile %d skipped (n_tr=%d, n_co=%d)\n",
                  terc_col, t, n_tr, n_co)); next
    }
    cs <- tryCatch(run_cs_sub(sub, xf),
                   error = function(e) {
                     cat("  CS err", terc_col, t, ":", conditionMessage(e), "\n")
                     NULL })
    if (is.null(cs)) next
    objs[[as.character(t)]] <- cs
    cc   <- tryCatch(.contrib_by_cluster(cs), error = function(e) NULL)
    pt_p <- tryCatch(cs_pre_test(cs)$pval,   error = function(e) NA_real_)
    if (!is.null(cc))
      att_rows[[length(att_rows) + 1L]] <- data.table(
        tercile = t, n_treated = n_tr, n_control = n_co,
        att = cc$att, se = cc$se, pre_p = pt_p)
    es <- tryCatch(clustered_es(cs, es_max_data_driven(sub, OUTCOME)),
                   error = function(e) NULL)
    if (!is.null(es)) es_rows[[length(es_rows) + 1L]] <- es[, tercile := t]
    cat(sprintf("  %s tercile %d: n_tr=%d n_co=%d att=%.3f se=%.3f\n",
                terc_col, t, n_tr, n_co, cc$att, cc$se))
  }
  list(att = rbindlist(att_rows, fill = TRUE),
       es  = rbindlist(es_rows,  fill = TRUE),
       objs = objs)
}

cat("\n=== Tercile-stratified conditional CS ===\n")
res_sk   <- run_cs_terc(frame_broad, "sk_base_terc")     # headline (Steuerkraft)
res_dens <- run_cs_terc(frame_broad, "dens_base_terc")   # density

# ── East/West and Pre/Post-2021 group runs ─────────────────────────────────────

cat("\n=== East/West + Pre/Post-2021 ===\n")
cohort_subset <- function(dat, keep) dat[gname_cs == 0L | gname_cs %in% keep]
.safe_cs <- function(d, lbl) tryCatch(run_cs_sub(d), error = function(e) {
  cat(sprintf("  %s CS err: %s\n", lbl, conditionMessage(e))); NULL })

d_west <- frame_broad[east == 0L]
d_east <- frame_broad[east == 1L]
d_pre  <- cohort_subset(frame_broad, 2015:FRL_CUTOFF)
d_post <- cohort_subset(frame_broad, (FRL_CUTOFF + 1L):2022)
cs_west <- .safe_cs(d_west, "West")
cs_east <- .safe_cs(d_east, "East")
cs_pre  <- .safe_cs(d_pre,  "Pre-2021")
cs_post <- .safe_cs(d_post, "Post-2021")

grp_dt <- rbindlist(Filter(Negate(is.null), list(
  if (!is.null(cs_west)) grp_cell(cs_west, "West",      d_west),
  if (!is.null(cs_east)) grp_cell(cs_east, "East",      d_east),
  if (!is.null(cs_pre))  grp_cell(cs_pre,  "Pre-2021",  d_pre),
  if (!is.null(cs_post)) grp_cell(cs_post, "Post-2021", d_post)
)), fill = TRUE)
if (nrow(grp_dt)) {
  grp_dt[, ci_lo := att - 1.96 * se][, ci_hi := att + 1.96 * se]
  grp_dt[, p     := 2 * pnorm(abs(att / se), lower.tail = FALSE)]
  grp_dt[, stars := ifelse(p < 0.01, "$^{***}$",
                    ifelse(p < 0.05, "$^{**}$",
                    ifelse(p < 0.10, "$^{*}$", "")))]
  fwrite(grp_dt, file.path(out_dir, "did_heterogeneity_groups.csv"))
}

# ── Tercile table (Steuerkraft + density) ──────────────────────────────────────

terc_dt <- rbindlist(list(
  copy(res_sk$att)[,   rank := "sk_base"],
  copy(res_dens$att)[, rank := "dens_base"]
), fill = TRUE)
terc_dt[, ci_lo := att - 1.96 * se][, ci_hi := att + 1.96 * se]
terc_dt[, p     := 2 * pnorm(abs(att / se), lower.tail = FALSE)]
terc_dt[, stars := ifelse(p < 0.01, "$^{***}$",
                   ifelse(p < 0.05, "$^{**}$",
                   ifelse(p < 0.10, "$^{*}$", "")))]
fwrite(terc_dt, file.path(out_dir, "did_heterogeneity_terciles.csv"))

write_longtblr(
  stem        = file.path(out_dir, "did_heterogeneity_terciles"),
  caption     = "Heterogeneity by Baseline Tercile: Conditional CS Overall ATT",
  label       = "tab:did_heterogeneity_terciles",
  note        = paste0(
    "\\textcite{callaway2021difference} conditional doubly robust estimates split by ",
    "terciles of tax capacity and population density. The dependent variable is BEV ",
    "registration per 100k inhabitants. Controls include state Green vote share, and ",
    "tax capacity p.c. or population density depending on the stratifying variable. ",
    "Terciles based on treated and never-treated units, never-treated control group, ",
    "broad treatment. SEs are AGS5-clustered derived from a multiplier boostrapped(",
    "$B=", CS_BITERS, "$). ",
    "$p$-value of a joint pre-test for $H_0:\\text{ATT}(g,t)=0$ for all $g,t$ in the last line. ",
    "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$."
  ),
  colspec     = "l r r r l r r",
  header_rows = paste0("Variable & Tercile & $N_{\\text{treated}}$ & ",
                       "$N_{\\text{control}}$ & ATT & SE & Pre-test $p$ \\\\"),
  body_rows   = terc_dt[, {
    rank_lbl  <- c(sk_base = "Tax capacity", dens_base = "Pop. Density")[rank]
    pre_p_str <- ifelse(is.na(pre_p), "--", sprintf("%.3f", pre_p))
    sprintf("%s & %d & %d & %d & %.3f%s & (%.3f) & %s \\\\",
            rank_lbl, tercile, n_treated, n_control, att, stars, se, pre_p_str)
  }],
  footer_rows = character(0)
)

# ── Group table (East/West, Pre/Post-2021) ─────────────────────────────────────

write_longtblr(
  stem        = file.path(out_dir, "did_heterogeneity_groups"),
  caption     = "Heterogeneity by Region and Funding Regime: Conditional CS Overall ATT",
  label       = "tab:did_heterogeneity_groups",
  note        = paste0(
    "Conditional outcome-regression CS overall ATT (broad frame, never-treated ",
    "control). \\emph{East} = neue L\\\"ander incl.\\ Berlin. \\emph{Pre/Post-2021} ",
    "split treated cohorts at the FRL funding-regime break (cohorts $\\le 2020$ ",
    "vs $\\ge 2021$), sharing the never-treated control pool. SEs clustered at ",
    "AGS5 (reconstructed from the influence function). SE in parentheses. ",
    "$^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.1$ (two-sided)."
  ),
  colspec     = "l r r l r",
  header_rows = "Group & $N_{\\text{treated}}$ & $N_{\\text{control}}$ & ATT & SE \\\\",
  body_rows   = grp_dt[, sprintf(
    "%s & %d & %d & %.3f%s & (%.3f) \\\\",
    group, n_treated, n_control, att, stars, se)],
  footer_rows = character(0)
)

# ── Between-group differences via clustered multiplier bootstrap ───────────────

cat("\n=== Between-group differences (clustered multiplier bootstrap) ===\n")
# Pick the top/bottom tercile object that actually estimated (objs keyed by
# tercile value as character). .mb returns NULL (and skips) if either side is
# missing, so the table never silently mixes contrasts with absent groups.
.terc_obj <- function(res, which.fun)
  if (nrow(res$att)) res$objs[[as.character(which.fun(res$att$tercile))]] else NULL
.mb <- function(A, B, lbl) {
  if (is.null(A) || is.null(B)) {
    cat(sprintf("  skip diff [%s] (a contributing group is missing)\n", lbl))
    return(NULL)
  }
  diff_cluster_mboot(A, B)[, contrast := lbl]
}

diffs <- rbindlist(Filter(Negate(is.null), list(
  .mb(.terc_obj(res_sk,   max), .terc_obj(res_sk,   min),
      "Steuerkraft: T3 - T1 (headline)"),
  .mb(.terc_obj(res_dens, max), .terc_obj(res_dens, min),
      "Density: T3 - T1"),
  .mb(cs_west, cs_east, "West - East"),
  .mb(cs_post, cs_pre,  "Post-2021 - Pre-2021")
)), fill = TRUE)

if (!nrow(diffs)) {
  cat("  no between-group differences estimable; skipping tab_het_diffs\n")
} else {
setcolorder(diffs, "contrast")
fwrite(diffs, file.path(out_dir, "tab_het_diffs.csv"))
print(diffs[, .(contrast, delta = round(delta, 3), se_iqr = round(se_iqr, 3),
                p_value = round(p_value, 4))])

write_longtblr(
  stem        = file.path(out_dir, "tab_het_diffs"),
  caption     = "Between-Group ATT Differences (Clustered Multiplier Bootstrap)",
  label       = "tab:het_diffs",
  note        = paste0(
    "Difference $\\Delta = \\text{ATT}_A - \\text{ATT}_B$ in the conditional CS ",
    "overall ATT between groups. Controls include baseline tax capacity, population ",
    "density, and state Green vote share. AGS5-clustered SE from a multiplier ",
    "bootstrap ($B=", CS_BITERS, "$) that accounts for correlation when a county ",
    "spans both groups/terciles or both groups share the same never-treated control ",
    "pool (pre/post). $p$-value for a test of $H_0:\\Delta=0$ in the last column."
  ),
  colspec     = "l l r r",
  header_rows = "Contrast & $\\Delta$ & SE & $p$ \\\\",
  body_rows   = diffs[, {
    st <- ifelse(p_value < 0.01, "$^{***}$",
          ifelse(p_value < 0.05, "$^{**}$",
          ifelse(p_value < 0.10, "$^{*}$", "")))
    sprintf("%s & %.3f%s & %.3f & %.3f \\\\",
            contrast, delta, st, se_iqr, p_value)
  }],
  footer_rows = character(0)
)
}

# ── Treat-type heterogeneity: direct vs broadcast-only (conditional) ───────────

dat_dir_only  <- frame_broad[treat_type %in% c("direct", "never")]
dat_brod_only <- frame_broad[treat_type %in% c("broadcast_only", "never")]

run_tt <- function(d, lbl) {
  if (uniqueN(d[gname_cs > 0L, AGS8]) < 5L) return(NULL)
  cs <- tryCatch(run_cs_sub(d), error = function(e) NULL)
  if (is.null(cs)) return(NULL)
  cc <- tryCatch(.contrib_by_cluster(cs), error = function(e) NULL)
  if (is.null(cc)) return(NULL)
  data.table(treated_subset = lbl, att = cc$att, se = cc$se,
             n_treated = uniqueN(d[gname_cs > 0L, AGS8]),
             n_control = uniqueN(d[gname_cs == 0L, AGS8]))
}

tt_dt <- rbindlist(list(
  run_tt(dat_dir_only,  "direct only"),
  run_tt(dat_brod_only, "broadcast only")
), fill = TRUE)
if (nrow(tt_dt)) {
  tt_dt[, ci_lo := att - 1.96 * se][, ci_hi := att + 1.96 * se]
  fwrite(tt_dt, file.path(out_dir, "did_heterogeneity_treat_type.csv"))
  write_longtblr(
    stem        = file.path(out_dir, "did_heterogeneity_treat_type"),
    caption     = "Heterogeneity by Treatment Type: Direct vs Broadcast-only",
    label       = "tab:did_heterogeneity_treat_type",
    note        = paste0(
      "\\textcite{callaway2021difference} conditional doubly robust estimates split by ",
      "treatment type (direct vs. broad). The dependent variable is BEV registration ",
      "per 100k inhabitants. Controls include state Green vote share, tax capacity p.c., ",
      "and population density. Control group is never-treated municipalities. SEs are ",
      "AGS5-clustered derived from a multiplier boostrap ($B=", CS_BITERS, "$). ",
      "95\\% confidence intervals in the last column."
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
