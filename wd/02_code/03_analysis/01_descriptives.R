# ───────────────────────────────────────────────────────────────────────────────
# 01_descriptives.R
#
# Inequality storyline: who is COVERED, not how many euros are received.
# Grant sizes barely vary; the distributive question is binary access.
#
# Outputs (04_results/01_descriptives/):
#   fig_concentration_coverage_sk.pdf       Lorenz-style curve, ranked by sk_base
#   fig_concentration_multi.pdf             Multi-ranking Lorenz curves (direct)
#   tab_concentration_index.csv             Concentration indices (broad/direct × ranking)
#   tab_balance.{tex,csv}                   means by treat_type + normalized diffs
#   tab_desc_means.{tex,csv}               Direct vs Never: means + t-test p-values
#   tab_desc_medians.{tex,csv}             Direct vs Never: medians + Mood p-values
#   fig_map_treatment.pdf                   VG250 choropleth of treat_type (optional)
#   tab_corr.csv                            Pairwise correlations on hazard frame
#   tab_vif.{tex,csv}                       VIFs on hazard frame channels
#   fig_map_bev_dual_2023.png              BEV stock p100k + BEV share maps (optional)
# ───────────────────────────────────────────────────────────────────────────────

library(data.table)
library(fixest)
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
  stop("Cannot determine script path. Run as: Rscript 01_descriptives.R")
}
root       <- dirname(dirname(dirname(self)))
code_dir   <- file.path(root, "02_code")
data_final <- file.path(root, "01_data", "03_final")
out_dir    <- file.path(root, "04_results", "01_descriptives")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(code_dir, "03_analysis", "_dict.R"))

# -- Load DiD-broad frame (carries baseline snapshot + treat_type) -------------

.read_frame <- function(p) fread(p,
  colClasses = list(character = c("AGS8", "AGS5", "AGS2", "treat_type")))
dat <- .read_frame(file.path(data_final, "frame_did_broad.csv"))
ph  <- .read_frame(file.path(data_final, "frame_hazard.csv"))

# AGS8-level staffing: AGS5 absolute VZE FTE distributed evenly across the
# Kreis's constituent municipalities (so a Kreis with many AGS8 isn't
# mechanically inflated). `n_vze_personal` is per-100k of AGS5 population,
# so AGS5 absolute FTE = n_vze_personal * (AGS5 pop / 1e5); then divide by
# n_municipalities_in_kreis.
ags5_year_agg <- dat[, .(
  ags5_pop        = sum(xbev, na.rm = TRUE),
  n_muni_in_kreis = uniqueN(AGS8)
), by = .(AGS5, year)]
dat <- merge(dat, ags5_year_agg, by = c("AGS5", "year"), all.x = TRUE)
dat[, staffing_ags8 := (n_vze_personal * ags5_pop / 1e5) / n_muni_in_kreis]

# Cross-section: one row per AGS8 (use the unit's earliest base-window year)
xs <- dat[, .SD[1L], by = AGS8]

# Staffing baseline snapshot per AGS8: earliest 2014-2016 year with non-NA
staffing_base_dt <- dat[year %in% 2014:2016 & !is.na(staffing_ags8),
                        .(AGS8, year, staffing_ags8)]
setorder(staffing_base_dt, AGS8, year)
staffing_base_dt <- staffing_base_dt[, .SD[1L], by = AGS8]
xs <- merge(xs, staffing_base_dt[, .(AGS8, staffing_base = staffing_ags8)],
            by = "AGS8", all.x = TRUE)

# -- 1. Coverage concentration curves ------------------------------------------
# Sort AGS8 ascending by base (kk_base, then sk_base in appendix); cumulative
# share of population vs cumulative share of TREATED population (broad / direct).
# Concentration index = 2 * area between curve and diagonal.

build_lorenz <- function(xs, rank_col, pop_col = "pop_base") {
  d <- xs[!is.na(get(rank_col)) & !is.na(get(pop_col))]
  setorderv(d, rank_col)
  d[, cum_pop := cumsum(get(pop_col)) / sum(get(pop_col))]
  d[, treated_broad  := as.integer(!is.na(first_treat_broad))]
  d[, treated_direct := as.integer(!is.na(first_treat_direct))]
  d[, cum_treat_b := cumsum(get(pop_col) * treated_broad)  / sum(get(pop_col) * treated_broad)]
  d[, cum_treat_d := cumsum(get(pop_col) * treated_direct) / sum(get(pop_col) * treated_direct)]
  rbind(
    data.table(x = d$cum_pop, y = d$cum_treat_b, line = "Broad coverage"),
    data.table(x = d$cum_pop, y = d$cum_treat_d, line = "Direct treatment"),
    data.table(x = c(0, 1),   y = c(0, 1),       line = "45°")
  )
}

# Concentration index: 1 - 2 * trapezoidal area under cumulative curve
# Sign convention: positive => richer-than-poor; negative => pro-poor.
conc_index <- function(d, rank_col, pop_col, treat_col) {
  d <- d[!is.na(get(rank_col)) & !is.na(get(pop_col)) & !is.na(get(treat_col))]
  setorderv(d, rank_col)
  x <- cumsum(d[[pop_col]]) / sum(d[[pop_col]])
  y <- cumsum(d[[pop_col]] * d[[treat_col]]) / sum(d[[pop_col]] * d[[treat_col]])
  # trapezoidal integral of y over x
  area <- sum(diff(x) * (y[-length(y)] + y[-1L]) / 2)
  1 - 2 * area
}

xs[, t_broad  := as.integer(!is.na(first_treat_broad))]
xs[, t_direct := as.integer(!is.na(first_treat_direct))]

# (a) SK Lorenz (kept as appendix; broad + direct overlaid)
ci_sk <- c(
  broad  = conc_index(xs, "sk_base", "pop_base", "t_broad"),
  direct = conc_index(xs, "sk_base", "pop_base", "t_direct")
)
cat("Concentration indices (sk_base):  broad ",
    round(ci_sk["broad"], 3), " | direct ", round(ci_sk["direct"], 3), "\n")

lz_sk <- build_lorenz(xs, "sk_base")

plot_lorenz <- function(lz, rank_label, ci) {
  ggplot(lz, aes(x = x, y = y, color = line, linetype = line)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("Broad coverage" = "#1b6ca8",
                                  "Direct treatment" = "#e67e22",
                                  "45°" = "grey50")) +
    scale_linetype_manual(values = c("Broad coverage" = "solid",
                                     "Direct treatment" = "solid",
                                     "45°" = "dashed")) +
    labs(
      x = paste0("Cumulative population share (Gemeinden ranked by ",
                 rank_label, ")"),
      y = "Cumulative share of treated population",
      color = NULL, linetype = NULL,
      caption = sprintf("Concentration index — broad: %.3f, direct: %.3f",
                        ci["broad"], ci["direct"])
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

ggsave(file.path(out_dir, "fig_concentration_coverage_sk.pdf"),
       plot_lorenz(lz_sk, "baseline log Steuerkraft", ci_sk),
       width = 7, height = 5)

# (b) Multi-ranking broad-coverage Lorenz: one line per ranking variable.
# Each line shows the cumulative share of the broad-treated population
# captured by the bottom-X% of the population when Gemeinden are sorted by
# the corresponding baseline variable.
ranking_vars <- c(
  sk_base       = "Log Steuerkraft",
  kk_base       = "Log Kaufkraft",
  green_base    = "Muni Grüne share",
  staffing_base = "Staffing (FTE)",
  bev_base      = "log1p BEV stock p100k",
  chg_base      = "log1p Charge pts p100k"
)

build_lorenz_one <- function(xs, rank_col, pop_col, treat_col) {
  d <- xs[!is.na(get(rank_col)) & !is.na(get(pop_col))]
  setorderv(d, rank_col)
  d[, cum_pop   := cumsum(get(pop_col)) / sum(get(pop_col))]
  d[, cum_treat := cumsum(get(pop_col) * get(treat_col)) /
                   sum(get(pop_col) * get(treat_col))]
  data.table(x = d$cum_pop, y = d$cum_treat)
}

lz_multi <- rbindlist(lapply(names(ranking_vars), function(rc) {
  d <- build_lorenz_one(xs, rc, "pop_base", "t_direct")
  d[, ranking := ranking_vars[[rc]]][]
}))

ci_multi <- vapply(names(ranking_vars), function(rc)
  conc_index(xs, rc, "pop_base", "t_direct"), numeric(1))
names(ci_multi) <- unname(ranking_vars)
cat("Concentration indices (direct, multi-ranking):\n")
for (nm in names(ci_multi)) {
  cat(sprintf("  %-26s %.3f\n", nm, ci_multi[[nm]]))
}

p_multi <- ggplot(lz_multi, aes(x = x, y = y, color = ranking)) +
  geom_abline(slope = 1, intercept = 0, color = "grey50",
              linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = "Cumulative population share (Gemeinden ranked by variable)",
    y = "Cumulative share of treated population (direct)",
    color = NULL,
    caption = paste0(
      "Concentration indices (direct): ",
      paste(sprintf("%s = %.3f", names(ci_multi), ci_multi),
            collapse = "; ")
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.caption = element_text(size = 8))

ggsave(file.path(out_dir, "fig_concentration_multi.pdf"),
       p_multi, width = 8, height = 5.5)

# CSV: SK both scopes + multi-ranking direct
ci_dt <- rbindlist(list(
  data.table(rank = "sk_base", scope = "broad",
             conc_index = ci_sk[["broad"]]),
  data.table(rank = "sk_base", scope = "direct",
             conc_index = ci_sk[["direct"]]),
  data.table(rank = names(ranking_vars), scope = "direct",
             conc_index = unname(ci_multi))
), use.names = TRUE, fill = TRUE)
fwrite(ci_dt, file.path(out_dir, "tab_concentration_index.csv"))

# -- 2. Balance: means by treat_type + normalized differences vs never --------

bal_vars <- c(
  log_pop_dens   = "Log pop density",
  log_steuerkraft= "log1p Steuerkraft",
  log_kaufkraft  = "Log Kaufkraft",
  bev_stock_p100k= "BEV stock p100k",
  ev_chargepoints_p100k = "Charge pts p100k",
  muni_gruene    = "Muni Grüne share",
  state_gruene   = "State Grüne share",
  fed_gruene     = "Fed Grüne share",
  n_vze_personal = "Personnel VZE p100k"
)

# Single pre-treatment snapshot per AGS8: earliest year in [2014,2016] with the
# variable observed (one row per (AGS8, variable)). Means by treat_type and
# normalized differences against the "never" cell.
bal_long <- rbindlist(lapply(names(bal_vars), function(v) {
  d <- dat[year %in% 2014:2016 & !is.na(get(v)),
           .(AGS8, treat_type, year, val = get(v))]
  setorder(d, AGS8, year)
  d <- d[, .SD[1L], by = AGS8]
  d[, var := v][]
}), fill = TRUE)

bal_means <- bal_long[, .(
  mean = mean(val, na.rm = TRUE),
  sd   = sd(val,   na.rm = TRUE),
  n    = .N
), by = .(var, treat_type)]

ref <- bal_means[treat_type == "never", .(var, mean_ref = mean, sd_ref = sd)]
bal_means <- merge(bal_means, ref, by = "var", all.x = TRUE)
bal_means[, norm_diff := (mean - mean_ref) / sqrt((sd^2 + sd_ref^2) / 2)]
bal_means[, var_label := vapply(var, function(v) bal_vars[[v]], character(1))]
setcolorder(bal_means, c("var", "var_label", "treat_type", "n", "mean", "sd", "norm_diff"))

fwrite(bal_means, file.path(out_dir, "tab_balance.csv"))

tex_lines <- c(
  "\\begin{tabular}{lllrrrr}",
  "\\hline",
  "Variable & Treat type & N & Mean & SD & Norm.\\ diff (vs never) \\\\",
  "\\hline",
  bal_means[, sprintf("%s & %s & %d & %.3f & %.3f & %.3f \\\\",
                      var_label, treat_type, n, mean, sd, norm_diff)],
  "\\hline",
  "\\end{tabular}"
)
writeLines(tex_lines, file.path(out_dir, "tab_balance.tex"))

# -- 2b. Direct vs Never: means + medians over 2015-latest --------------------
# Per-AGS8 within-unit mean over the 2015-latest window for each variable,
# then group-level mean (SE) with two-sided t-test and group-level median
# with Mood's median test (chi-squared on the 2x2 above/below pooled-median
# table). Sample restriction: treat_type ∈ {direct, never}; Stadtstaaten
# kept. N at the bottom of each table is the AGS8 count per group.

dat[, muni_gruene_pp  := muni_gruene  * 100]
dat[, state_gruene_pp := state_gruene * 100]
dat[, fed_gruene_pp   := fed_gruene   * 100]
dat[, xbev_100k       := xbev / 1e5]

desc_vars <- c(
  xbev_100k             = "Population (100k)",
  log_pop_dens          = "Log population density",
  log_steuerkraft       = "Log Steuerkraft p.c.",
  log_kaufkraft         = "Log Kaufkraft p.c.",
  muni_gruene_pp        = "Muni.\\ Grüne share (pp)",
  state_gruene_pp       = "State Grüne share (pp)",
  fed_gruene_pp         = "Fed.\\ Grüne share (pp)",
  staffing_ags8         = "Staffing (FTE, AGS5/N muni)",
  bev_stock_p100k       = "BEV stock p100k",
  ev_chargepoints_p100k = "Charge points p100k"
)

dat_dn <- dat[treat_type %in% c("direct", "never") & year >= 2015]
unit_means <- dat_dn[, lapply(.SD, mean, na.rm = TRUE),
                     by = .(AGS8, treat_type),
                     .SDcols = names(desc_vars)]

N_direct <- uniqueN(unit_means[treat_type == "direct", AGS8])
N_never  <- uniqueN(unit_means[treat_type == "never",  AGS8])

mood_p <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L) return(NA_real_)
  m <- median(c(x, y))
  ct <- matrix(c(sum(x > m), sum(x <= m),
                 sum(y > m), sum(y <= m)), nrow = 2L)
  if (any(rowSums(ct) == 0L) || any(colSums(ct) == 0L)) return(NA_real_)
  suppressWarnings(chisq.test(ct, correct = FALSE))$p.value
}

mean_rows <- rbindlist(lapply(names(desc_vars), function(v) {
  x <- unit_means[treat_type == "direct", get(v)]; x <- x[is.finite(x)]
  y <- unit_means[treat_type == "never" , get(v)]; y <- y[is.finite(y)]
  tt <- suppressWarnings(t.test(x, y, alternative = "two.sided"))
  data.table(
    var = v, label = desc_vars[[v]],
    direct_mean = mean(x),
    direct_se   = sd(x) / sqrt(length(x)),
    direct_n    = length(x),
    never_mean  = mean(y),
    never_se    = sd(y) / sqrt(length(y)),
    never_n     = length(y),
    p_value     = tt$p.value
  )
}))
fwrite(mean_rows, file.path(out_dir, "tab_desc_means.csv"))

med_rows <- rbindlist(lapply(names(desc_vars), function(v) {
  x <- unit_means[treat_type == "direct", get(v)]; x <- x[is.finite(x)]
  y <- unit_means[treat_type == "never" , get(v)]; y <- y[is.finite(y)]
  data.table(
    var = v, label = desc_vars[[v]],
    direct_median = median(x), direct_n = length(x),
    never_median  = median(y), never_n  = length(y),
    p_value       = mood_p(x, y)
  )
}))
fwrite(med_rows, file.path(out_dir, "tab_desc_medians.csv"))

# LaTeX twin: means
tex_means <- c(
  "\\begin{tabular}{lccc}",
  "\\hline",
  "Variable & Direct & Never & p-value (t-test) \\\\",
  "         & Mean (SE) & Mean (SE) &   \\\\",
  "\\hline",
  mean_rows[, sprintf("%s & %.3f (%.3f) & %.3f (%.3f) & %.3f \\\\",
                      label, direct_mean, direct_se,
                      never_mean, never_se, p_value)],
  "\\hline",
  sprintf("N (AGS8) & %d & %d & \\\\", N_direct, N_never),
  "\\hline",
  "\\end{tabular}"
)
writeLines(tex_means, file.path(out_dir, "tab_desc_means.tex"))

# LaTeX twin: medians
tex_meds <- c(
  "\\begin{tabular}{lccc}",
  "\\hline",
  "Variable & Direct & Never & p-value (Mood's median) \\\\",
  "\\hline",
  med_rows[, sprintf("%s & %.3f & %.3f & %.3f \\\\",
                     label, direct_median, never_median, p_value)],
  "\\hline",
  sprintf("N (AGS8) & %d & %d & \\\\", N_direct, N_never),
  "\\hline",
  "\\end{tabular}"
)
writeLines(tex_meds, file.path(out_dir, "tab_desc_medians.tex"))

cat(sprintf(
  "Direct vs Never descriptives: N_direct=%d, N_never=%d -> %s\n",
  N_direct, N_never, out_dir))

# -- 3. Treatment choropleth ---------------------------------------------------
# Optional: requires sf + the VG250 gpkg. Skips silently if either is missing.

gpkg_path <- file.path(root, "01_data", "04_shapefiles",
                       "vg250_01-01.utm32s.gpkg.ebenen",
                       "vg250_ebenen_0101", "DE_VG250.gpkg")
make_map <- function() {
  gem <- sf::st_read(gpkg_path, layer = "vg250_gem", quiet = TRUE)
  # Detect AGS column (could be "AGS", "ARS", or "ags")
  ags_col <- intersect(c("AGS", "AGS_0", "ARS", "ags"), names(gem))[1L]
  if (is.na(ags_col)) stop("No AGS column found in VG250 layer")
  if ("GF" %in% names(gem)) gem <- gem[gem$GF == 4, ]
  gem$AGS8 <- sprintf("%08d", as.integer(gem[[ags_col]]))
  # Use sf's tidy join via dplyr-style merge so geometry survives
  lkp <- as.data.frame(xs[, .(AGS8, treat_type)])
  gem <- merge(gem, lkp, by = "AGS8", all.x = TRUE)
  gem$treat_type[is.na(gem$treat_type)] <- "never"
  p_map <- ggplot(gem) +
    geom_sf(aes(fill = treat_type), color = NA) +
    scale_fill_manual(values = c(
      direct          = "#d73027",
      broadcast_only  = "#fdae61",
      never           = "#e0e0e0"
    )) +
    labs(fill = "Treatment type") +
    theme_void(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "fig_map_treatment.pdf"), p_map,
         width = 7, height = 9)
}
if (requireNamespace("sf", quietly = TRUE) && file.exists(gpkg_path)) {
  res <- tryCatch(make_map(), error = function(e) {
    cat("Skipped fig_map_treatment.pdf (",
        conditionMessage(e), ")\n", sep = "")
  })
} else {
  cat("Skipped fig_map_treatment.pdf (sf or VG250 gpkg not available)\n")
}

# -- 4. Correlation / VIF diagnostics on the hazard frame ----------------------
# z-scored on the hazard frame; we re-z here for consistency with 02_hazard.R.

ph[, log_dens_z     := z(log_pop_dens)]
ph[, sk_z           := z(log_steuerkraft_L1)]
ph[, kk_z           := z(log_kaufkraft_L1)]
ph[, bev_z          := z(log1p(pmax(bev_stock_p100k_L1, 0)))]
ph[, chg_z          := z(log1p(pmax(ev_chargepoints_p100k_L1, 0)))]
ph[, pers_z         := z(log1p(pmax(n_vze_personal_L1, 0)))]
ph[, muni_gruene_z  := z(muni_gruene_L1)]

X_cols <- c("log_dens_z", "sk_z", "kk_z", "bev_z", "chg_z", "pers_z", "muni_gruene_z")
X <- ph[, ..X_cols]
X <- X[complete.cases(X)]

corr_mat <- cor(X, use = "complete.obs")
fwrite(as.data.table(corr_mat, keep.rownames = "var"),
       file.path(out_dir, "tab_corr.csv"))

# VIF via auxiliary OLS regressions
vif_vals <- sapply(X_cols, function(v) {
  others <- setdiff(X_cols, v)
  fit <- lm(as.formula(sprintf("%s ~ %s", v, paste(others, collapse = " + "))),
            data = X)
  1 / (1 - summary(fit)$r.squared)
})
vif_dt <- data.table(var = X_cols, VIF = round(vif_vals, 3))
fwrite(vif_dt, file.path(out_dir, "tab_vif.csv"))

tex_lines <- c(
  "\\begin{tabular}{lr}",
  "\\hline",
  "Variable & VIF \\\\",
  "\\hline",
  vif_dt[, sprintf("%s & %.2f \\\\", var, VIF)],
  "\\hline",
  "\\end{tabular}"
)
writeLines(tex_lines, file.path(out_dir, "tab_vif.tex"))

# -- 5. Dual BEV choropleth, 2023 -----------------------------------------------
# Left panel:  BEV stock per 100k population, winsorised at 99th pct.
# Right panel: BEV + hybrid share of total vehicle stock, winsorised at 99th pct.
# Same sequential blue palette (UBonn blue anchor). Sized for ~1/3 of an A4
# page. Rendered via ragg with SF Pro Display Medium where available.

make_bev_dual_map <- function() {
  gem <- sf::st_read(gpkg_path, layer = "vg250_gem", quiet = TRUE)
  ags_col <- intersect(c("AGS", "AGS_0", "ARS", "ags"), names(gem))[1L]
  if (is.na(ags_col)) stop("No AGS column found in VG250 layer")
  if ("GF" %in% names(gem)) gem <- gem[gem$GF == 4, ]
  gem$AGS8 <- sprintf("%08d", as.integer(gem[[ags_col]]))

  # KBA stock data
  kba <- fread(file.path(root, "01_data", "02_intermediate", "kba",
                          "kba_ags8_panel.csv"),
               select = c("AGS8", "year", "B_elektro_overall",
                          "B_hybrid_overall", "B_total_overall"),
               colClasses = list(character = "AGS8"))
  kba23 <- kba[year == 2023L]

  # Population from the main panel (needed for p100k)
  pan23 <- fread(file.path(root, "01_data", "03_final", "emk_inkar_panel_ags8.csv"),
                 select = c("AGS8", "year", "xbev"),
                 colClasses = list(character = "AGS8"))[year == 2023L]

  dat23 <- merge(kba23, pan23[, .(AGS8, xbev)], by = "AGS8", all.x = TRUE)

  # Left: BEV stock p100k, winsorised at 99th pct
  dat23[, bev_p100k := B_elektro_overall / xbev * 100000]
  cap_p100k <- quantile(dat23$bev_p100k, 0.99, na.rm = TRUE)
  dat23[, bev_p100k_w := pmin(bev_p100k, cap_p100k)]
  # Count-weighted national average
  natl_avg_p100k <- dat23[, sum(B_elektro_overall, na.rm = TRUE) /
                              sum(xbev, na.rm = TRUE) * 100000]

  # Right: BEV share (pure electric only), winsorised at 99th pct
  dat23[, alt_share := B_elektro_overall / B_total_overall]
  cap_share <- quantile(dat23$alt_share, 0.99, na.rm = TRUE)
  dat23[, alt_share_w := pmin(alt_share, cap_share)]
  natl_avg_share <- dat23[, sum(B_elektro_overall, na.rm = TRUE) /
                              sum(B_total_overall, na.rm = TRUE)]

  gem <- merge(gem, as.data.frame(dat23[, .(AGS8, bev_p100k_w, alt_share_w)]),
               by = "AGS8", all.x = TRUE)

  ubonn_blue <- "#004CFF"
  pal <- c("#f7fbff", "#deebf7", "#9ecae1", "#4292c6", "#08519c", ubonn_blue)

  font_family <- "sans"
  sf_pro_path <- "/Library/Fonts/SF-Pro-Display-Medium.otf"
  if (requireNamespace("systemfonts", quietly = TRUE) && file.exists(sf_pro_path)) {
    systemfonts::register_font(name = "SFProDisplayMedium", plain = sf_pro_path)
    font_family <- "SFProDisplayMedium"
  }

  # ~1/3 of A4 page (text-width × 1/3 page height)
  fig_width_in  <- 6.5
  fig_height_in <- 4.0
  bar_height_in <- fig_height_in * 0.52

  # Annotation centered under each map in data coordinates; clip = "off" lets
  # it render below the bounding box. Bottom plot.margin reserves the space.
  bbox   <- sf::st_bbox(gem)
  capt_x <- (bbox["xmin"] + bbox["xmax"]) / 2
  capt_y <- bbox["ymin"] - 0.055 * (bbox["ymax"] - bbox["ymin"])

  base_theme <- theme_void(base_size = 8, base_family = font_family) +
    theme(
      legend.position    = "right",
      legend.box.spacing = grid::unit(0.15, "cm"),
      legend.title       = element_text(size = 6.5, lineheight = 1.15),
      legend.text        = element_text(size = 6),
      plot.margin        = margin(t = 1, r = 2, b = 16, l = 2),
      plot.title         = element_text(size = 7.5, hjust = 0.5,
                                        margin = margin(b = 2))
    )

  p_left <- ggplot(gem) +
    labs(title = "BEV stock per 100k population") +
    geom_sf(aes(fill = bev_p100k_w), color = NA) +
    scale_fill_gradientn(
      colours  = pal,
      na.value = "grey85",
      name     = "BEV stock\nper 100k\n(2023)",
      guide    = guide_colorbar(
        barheight = grid::unit(bar_height_in, "in"),
        barwidth  = grid::unit(0.18, "in")
      )
    ) +
    annotate("text", x = capt_x, y = capt_y,
             label = sprintf("Germany avg.: %.0f per 100k", natl_avg_p100k),
             family = font_family, size = 5.5 / .pt, colour = "black",
             hjust = 0.5, vjust = 1) +
    coord_sf(clip = "off") +
    base_theme

  p_right <- ggplot(gem) +
    labs(title = "BEV share of vehicle stock") +
    geom_sf(aes(fill = alt_share_w), color = NA) +
    scale_fill_gradientn(
      colours  = pal,
      labels   = scales::percent,
      na.value = "grey85",
      name     = "BEV share of\nvehicle stock\n(2023)",
      guide    = guide_colorbar(
        barheight = grid::unit(bar_height_in, "in"),
        barwidth  = grid::unit(0.18, "in")
      )
    ) +
    annotate("text", x = capt_x, y = capt_y,
             label = sprintf("Germany avg.: %s",
                             scales::percent(natl_avg_share, accuracy = 0.1)),
             family = font_family, size = 5.5 / .pt, colour = "black",
             hjust = 0.5, vjust = 1) +
    coord_sf(clip = "off") +
    base_theme

  p_combined <- p_left + p_right +
    patchwork::plot_layout(ncol = 2) &
    theme(plot.background = element_rect(fill = "white", colour = NA))

  ggsave(file.path(out_dir, "fig_map_bev_dual_2023.png"), p_combined,
         width = fig_width_in, height = fig_height_in, dpi = 300,
         device = ragg::agg_png)
  cat("Saved fig_map_bev_dual_2023.png\n")
}
if (requireNamespace("sf", quietly = TRUE) &&
    requireNamespace("scales", quietly = TRUE) &&
    requireNamespace("ragg", quietly = TRUE) &&
    requireNamespace("patchwork", quietly = TRUE) &&
    file.exists(gpkg_path)) {
  tryCatch(make_bev_dual_map(), error = function(e) {
    cat("Skipped fig_map_bev_dual_2023.png (",
        conditionMessage(e), ")\n", sep = "")
  })
} else {
  cat("Skipped fig_map_bev_dual_2023.png (sf/scales/ragg/patchwork or VG250 gpkg not available)\n")
}

cat(sprintf("\nDescriptives written -> %s\n", out_dir))
