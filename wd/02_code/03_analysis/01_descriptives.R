# ───────────────────────────────────────────────────────────────────────────────
# 01_descriptives.R
#
# Inequality storyline: who is COVERED, not how many euros are received.
# Grant sizes barely vary; the distributive question is binary access.
#
# Outputs (04_results/01_descriptives/):
#   fig_concentration_coverage.pdf          headline Lorenz-style curve (kk_base)
#   fig_concentration_coverage_sk.pdf       same, ranked by sk_base (appendix)
#   tab_receipt_by_quintile.{tex,csv}       receipt prob × quintile × {broad,direct}
#   tab_balance.{tex,csv}                   means by treat_type + normalized diffs
#   fig_map_treatment.pdf                   VG250 choropleth of treat_type
#   tab_corr_vif.{tex,csv}                  pairwise corr + VIF on hazard frame
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

dat <- readRDS(file.path(data_final, "frame_did_broad.rds"))
ph  <- readRDS(file.path(data_final, "frame_hazard.rds"))

# Cross-section: one row per AGS8 (use the unit's earliest base-window year)
xs <- dat[, .SD[1L], by = AGS8]

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

ci_kk <- c(
  broad  = conc_index(xs, "kk_base", "pop_base", "t_broad"),
  direct = conc_index(xs, "kk_base", "pop_base", "t_direct")
)
ci_sk <- c(
  broad  = conc_index(xs, "sk_base", "pop_base", "t_broad"),
  direct = conc_index(xs, "sk_base", "pop_base", "t_direct")
)

cat("Concentration indices (kk_base):  broad ",
    round(ci_kk["broad"], 3), " | direct ", round(ci_kk["direct"], 3), "\n")
cat("Concentration indices (sk_base):  broad ",
    round(ci_sk["broad"], 3), " | direct ", round(ci_sk["direct"], 3), "\n")

lz_kk <- build_lorenz(xs, "kk_base")
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
      x = paste0("Cumulative population share (Gemeinden ranked by ", rank_label, ")"),
      y = "Cumulative share of treated population",
      color = NULL, linetype = NULL,
      caption = sprintf("Concentration index — broad: %.3f, direct: %.3f",
                        ci["broad"], ci["direct"])
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

ggsave(file.path(out_dir, "fig_concentration_coverage.pdf"),
       plot_lorenz(lz_kk, "baseline log Kaufkraft", ci_kk),
       width = 7, height = 5)
ggsave(file.path(out_dir, "fig_concentration_coverage_sk.pdf"),
       plot_lorenz(lz_sk, "baseline log Steuerkraft", ci_sk),
       width = 7, height = 5)

fwrite(data.table(
  rank = c(rep("kk_base", 2L), rep("sk_base", 2L)),
  scope = rep(c("broad", "direct"), 2L),
  conc_index = c(ci_kk, ci_sk)
), file.path(out_dir, "tab_concentration_index.csv"))

# -- 2. Receipt probability by quintile ---------------------------------------

receipt_tab <- function(xs, q_col, treat_col) {
  d <- xs[!is.na(get(q_col))]
  d[, .(
    n_units    = .N,
    pop_share  = sum(pop_base, na.rm = TRUE),
    receipt    = mean(get(treat_col)),
    covered_pop = sum(pop_base * get(treat_col), na.rm = TRUE)
  ), by = q_col][order(get(q_col))]
}

receipt_rows <- rbindlist(list(
  receipt_tab(xs, "kk_base_q5", "t_broad" )[, `:=`(rank = "kk_base", scope = "broad" )],
  receipt_tab(xs, "kk_base_q5", "t_direct")[, `:=`(rank = "kk_base", scope = "direct")],
  receipt_tab(xs, "sk_base_q5", "t_broad" )[, `:=`(rank = "sk_base", scope = "broad" )],
  receipt_tab(xs, "sk_base_q5", "t_direct")[, `:=`(rank = "sk_base", scope = "direct")]
), use.names = TRUE, fill = TRUE)

# Total covered-population share within scope
receipt_rows[, pop_cov_share := covered_pop / sum(covered_pop), by = .(rank, scope)]

fwrite(receipt_rows, file.path(out_dir, "tab_receipt_by_quintile.csv"))

# Compact LaTeX twin: receipt prob × scope, rank
tex_lines <- c(
  "\\begin{tabular}{llrrrr}",
  "\\hline",
  "Rank & Q & Scope & N units & P(receipt) & Pop cov.\\ share \\\\",
  "\\hline",
  receipt_rows[, sprintf("%s & %s & %s & %d & %.3f & %.3f \\\\",
                         rank,
                         ifelse(is.na(kk_base_q5), as.character(sk_base_q5),
                                as.character(kk_base_q5)),
                         scope, n_units, receipt, pop_cov_share)],
  "\\hline",
  "\\end{tabular}"
)
writeLines(tex_lines, file.path(out_dir, "tab_receipt_by_quintile.tex"))

# -- 3. Balance: means by treat_type + normalized differences vs never --------

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

# -- 4. Treatment choropleth ---------------------------------------------------
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

# -- 5. Correlation / VIF diagnostics on the hazard frame ----------------------
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

cat(sprintf("\nDescriptives written -> %s\n", out_dir))
