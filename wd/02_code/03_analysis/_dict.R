# Shared variable / label conventions for analysis scripts.
# Sourced from each NN_*.R via `source(file.path(code_dir, "03_analysis", "_dict.R"))`.

# Variable labels for etable / modelsummary
dict <- c(
  # IDs / FEs
  AGS8 = "AGS8 (Municipality)",
  AGS5 = "AGS5 (County)",
  AGS2 = "AGS2 (State)",
  year = "Year",

  # Hazard channels (lagged, z-scored on estimation sample)
  log_dens_z      = "Log pop. density (z)",
  log_dens        = "Log pop. density",
  sk_z            = "Log tax capacity (z, L1)",
  bev_z           = "$\\sqrt{\\cdot}$ BEV stock p100k (z, L1)",
  chg_z           = "$\\sqrt{\\cdot}$ EV chargers p100k (z, L1)",
  pers_z          = "Log Personnel FTE p100k (z, L1)",
  muni_gruene_z   = "Muni Green share (z, L1)",
  state_gruene_z  = "State Green share (z, L1)",
  fed_gruene_z    = "Fed Green share (z, L1)",
  kreis_funded    = "County-funded (strict past)",

  # Baseline snapshot z-scores (DiD covariates)
  sk_base_z          = "Baseline log tax capacity (z)",
  dens_base_z        = "Baseline log pop. density (z)",
  green_base_z       = "Baseline muni Green share (z)",
  state_green_base_z = "Baseline state Green share (z)",
  bev_base_z         = "Baseline $\\sqrt{\\cdot}$ BEV stock p100k (z)",
  chg_base_z         = "Baseline $\\sqrt{\\cdot}$ EV chargers p100k (z)",

  # Heterogeneity grouping
  east               = "East Germany (1/0)",

  # Events
  onset_direct = "Direct-treatment onset",
  onset_broad  = "Coverage-event onset",

  # Outcomes
  bev_neuzulassungen_p100k = "BEV new registrations (p100k)",
  bev_corporate_p100k      = "BEV new registrations, corporate (p100k)",
  bev_private_p100k        = "BEV new registrations, private (p100k)",
  bev_stock_p100k          = "BEV stock (p100k)",
  ice_neuzulassungen_p100k = "ICE new registrations (p100k) — placebo",
  bev_share_pct            = "BEV market share (pp, win.)"
)

# Outcome -> human label
OUTCOME_LABELS <- c(
  bev_neuzulassungen_p100k = "BEV new registrations per 100k population",
  bev_corporate_p100k      = "BEV new registrations (corporate) per 100k",
  bev_private_p100k        = "BEV new registrations (private) per 100k",
  bev_stock_p100k          = "BEV stock per 100k population",
  ice_neuzulassungen_p100k = "ICE new registrations per 100k (placebo)",
  bev_share_pct            = "BEV market share (percentage points, win. 99th pct)"
)

# Event-study horizon for plots / aggte. Each downstream script computes a
# data-driven cap via `es_max_data_driven()` below and passes it to aggte().
ES_MIN <- -4L
ES_MAX <-  7L

# Stadtstaaten (n_vze_personal conflates municipal / Länder roles)
STADTSTAATEN <- c("02", "04", "11")

# East Germany split (neue Länder). Berlin (11) is treated as East here; to run
# the contrast excluding city-states, drop AGS2 %in% STADTSTAATEN first. Used by
# the East/West heterogeneity cut in 06_heterogeneity.R.
EAST_AGS2 <- c("11", "12", "13", "14", "15", "16")

# FRL / EmoG funding-regime break. Treated cohorts <= FRL_CUTOFF (2020) vs
# >= 2021. Used by the pre/post-2021 heterogeneity cut in 06_heterogeneity.R.
FRL_CUTOFF <- 2020L

# Minimum treated-unit count per cohort. Cohorts below this threshold are
# dropped entirely from DiD frames in 00_prep_analysis.R.
COHORT_MIN <- 5L

# Minimum number of cohorts required at event time e for es_max_data_driven()
# to include e in the display horizon.
MIN_COHORTS_PER_E <- 3L

# Canonical CS-dr xformla — used by 04, 05, 06, 07. Baseline (z'd) covariates
# computed as per-AGS8 means over CS_BASE_WINDOW = 2013:2015 (set in
# 00_prep_analysis.R). bev_base_z and chg_base_z excluded: baseline 2013-15
# BEV/charge counts are mass-zero across AGS8 and collapse DR design matrices.
XFORMLA_CS <- ~ sk_base_z + state_green_base_z + dens_base_z

# Data-driven upper horizon for aggte(type = "dynamic"). Returns the largest e
# (>= 1) such that at least MIN_COHORTS_PER_E eligible cohorts have
# post_avail >= e. Reads cohorts from gname_col. Falls back to 1L if fewer
# than MIN_COHORTS_PER_E cohorts exist at e = 1.
es_max_data_driven <- function(dat, yname, gname_col = "gname_cs") {
  yr_max <- dat[!is.na(get(yname)), max(year)]
  cohorts <- sort(unique(dat[get(gname_col) > 0L, get(gname_col)]))
  if (!is.finite(yr_max) || length(cohorts) == 0L) return(1L)
  post_avail <- yr_max - cohorts                      # vector, one per cohort
  # for each candidate e, count cohorts with post_avail >= e
  max_possible <- max(post_avail)
  if (max_possible < 1L) return(1L)
  counts <- vapply(seq_len(max_possible),
                   function(e) sum(post_avail >= e), integer(1))
  ok <- which(counts >= MIN_COHORTS_PER_E)
  if (length(ok) == 0L) return(1L)
  max(1L, max(ok))
}

# Standard z helper, used after every sample restriction so each model's
# regressors are standardised on its own estimation sample.
z <- function(x) as.numeric(scale(x))

# Population-weighted quantile cuts. Returns the q-th sample quantiles of x
# weighted by w (NA-safe).
wq <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x  <- x[ok]; w <- w[ok]
  o  <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1L]], numeric(1))
}

# CSV+TeX twin writer for fixest etable. Writes a .tex and a .csv with the
# same coefficient/SE numbers, plus a manifest row.
write_estimates_csv <- function(models, file) {
  rows <- lapply(seq_along(models), function(i) {
    m  <- models[[i]]
    nm <- names(models)[i]
    co <- coef(m); se <- sqrt(diag(vcov(m)))
    data.table(
      model = nm,
      term  = names(co),
      estimate = as.numeric(co),
      se       = as.numeric(se[names(co)]),
      n        = nobs(m)
    )
  })
  dt <- data.table::rbindlist(rows)
  dt[, ci_lo := estimate - 1.96 * se]
  dt[, ci_hi := estimate + 1.96 * se]
  data.table::fwrite(dt, file)
  invisible(dt)
}

# Resolve project root from the running script path. Drop-in for the
# argv/self_flag block at the top of every analysis script.
resolve_root <- function() {
  argv      <- commandArgs(trailingOnly = FALSE)
  self_flag <- grep("--file=", argv, value = TRUE)
  self <- if (length(self_flag)) {
    normalizePath(sub("--file=", "", self_flag))
  } else if (
    requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()
  ) {
    normalizePath(rstudioapi::getSourceEditorContext()$path)
  } else {
    stop("Cannot determine script path. Run as: Rscript <script>.R")
  }
  dirname(dirname(dirname(self)))
}

# Standard output directory for a script: 03_output/<script_stem>/
results_dir <- function(root, stem) {
  d <- file.path(root, "03_output", stem)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

# ── Table writers ──────────────────────────────────────────────────────────────

# Core table writer (tabularray package).
#
# env = "tblr" (default): single-page table wrapped in a \begin{table} float
#   with \centering, \caption, \label, and a threeparttable+tablenotes block
#   for notes below the tblr body. Use for all coefficient and event-study
#   tables that fit on one page.
#
# env = "longtblr": self-contained long-table float; caption, label, and note{}
#   live in the outer [...] spec. Use only for HR/AME joint tables that span
#   multiple columns with \SetCell spanning (these cannot live inside a box).
write_longtblr <- function(stem, caption, label, note, colspec,
                            header_rows, body_rows, footer_rows,
                            rowsep = "0pt", env = "tblr",
                            resize = FALSE) {
  tblr_body <- c(
    paste0("\\begin{tblr}{"),
    paste0("  colspec = {", colspec, "},"),
    paste0("  rowsep  = ", rowsep, ","),
    "}",
    "\\hline\\hline",
    header_rows,
    "\\hline",
    body_rows,
    if (length(footer_rows)) c("\\hline", footer_rows) else NULL,
    "\\hline\\hline",
    "\\end{tblr}"
  )

  if (env == "tblr") {
    if (resize) {
      tblr_inner <- tblr_body[-length(tblr_body)]  # drop \end{tblr}
      tbl_body_tex <- c(
        "\\resizebox{\\textwidth}{!}{%",
        tblr_inner,
        "\\end{tblr}%",
        "}"
      )
    } else {
      tbl_body_tex <- tblr_body
    }
    note_block <- if (nzchar(paste(note, collapse = ""))) {
      c("\\begin{tablenotes}", note, "\\end{tablenotes}")
    } else character(0)
    tex <- c(
      "\\begin{table}[!ht]",
      "\\centering",
      paste0("\\caption{", caption, "}"),
      paste0("\\label{", label, "}"),
      tbl_body_tex,
      note_block,
      "\\end{table}"
    )
  } else {
    tex <- c(
      sprintf("\\begin{%s}[", env),
      paste0("    caption = {", caption, "},"),
      paste0("    label = {", label, "},"),
      paste0("    note{} = {\\small ", note, "},"),
      "]{",
      paste0("  colspec = {", colspec, "},"),
      paste0("  rowsep  = ", rowsep, ","),
      "}",
      "\\hline\\hline",
      header_rows,
      "\\hline",
      body_rows,
      if (length(footer_rows)) c("\\hline", footer_rows) else NULL,
      "\\hline\\hline",
      sprintf("\\end{%s}", env)
    )
  }
  writeLines(tex, paste0(stem, ".tex"))
  invisible(tex)
}

# Companion .tex for a figure: a LaTeX figure float (caption + label + a
# scriptsize note) pointing at a rendered image, mirroring the table convention
# so captions/notes live in TeX rather than baked into the plot. `img_file` is
# the full path to the .png/.pdf; the .tex is written to the same stem;
# \includegraphics resolves <graphics_prefix>/<basename>. `note` defaults to
# empty — callers pass the methodological note explicitly.
write_fig_tex <- function(img_file, caption, label, note = "",
                          graphics_prefix = "0_3_figures",
                          width = "\\textwidth") {
  base <- basename(img_file)
  stem <- sub("\\.(png|pdf)$", "", img_file)
  tex <- c(
    "\\begin{figure}[!ht]",
    "    \\centering",
    sprintf("    \\includegraphics[width=%s]{%s/%s}",
            width, graphics_prefix, base),
    sprintf("    \\caption{%s}", caption),
    sprintf("    \\label{%s}", label),
    if (nzchar(note)) c(
      "    \\begin{figurenotes}",
      sprintf("        %s", note),
      "    \\end{figurenotes}"
    ) else NULL,
    "\\end{figure}"
  )
  writeLines(tex, paste0(stem, ".tex"))
  invisible(tex)
}

# Coefficient table writer for fixest models. Builds body from a `groups` list
# (names = section labels used as blank-line separators, values = char vector
# of term names; NULL value inserts a \\hline separator).
# event_counts: named list model_name -> integer, for "Onset events" footer row.
write_coef_longtblr <- function(
  models, stem, caption, label, note, groups, var_labels,
  event_counts = NULL, show_pr2 = TRUE, env = "tblr", resize = FALSE
) {
  col_names <- names(models)
  n_spec    <- length(col_names)

  .stars <- function(p) ifelse(is.na(p), "",
    ifelse(p < 0.01, "$^{***}$", ifelse(p < 0.05, "$^{**}$", ifelse(p < 0.1, "$^{*}$", ""))))

  .cell <- function(nm, term) {
    ct <- tryCatch(coeftable(models[[nm]]), error = function(e) NULL)
    if (is.null(ct) || !term %in% rownames(ct)) return(c(est = "", se = ""))
    b <- ct[term, "Estimate"]; s <- ct[term, "Std. Error"]; p <- ct[term, ncol(ct)]
    c(est = sprintf("%.3f%s", b, .stars(p)), se = sprintf("(%.3f)", s))
  }

  .coef_rows <- function(term) {
    lab   <- var_labels[term]; if (is.na(lab)) lab <- term
    cells <- lapply(col_names, .cell, term = term)
    r1    <- c(paste0("\\quad ", lab), vapply(cells, `[[`, character(1), "est"))
    r2    <- c("",                     vapply(cells, `[[`, character(1), "se"))
    c(paste(r1, collapse = " & "), " \\\\",
      paste(r2, collapse = " & "), " \\\\")
  }

  body_rows <- c()
  for (i in seq_along(groups)) {
    terms <- groups[[i]]
    if (is.null(terms)) {
      body_rows <- c(body_rows, "\\hline")
    } else {
      if (i > 1L && !is.null(groups[[i - 1L]])) body_rows <- c(body_rows, "")
      body_rows <- c(body_rows, unlist(lapply(terms, .coef_rows)))
    }
  }

  header_rows <- paste0(paste(c("", col_names), collapse = " & "), " \\\\")

  fmt_n <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)
  n_vals <- setNames(
    vapply(col_names, function(nm)
      tryCatch(nobs(models[[nm]]), error = function(e) NA_integer_), integer(1)),
    col_names)
  pr2_vals <- if (show_pr2)
    setNames(vapply(col_names, function(nm)
      tryCatch(as.numeric(r2(models[[nm]], "pr2")), error = function(e) NA_real_),
      numeric(1)), col_names)
  else NULL

  footer_rows <- c(
    paste0(paste(c("Year/State FE", rep("Yes", n_spec)), collapse = " & "), " \\\\"),
    paste0(paste(c("Observations",  fmt_n(n_vals)),       collapse = " & "), " \\\\")
  )
  if (!is.null(event_counts)) {
    ev_vals <- vapply(col_names, function(nm) fmt_n(event_counts[[nm]]), character(1))
    footer_rows <- c(footer_rows,
      paste0(paste(c("Onset events", ev_vals), collapse = " & "), " \\\\"))
  }
  if (!is.null(pr2_vals)) {
    footer_rows <- c(footer_rows,
      paste0(paste(c("Pseudo-$R^2$", sprintf("%.3f", pr2_vals)), collapse = " & "), " \\\\"))
  }

  write_longtblr(stem = stem, caption = caption, label = label, note = note,
                 colspec = sprintf("l *{%d}{l}", n_spec),
                 header_rows = header_rows,
                 body_rows   = body_rows,
                 footer_rows = footer_rows,
                 env         = env,
                 resize      = resize)
  write_estimates_csv(models, paste0(stem, ".csv"))
  invisible(models)
}
