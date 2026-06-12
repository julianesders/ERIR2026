# ───────────────────────────────────────────────────────────────────────────────
# 07_assemble.R
#
# Walk 04_results/ and build a manifest of every CSV / TeX / PDF written by
# the analysis pipeline. One row per file: file path, script that produced it
# (inferred from the parent directory), content type, last-modified timestamp.
#
# Output: 04_results/manifest.csv
# ───────────────────────────────────────────────────────────────────────────────

library(data.table)

argv      <- commandArgs(trailingOnly = FALSE)
self_flag <- grep("--file=", argv, value = TRUE)
self <- if (length(self_flag)) {
  normalizePath(sub("--file=", "", self_flag))
} else if (
  requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()
) {
  normalizePath(rstudioapi::getSourceEditorContext()$path)
} else {
  stop("Cannot determine script path. Run as: Rscript 07_assemble.R")
}
root        <- dirname(dirname(dirname(self)))
results_dir <- file.path(root, "04_results")

files <- list.files(results_dir, recursive = TRUE, full.names = TRUE,
                    pattern = "\\.(csv|tex|pdf|txt)$")
files <- files[!grepl("manifest\\.csv$", files)]

manifest <- data.table(
  file       = sub(paste0("^", results_dir, "/"), "", files),
  script_dir = basename(dirname(files)),
  ext        = tools::file_ext(files),
  mtime      = file.info(files)$mtime,
  size_kb    = round(file.info(files)$size / 1024, 1)
)

manifest[, content := fcase(
  ext == "csv" & grepl("^est_",  basename(file)), "estimates",
  ext == "csv" & grepl("^tab_",  basename(file)), "table",
  ext == "tex",                                   "table_tex",
  ext == "pdf" & grepl("^fig_",  basename(file)), "figure",
  ext == "pdf" & grepl("^es_",   basename(file)), "event_study",
  ext == "txt",                                    "log",
  default = "other"
)]

setorder(manifest, script_dir, file)
fwrite(manifest, file.path(results_dir, "manifest.csv"))

cat(sprintf("Manifest: %d files indexed.\n", nrow(manifest)))
cat(sprintf("By content:\n"))
print(manifest[, .N, by = content])
cat(sprintf("\nWritten -> %s/manifest.csv\n", results_dir))
