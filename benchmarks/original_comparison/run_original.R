#!/usr/bin/env Rscript

parse_options <- function(arguments) {
  if (length(arguments) %% 2L != 0L) {
    stop("Every option must be followed by one value", call. = FALSE)
  }
  result <- list()
  if (!length(arguments)) return(result)
  for (index in seq.int(1L, length(arguments), by = 2L)) {
    name <- sub("^--", "", arguments[[index]])
    result[[name]] <- arguments[[index + 1L]]
  }
  result
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
required <- c(
  "rda", "bulk", "decepticonx-r", "decepticon-methods-r", "helper-dir",
  "output"
)
missing <- required[!required %in% names(options)]
if (length(missing)) {
  stop("Missing option(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

normalise_file <- function(path) normalizePath(path, mustWork = TRUE)
rda_path <- normalise_file(options[["rda"]])
bulk_path <- normalise_file(options[["bulk"]])
decepticonx_r <- normalise_file(options[["decepticonx-r"]])
decepticon_methods_r <- normalise_file(options[["decepticon-methods-r"]])
helper_dir <- normalizePath(options[["helper-dir"]], mustWork = TRUE)
output_dir <- normalizePath(
  options[["output"]], mustWork = FALSE
)
seed <- if ("seed" %in% names(options)) as.integer(options[["seed"]]) else 20260814L
if (is.na(seed)) stop("--seed must be an integer", call. = FALSE)

if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory must be absent or empty: ", output_dir, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "custom_signature_matrix"), showWarnings = FALSE)
dir.create(file.path(output_dir, "res"), showWarnings = FALSE)

timing_path <- file.path(output_dir, "timings.tsv")
warning_path <- file.path(output_dir, "warnings.tsv")
writeLines("stage\tstatus\telapsed_seconds\tuser_seconds\tsystem_seconds\terror", timing_path)
writeLines("stage\twarning", warning_path)

clean_field <- function(value) {
  value <- gsub("[\t\r\n]+", " ", as.character(value))
  trimws(value)
}

append_warning <- function(stage, message) {
  cat(
    clean_field(stage), clean_field(message), sep = "\t", file = warning_path,
    append = TRUE
  )
  cat("\n", file = warning_path, append = TRUE)
}

run_stage <- function(name, expression) {
  warnings_seen <- character()
  error_text <- ""
  status <- "ok"
  value <- NULL
  gc(verbose = FALSE)
  timing <- system.time({
    value <- tryCatch(
      withCallingHandlers(
        eval.parent(substitute(expression)),
        warning = function(condition) {
          warnings_seen <<- c(warnings_seen, conditionMessage(condition))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(condition) {
        status <<- "error"
        error_text <<- conditionMessage(condition)
        NULL
      }
    )
  })
  if (length(warnings_seen)) {
    for (message in unique(warnings_seen)) append_warning(name, message)
  }
  fields <- c(
    name,
    status,
    unname(timing[["elapsed"]]),
    unname(timing[["user.self"]]),
    unname(timing[["sys.self"]]),
    clean_field(error_text)
  )
  cat(paste(fields, collapse = "\t"), "\n", sep = "", file = timing_path, append = TRUE)
  if (status != "ok") stop("Stage ", name, " failed: ", error_text, call. = FALSE)
  invisible(value)
}

sha256 <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  if (!length(output)) return(NA_character_)
  strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
}

canonical_mapping_sha256 <- function(mapping) {
  if (!is.numeric(mapping) || is.null(names(mapping)) || any(!nzchar(names(mapping)))) {
    stop("EPIC mRNA mapping must be a named numeric vector", call. = FALSE)
  }
  if (anyDuplicated(names(mapping)) || anyNA(mapping) || any(!is.finite(mapping))) {
    stop("EPIC mRNA mapping contains duplicate, missing, or non-finite values", call. = FALSE)
  }
  keys <- sort(names(mapping), method = "radix")
  canonical <- paste(
    sprintf("%s\t%.17g", keys, as.numeric(mapping[keys])),
    collapse = "\n"
  )
  path <- tempfile("epic-mapping-")
  on.exit(unlink(path), add = TRUE)
  writeChar(canonical, path, eos = NULL, useBytes = TRUE)
  sha256(path)
}

git_head <- function(path) {
  directory <- normalizePath(dirname(path), mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(directory, ".git")) || file.exists(file.path(directory, ".git"))) break
    parent <- dirname(directory)
    if (identical(parent, directory)) return(NA_character_)
    directory <- parent
  }
  result <- system2("git", c("-C", shQuote(directory), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)
  if (length(result)) result[[1]] else NA_character_
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(seed)
started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
setwd(output_dir)

run_stage("input_load", {
  suppressPackageStartupMessages(library(Biobase))
  loaded_names <- load(rda_path, envir = environment())
  if (length(loaded_names) != 1L) stop("RDA must contain exactly one object")
  expression_set <- get(loaded_names[[1]], envir = environment())
  if (!methods::is(expression_set, "ExpressionSet")) stop("RDA object is not an ExpressionSet")
  single_cell <- Biobase::exprs(expression_set)
  metadata <- Biobase::pData(expression_set)
  if (!"cell_type" %in% colnames(metadata)) stop("ExpressionSet has no cell_type metadata")
  subtype <- as.matrix(metadata[, "cell_type", drop = FALSE])
  rownames(subtype) <- colnames(single_cell)
  keep <- !is.na(subtype[, 1L]) & subtype[, 1L] != "Unknown"
  single_cell <- single_cell[, keep, drop = FALSE]
  subtype <- subtype[keep, , drop = FALSE]
  bulk <- utils::read.table(
    bulk_path, row.names = 1, header = TRUE, sep = "\t", check.names = FALSE
  )
  if (anyNA(bulk) || any(!is.finite(as.matrix(bulk)))) stop("Bulk contains missing/non-finite values")
})

source(decepticonx_r, local = globalenv())
source(decepticon_methods_r, local = globalenv())

epic_mrna_cell <- getExportedValue("EPIC", "mRNA_cell_default")
epic_mrna_cell_sha256 <- canonical_mapping_sha256(epic_mrna_cell)
epic_mrna_cell_entries <- length(epic_mrna_cell)

run_stage("ref_bayesprism", BayesPrism_base(bulk, single_cell, subtype))
run_stage("ref_monocle3", Monocle3_base(single_cell, subtype))
run_stage("ref_music2_legacy", MuSiC2_base(single_cell, bulk, subtype))

signature_paths <- file.path(
  "custom_signature_matrix",
  c("Monocle3_base.txt", "BayesPrism_base.txt", "MuSiC2_base.txt")
)
if (!all(file.exists(signature_paths))) stop("One or more original reference files are missing")

helpers <- c("CIBERSORT.R", "music_prop.R", "music_basis.R")
helper_sources <- file.path(helper_dir, helpers)
if (!all(file.exists(helper_sources))) {
  stop("Helper directory is missing: ", paste(helpers[!file.exists(helper_sources)], collapse = ", "))
}
if (!all(file.copy(helper_sources, file.path(output_dir, helpers), overwrite = FALSE))) {
  stop("Could not stage one or more runtime helper files")
}
runtime_helpers <- file.path(output_dir, helpers)
tryCatch(
  {
    set.seed(seed)
    run_stage("deconv_cibersort", DECEPTICON_ciber_custom(bulk_path, signature_paths))
    set.seed(seed)
    run_stage("deconv_cibersort_abs", DECEPTICON_ciber_abs_custom(bulk_path, signature_paths))
    set.seed(seed)
    run_stage("deconv_epic", DECEPTICON_epic_custom(bulk_path, signature_paths))
    set.seed(seed)
    run_stage("deconv_deconrnaseq", DECEPTICON_decon_custom(bulk_path, signature_paths))

    bulk_eset <- run_stage("music_bulk_eset", {
      suppressPackageStartupMessages(library(SCDC))
      fdata <- rownames(bulk)
      pdata <- cbind(
        sampleID = seq_len(ncol(bulk)),
        subjectname = rep("num1", ncol(bulk)),
        celltypeID = seq_len(ncol(bulk))
      )
      SCDC::getESET(bulk, fdata = fdata, pdata = pdata)
    })
    set.seed(seed)
    run_stage("deconv_music", DECEPTICON_music_custom(bulk_eset, signature_paths))
  },
  finally = unlink(runtime_helpers)
)

manifest <- data.frame(
  key = c(
    "started_at_utc", "finished_at_utc", "seed", "rda_sha256", "bulk_sha256",
    "rda_genes", "rda_cells_after_exclusion", "bulk_genes", "bulk_samples",
    "decepticonx_git", "decepticon_git", "cibersort_helper_sha256",
    "music_prop_helper_sha256", "music_basis_helper_sha256",
    "epic_mrna_cell_sha256", "epic_mrna_cell_entries", "epic_version"
  ),
  value = c(
    started_at, format(Sys.time(), tz = "UTC", usetz = TRUE), seed,
    sha256(rda_path), sha256(bulk_path), nrow(single_cell), ncol(single_cell),
    nrow(bulk), ncol(bulk), git_head(decepticonx_r), git_head(decepticon_methods_r),
    sha256(file.path(helper_dir, "CIBERSORT.R")),
    sha256(file.path(helper_dir, "music_prop.R")),
    sha256(file.path(helper_dir, "music_basis.R")),
    epic_mrna_cell_sha256, epic_mrna_cell_entries,
    as.character(utils::packageVersion("EPIC"))
  ),
  stringsAsFactors = FALSE
)
utils::write.table(
  manifest, file.path(output_dir, "manifest.tsv"), sep = "\t", row.names = FALSE,
  quote = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
