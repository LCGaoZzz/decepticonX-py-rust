#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(Matrix))

usage <- function() {
  cat(
    "Usage: Rscript export_music2_official.R --data-dir DIR --output-dir DIR\n",
    "\n",
    "Export the pinned MuSiC2 R objects into portable benchmark tables and a\n",
    "cells-by-genes Matrix Market reference.\n",
    sep = ""
  )
}

parse_args <- function(args) {
  if ("--help" %in% args || "-h" %in% args) {
    usage()
    quit(status = 0)
  }
  allowed <- c("--data-dir", "--output-dir")
  values <- list()
  position <- 1
  while (position <= length(args)) {
    flag <- args[[position]]
    if (!(flag %in% allowed)) {
      stop("Unknown argument: ", flag, call. = FALSE)
    }
    if (position == length(args)) {
      stop("Missing value for ", flag, call. = FALSE)
    }
    values[[substring(flag, 3)]] <- args[[position + 1]]
    position <- position + 2
  }
  missing <- setdiff(c("data-dir", "output-dir"), names(values))
  if (length(missing) > 0) {
    stop(
      "Missing required argument(s): ",
      paste(paste0("--", missing), collapse = ", "),
      call. = FALSE
    )
  }
  values
}

options <- parse_args(commandArgs(trailingOnly = TRUE))
data_dir <- options[["data-dir"]]
out_dir <- options[["output-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_paths <- c(
  bulk = file.path(data_dir, "music2_bulk_eset.rds"),
  single_cell = file.path(data_dir, "music2_EMTABsce_healthy.rds"),
  truth = file.path(data_dir, "music2_true_proportion.RData")
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0) {
  stop(
    "Missing downloaded input(s): ", paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

bulk <- readRDS(input_paths[["bulk"]])
sce <- readRDS(input_paths[["single_cell"]])
truth_env <- new.env(parent = emptyenv())
loaded <- load(input_paths[["truth"]], envir = truth_env)
if (!identical(loaded, "prop_all")) {
  stop("Unexpected truth objects: ", paste(loaded, collapse = ", "))
}

target_types <- c("acinar", "alpha", "beta", "delta", "ductal", "gamma")
cell_metadata <- as.data.frame(colData(sce))
required_cell_columns <- c("cellType", "sampleID", "SubjectName")
missing_cell_columns <- setdiff(required_cell_columns, colnames(cell_metadata))
if (length(missing_cell_columns) > 0) {
  stop("Single-cell metadata is missing: ", paste(missing_cell_columns, collapse = ", "))
}
keep <- as.character(cell_metadata$cellType) %in% target_types
reference_counts <- assay(sce, "counts")[, keep, drop = FALSE]
cell_metadata <- cell_metadata[keep, , drop = FALSE]
if (anyDuplicated(rownames(reference_counts))) {
  stop("Reference gene symbols are not unique")
}

reference_sparse <- as(Matrix(t(reference_counts), sparse = TRUE), "dgCMatrix")
writeMM(reference_sparse, file.path(out_dir, "reference_cells_by_genes.mtx"))
write.table(
  data.frame(gene = rownames(reference_counts), check.names = FALSE),
  file.path(out_dir, "reference_genes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(
    cell_id = colnames(sce)[keep],
    donor = as.character(cell_metadata$sampleID),
    subject = as.character(cell_metadata$SubjectName),
    cell_type = as.character(cell_metadata$cellType),
    check.names = FALSE
  ),
  file.path(out_dir, "reference_obs.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

bulk_matrix <- exprs(bulk)
if (anyDuplicated(rownames(bulk_matrix))) {
  stop("Bulk gene symbols are not unique")
}
if (!identical(colnames(bulk_matrix), sampleNames(bulk))) {
  stop("Bulk sample names and ExpressionSet sampleNames differ")
}
if (!("group" %in% colnames(pData(bulk))) || !("sampleID" %in% colnames(pData(bulk)))) {
  stop("Bulk sample metadata must contain sampleID and group")
}
write.table(
  data.frame(gene = rownames(bulk_matrix), bulk_matrix, check.names = FALSE),
  file.path(out_dir, "bulk.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(
    sample = sampleNames(bulk),
    source_sample_id = as.character(pData(bulk)$sampleID),
    group = as.character(pData(bulk)$group),
    check.names = FALSE
  ),
  file.path(out_dir, "sample_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

truth_long <- truth_env$prop_all
required_truth_columns <- c("sampleID", "celltype", "proportion")
missing_truth_columns <- setdiff(required_truth_columns, colnames(truth_long))
if (length(missing_truth_columns) > 0) {
  stop("Truth table is missing: ", paste(missing_truth_columns, collapse = ", "))
}
truth_long$sampleID <- as.character(truth_long$sampleID)
truth_long$celltype <- as.character(truth_long$celltype)
truth_long$proportion <- as.numeric(truth_long$proportion)
if (!setequal(unique(truth_long$sampleID), sampleNames(bulk))) {
  stop("Truth and bulk sample IDs differ")
}
if (!setequal(unique(truth_long$celltype), target_types)) {
  stop("Truth cell types differ from the six-type benchmark contract")
}
truth <- matrix(
  NA_real_,
  nrow = ncol(bulk_matrix),
  ncol = length(target_types),
  dimnames = list(colnames(bulk_matrix), target_types)
)
for (row in seq_len(nrow(truth_long))) {
  truth[truth_long$sampleID[[row]], truth_long$celltype[[row]]] <-
    truth_long$proportion[[row]]
}
if (any(!is.finite(truth)) || any(truth < 0)) {
  stop("Truth matrix contains invalid values")
}
if (max(abs(rowSums(truth) - 1)) > 1e-8) {
  stop("Truth rows do not sum to one")
}
write.table(
  data.frame(sample = rownames(truth), truth, check.names = FALSE),
  file.path(out_dir, "truth.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

reference_values <- as.vector(reference_counts)
summary <- data.frame(
  item = c(
    "reference_genes", "reference_cells", "reference_donors",
    "reference_fractional_entries", "reference_max_fractional_deviation",
    "bulk_genes", "bulk_samples", "healthy_samples", "t2d_samples",
    "truth_cell_types", "truth_max_row_sum_error"
  ),
  value = c(
    nrow(reference_counts), ncol(reference_counts),
    length(unique(as.character(cell_metadata$sampleID))),
    sum(abs(reference_values - round(reference_values)) > 1e-8),
    max(abs(reference_values - round(reference_values))),
    nrow(bulk_matrix), ncol(bulk_matrix),
    sum(tolower(as.character(pData(bulk)$group)) == "healthy"),
    sum(tolower(as.character(pData(bulk)$group)) == "t2d"),
    ncol(truth), max(abs(rowSums(truth) - 1))
  )
)
write.table(
  summary, file.path(out_dir, "export_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Exported pinned MuSiC2 benchmark inputs to", out_dir, "\n")
print(summary)
cat("Reference cell types:\n")
print(table(cell_metadata$cellType))
cat("Bulk groups:\n")
print(table(pData(bulk)$group))
