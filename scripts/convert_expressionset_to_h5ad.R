#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
allow_normalized <- "--allow-normalized" %in% args
args <- args[args != "--allow-normalized"]

if (length(args) != 2L) {
  stop(
    "Usage: convert_expressionset_to_h5ad.R INPUT.rda OUTPUT.h5ad ",
    "[--allow-normalized]",
    call. = FALSE
  )
}

required <- c("Biobase", "Matrix", "reticulate")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

# Import Python before loading a potentially multi-gigabyte RDA so an
# incompatible reticulate/NumPy environment fails quickly.  Reticulate builds
# linked against the NumPy 1.x ABI should use a Python environment with
# NumPy < 2 for this one-time converter.
py_modules <- tryCatch(
  list(
    ad = reticulate::import("anndata", convert = FALSE),
    np = reticulate::import("numpy", convert = FALSE),
    sp = reticulate::import("scipy.sparse", convert = FALSE),
    pd = reticulate::import("pandas", convert = FALSE)
  ),
  error = function(error) {
    stop(
      "Could not import anndata/numpy/scipy/pandas through reticulate: ",
      conditionMessage(error),
      call. = FALSE
    )
  }
)

input <- normalizePath(args[[1]], mustWork = TRUE)
output <- normalizePath(dirname(args[[2]]), mustWork = TRUE)
output <- file.path(output, basename(args[[2]]))

loaded <- load(input, envir = environment())
if (length(loaded) != 1L) {
  stop("Expected exactly one object in the RDA; found: ", paste(loaded, collapse = ", "))
}
eset <- get(loaded[[1]], envir = environment())
if (!methods::is(eset, "ExpressionSet")) {
  stop("The RDA object is not a Biobase ExpressionSet", call. = FALSE)
}

x <- Biobase::exprs(eset)
if (anyNA(x) || any(!is.finite(x)) || any(x < 0)) {
  stop("Expression assay contains NA/Inf or negative values", call. = FALSE)
}
if (!length(x) || !nrow(x) || !ncol(x)) {
  stop("Expression assay must contain genes and cells", call. = FALSE)
}

# Check the full assay without allocating another full-sized dense matrix.
# Sampling can incorrectly classify a mostly-integer transformed assay.
integer_like <- TRUE
for (first in seq.int(1L, ncol(x), by = 64L)) {
  last <- min(first + 63L, ncol(x))
  block <- x[, first:last, drop = FALSE]
  if (any(abs(block - round(block)) > 1e-8)) {
    integer_like <- FALSE
    break
  }
}
if (!integer_like && !allow_normalized) {
  stop(
    "Expression assay is not integer-like raw counts. Re-run with ",
    "--allow-normalized only for an explicitly marked compatibility input.",
    call. = FALSE
  )
}

obs <- Biobase::pData(eset)
var <- Biobase::fData(eset)
rownames(obs) <- Biobase::sampleNames(eset)
rownames(var) <- Biobase::featureNames(eset)
if (!"gene_symbol" %in% colnames(var)) {
  var$gene_symbol <- rownames(var)
}

# Factors cross the reticulate boundary more reliably as strings. Numeric
# embeddings are restored below in obsm.
for (name in colnames(obs)) {
  if (is.factor(obs[[name]])) obs[[name]] <- as.character(obs[[name]])
}
for (name in colnames(var)) {
  if (is.factor(var[[name]])) var[[name]] <- as.character(var[[name]])
}

# AnnData is cells x genes. A sparse transpose keeps the on-disk h5ad compact;
# the conversion can still require several GB of peak memory for a large dense
# ExpressionSet, so run it in an environment with adequate RAM.
x_csr_source <- Matrix::Matrix(t(x), sparse = TRUE)
rm(x)
gc()

ad <- py_modules$ad
np <- py_modules$np
sp <- py_modules$sp
pd <- py_modules$pd

py_x <- sp$csr_matrix(reticulate::r_to_py(x_csr_source, convert = TRUE))
py_obs <- pd$DataFrame(reticulate::r_to_py(obs, convert = TRUE))
py_var <- pd$DataFrame(reticulate::r_to_py(var, convert = TRUE))
py_obs$index <- reticulate::r_to_py(rownames(obs))
py_var$index <- reticulate::r_to_py(rownames(var))

adata <- ad$AnnData(X = py_x, obs = py_obs, var = py_var)
if (all(c("tsneX1", "tsneX2") %in% colnames(obs))) {
  tsne <- cbind(as.numeric(obs$tsneX1), as.numeric(obs$tsneX2))
  adata$obsm[["X_tsne"]] <- np$array(tsne, dtype = "float64")
}
adata$uns[["source_format"]] <- "Biobase::ExpressionSet"
adata$uns[["source_object"]] <- loaded[[1]]
adata$uns[["expression_scale"]] <- if (integer_like) "raw_counts" else "normalized_noninteger"
invisible(adata$write_h5ad(output, compression = "gzip"))

cat("Wrote", output, "\n")
cat("shape:", nrow(obs), "cells x", nrow(var), "genes\n")
cat("expression_scale:", if (integer_like) "raw_counts" else "normalized_noninteger", "\n")
