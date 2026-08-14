# Generate the checked-in selector oracle using base R only.
#
# Usage:
#   Rscript tests/oracles/generate_r_selector_oracle.R [output.tsv]
#
# The implementation mirrors DECEPTICON's row-by-row flattening, stable
# decreasing order, and selection of positions 1 and 3 among the top four
# symmetric correlation entries.  It intentionally tests the selector only;
# DECEPTICON_custom_output's 15-strategy normalization block is undefined.

optimal_pairs <- function(matrix, n = 2L) {
  flattened_by_row <- as.numeric(t(matrix))
  ordered <- order(flattened_by_row, decreasing = TRUE)
  selected <- ordered[seq.int(1L, 2L * n, by = 2L)]
  cbind(
    left = ceiling(selected / nrow(matrix)),
    right = ((selected - 1L) %% nrow(matrix)) + 1L
  )
}

cases <- list(
  top_pairs = matrix(
    c(
      0.0, 0.9, 0.1, 0.8,
      0.9, 0.0, 0.7, 0.2,
      0.1, 0.7, 0.0, 0.95,
      0.8, 0.2, 0.95, 0.0
    ),
    nrow = 4L,
    byrow = TRUE
  ),
  masked_negative = matrix(
    c(
      0.0, -0.5, -0.5,
      -0.5, 0.0, -0.5,
      -0.5, -0.5, 0.0
    ),
    nrow = 3L,
    byrow = TRUE
  )
)

rows <- do.call(
  rbind,
  lapply(names(cases), function(case_name) {
    pairs <- optimal_pairs(cases[[case_name]], 2L)
    data.frame(
      case = case_name,
      pair_rank = seq_len(nrow(pairs)),
      left = pairs[, "left"],
      right = pairs[, "right"],
      stringsAsFactors = FALSE
    )
  })
)

arguments <- commandArgs(trailingOnly = TRUE)
destination <- if (length(arguments)) arguments[[1L]] else ""
write.table(
  rows,
  file = destination,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
