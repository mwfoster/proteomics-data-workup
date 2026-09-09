protein_knn_numeric_matrix <- function(report, abundance_columns) {
  abundance_columns <- intersect(as.character(abundance_columns), colnames(report))
  if (!length(abundance_columns)) stop("No protein abundance columns were found for kNN imputation.", call. = FALSE)
  values <- as.data.frame(report[, abundance_columns, drop = FALSE], stringsAsFactors = FALSE, check.names = FALSE)
  values[] <- lapply(values, function(column) suppressWarnings(as.numeric(as.character(column))))
  matrix_values <- as.matrix(values)
  matrix_values[!is.finite(matrix_values) | matrix_values <= 0] <- NA_real_
  matrix_values
}

protein_knn_fill_log2 <- function(log2_matrix, k) {
  output <- log2_matrix
  missing_cells <- which(is.na(log2_matrix), arr.ind = TRUE)
  if (!nrow(missing_cells)) return(output)
  k <- max(1L, suppressWarnings(as.integer(k)[1L]))

  for (cell_index in seq_len(nrow(missing_cells))) {
    row_index <- missing_cells[cell_index, 1L]
    column_index <- missing_cells[cell_index, 2L]
    candidate_rows <- which(!is.na(log2_matrix[, column_index]))
    candidate_rows <- setdiff(candidate_rows, row_index)
    distances <- vapply(candidate_rows, function(candidate_index) {
      comparison_columns <- setdiff(seq_len(ncol(log2_matrix)), column_index)
      shared <- comparison_columns[
        is.finite(log2_matrix[row_index, comparison_columns]) &
          is.finite(log2_matrix[candidate_index, comparison_columns])
      ]
      if (!length(shared)) return(Inf)
      sqrt(mean((log2_matrix[row_index, shared] - log2_matrix[candidate_index, shared])^2))
    }, numeric(1))
    usable <- which(is.finite(distances))

    if (length(usable)) {
      usable <- usable[order(distances[usable], candidate_rows[usable])]
      usable <- usable[seq_len(min(k, length(usable)))]
      chosen_rows <- candidate_rows[usable]
      chosen_distances <- distances[usable]
      zero_distance <- chosen_distances <= sqrt(.Machine$double.eps)
      estimate <- if (any(zero_distance)) {
        mean(log2_matrix[chosen_rows[zero_distance], column_index])
      } else {
        weights <- 1 / chosen_distances
        stats::weighted.mean(log2_matrix[chosen_rows, column_index], weights)
      }
    } else {
      estimate <- NA_real_
    }

    if (!is.finite(estimate)) estimate <- stats::median(log2_matrix[row_index, ], na.rm = TRUE)
    if (!is.finite(estimate)) estimate <- stats::median(log2_matrix[, column_index], na.rm = TRUE)
    if (!is.finite(estimate)) estimate <- stats::median(log2_matrix, na.rm = TRUE)
    output[row_index, column_index] <- estimate
  }
  output
}

knn_impute_protein_report <- function(report, abundance_columns, k = 10L,
                                      max_missing_percent = 50,
                                      column_groups = NULL) {
  report <- as.data.frame(report, stringsAsFactors = FALSE, check.names = FALSE)
  abundance_columns <- intersect(as.character(abundance_columns), colnames(report))
  matrix_values <- protein_knn_numeric_matrix(report, abundance_columns)
  max_missing_percent <- suppressWarnings(as.numeric(max_missing_percent)[1L])
  if (!is.finite(max_missing_percent) || max_missing_percent < 0 || max_missing_percent > 100) {
    stop("Maximum missingness must be between 0 and 100 percent.", call. = FALSE)
  }
  missing_percent <- 100 * rowMeans(is.na(matrix_values))
  keep <- missing_percent <= max_missing_percent
  dropped_features <- sum(!keep)
  report <- report[keep, , drop = FALSE]
  matrix_values <- matrix_values[keep, , drop = FALSE]
  missing_before <- sum(is.na(matrix_values))

  if (is.null(column_groups)) column_groups <- list(All_samples = abundance_columns)
  column_groups <- lapply(column_groups, function(columns) intersect(as.character(columns), abundance_columns))
  column_groups <- column_groups[lengths(column_groups) > 0L]
  covered <- unique(unlist(column_groups, use.names = FALSE))
  if (!setequal(covered, abundance_columns)) {
    stop("Every abundance column must belong to an imputation group.", call. = FALSE)
  }

  imputed <- matrix_values
  for (columns in column_groups) {
    indexes <- match(columns, abundance_columns)
    group_log2 <- log2(matrix_values[, indexes, drop = FALSE])
    imputed[, indexes] <- 2^protein_knn_fill_log2(group_log2, k)
  }
  report[, abundance_columns] <- as.data.frame(imputed, check.names = FALSE)

  list(
    data = report,
    dropped_features = as.integer(dropped_features),
    missing_before = as.integer(missing_before),
    missing_after = as.integer(sum(is.na(imputed))),
    k = max(1L, suppressWarnings(as.integer(k)[1L])),
    max_missing_percent = max_missing_percent
  )
}

resolve_protein_s3_source <- function(method, uploaded_s3 = NULL, cached_s3 = NULL,
                                      generated_s3 = NULL, s2 = NULL) {
  method <- as.character(method)[1L]
  if (is.na(method) || !nzchar(method)) method <- "spectronaut"
  switch(
    method,
    knn = if (!is.null(generated_s3)) generated_s3 else cached_s3,
    s2 = s2,
    if (!is.null(uploaded_s3)) uploaded_s3 else cached_s3
  )
}
