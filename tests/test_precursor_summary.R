source(if (file.exists("project_duckdb.R")) "project_duckdb.R" else "../project_duckdb.R")

report <- data.frame(
  PG.ProteinGroups = c("P1", "P2", "P3"),
  A_quantified_precursors = c(1, 0, NA),
  B_quantified_precursors = c(3, NA, NA),
  C_quantified_precursors = c(99, 99, 99),
  check.names = FALSE
)
counts <- proteomics_precursor_count_matrix(
  report,
  precursor_sample_ids = c("A", "B", "C"),
  samples = c("B", "A", "D"),
  rows = c(1L, 2L, 3L)
)
stopifnot(identical(dim(counts), c(3L, 3L)))
stopifnot(identical(colnames(counts), c("B", "A", "D")))
stopifnot(identical(unname(counts[1, ]), c(3, 1, NA_real_)))
stopifnot(identical(unname(counts[2, ]), c(NA_real_, 0, NA_real_)))
stopifnot(identical(proteomics_mean_precursor_counts(counts), c(2, 0, NA_real_)))

table <- data.frame(
  PG.ProteinGroups = c("P1", "P2"),
  SPQC_percent_CV = c(10, 20),
  contrast_log2_fold_change = c(2, 3),
  A_Protein_group_abundance = c(100, 200),
  check.names = FALSE
)
inserted <- proteomics_insert_precursor_average(table, c(2, 0))
stopifnot(identical(colnames(inserted), c(
  "PG.ProteinGroups", "Average_quantified_precursor_count", "SPQC_percent_CV",
  "contrast_log2_fold_change", "A_Protein_group_abundance"
)))
stopifnot(identical(inserted$Average_quantified_precursor_count, c(2, 0)))
stopifnot(is.na(proteomics_measurement_kind("Average_quantified_precursor_count")))
stopifnot(identical(proteomics_annotation_columns(colnames(inserted)), "PG.ProteinGroups"))

without_cv <- table[, -2, drop = FALSE]
inserted_without_cv <- proteomics_insert_precursor_average(without_cv, c(2, 0))
stopifnot(identical(colnames(inserted_without_cv)[1:3], c(
  "PG.ProteinGroups", "Average_quantified_precursor_count", "contrast_log2_fold_change"
)))

cat("Precursor summary tests passed.\n")
