app_dir <- if (file.exists("app.R")) getwd() else normalizePath("..", mustWork = TRUE)
source(file.path(app_dir, "imputation_helpers.R"))

report <- data.frame(
  Protein = c("P1", "P2", "P3", "P4"),
  S1_Protein_group_abundance = c(2, 2, 32, NA),
  S2_Protein_group_abundance = c(4, 4, 32, NA),
  S3_Protein_group_abundance = c(NA, 8, 64, 16),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

result <- knn_impute_protein_report(
  report,
  abundance_columns = grep("_Protein_group_abundance$", colnames(report), value = TRUE),
  k = 1,
  max_missing_percent = 50
)

# Removing or breaking nearest-protein selection makes P1/S3 differ from 8.
stopifnot(identical(result$data$Protein, c("P1", "P2", "P3")))
stopifnot(isTRUE(all.equal(result$data$S3_Protein_group_abundance[1], 8, tolerance = 1e-10)))
# Changing observed values during imputation is data corruption.
stopifnot(identical(result$data$S1_Protein_group_abundance, c(2, 2, 32)))
stopifnot(identical(result$dropped_features, 1L))
stopifnot(identical(result$missing_before, 1L))
stopifnot(identical(result$missing_after, 0L))

grouped <- data.frame(
  Protein = c("P1", "P2", "P3"),
  A1_Protein_group_abundance = c(2, 2, 32),
  A2_Protein_group_abundance = c(NA, 4, 32),
  B1_Protein_group_abundance = c(128, 8, 64),
  B2_Protein_group_abundance = c(NA, 16, 64),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
group_columns <- list(
  A = c("A1_Protein_group_abundance", "A2_Protein_group_abundance"),
  B = c("B1_Protein_group_abundance", "B2_Protein_group_abundance")
)
group_result <- knn_impute_protein_report(
  grouped,
  abundance_columns = unlist(group_columns, use.names = FALSE),
  k = 1,
  max_missing_percent = 60,
  column_groups = group_columns
)

# Removing group isolation would let the very different B values affect A2.
stopifnot(isTRUE(all.equal(group_result$data$A2_Protein_group_abundance[1], 4, tolerance = 1e-10)))
stopifnot(isTRUE(all.equal(group_result$data$B2_Protein_group_abundance[1], 64, tolerance = 1e-10)))

uploaded_s3 <- data.frame(Protein = "uploaded", stringsAsFactors = FALSE)
cached_s3 <- data.frame(Protein = "cached", stringsAsFactors = FALSE)
generated_s3 <- data.frame(Protein = "knn", stringsAsFactors = FALSE)
s2 <- data.frame(Protein = "s2", stringsAsFactors = FALSE)

# These assertions catch S3 silently reverting to a Spectronaut upload.
stopifnot(identical(resolve_protein_s3_source("spectronaut", uploaded_s3, cached_s3, generated_s3, s2), uploaded_s3))
stopifnot(identical(resolve_protein_s3_source("spectronaut", NULL, cached_s3, generated_s3, s2), cached_s3))
stopifnot(identical(resolve_protein_s3_source("knn", uploaded_s3, cached_s3, generated_s3, s2), generated_s3))
stopifnot(identical(resolve_protein_s3_source("s2", uploaded_s3, cached_s3, generated_s3, s2), s2))

message("Protein kNN imputation tests passed.")
