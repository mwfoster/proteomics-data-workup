source(if (file.exists("project_duckdb.R")) "project_duckdb.R" else "../project_duckdb.R")

restored_s3 <- data.frame(
  Protein = c("P1", "P2"),
  A_Protein_group_abundance = c(10, 20),
  B_Protein_group_abundance = c(15, 25),
  check.names = FALSE
)

stopifnot(!proteomics_is_no_impute_source(NULL))
stopifnot(!proteomics_is_no_impute_source(character(0)))
stopifnot(!proteomics_is_no_impute_source(NA_character_))
stopifnot(proteomics_is_no_impute_source("S2"))
stopifnot(proteomics_is_no_impute_source("no_impute"))

resolved <- resolve_proteomics_project_table(
  uploaded_table = NULL,
  cache = list(processed_s3 = restored_s3),
  cache_name = "processed_s3"
)
stopifnot(identical(resolved, restored_s3))

uploaded_s3 <- data.frame(Protein = "uploaded", A.PG.Quantity = 99, check.names = FALSE)
resolved_upload <- resolve_proteomics_project_table(
  uploaded_table = uploaded_s3,
  cache = list(processed_s3 = restored_s3),
  cache_name = "processed_s3"
)
stopifnot(identical(resolved_upload, uploaded_s3))

stopifnot(use_cached_proteomics_table(restored_s3, 0L))
stopifnot(!use_cached_proteomics_table(restored_s3, 1L))
stopifnot(!use_cached_proteomics_table(NULL, 0L))
upload_descriptor <- data.frame(
  name = "report.tsv", size = 123, type = "text/tab-separated-values",
  datapath = tempfile(fileext = ".tsv"), stringsAsFactors = FALSE
)
stopifnot(proteomics_is_upload_descriptor(upload_descriptor))
# Treating this descriptor as a cached protein table caused a one-row S2 to be saved.
stopifnot(!use_cached_proteomics_table(upload_descriptor, 0L))
stopifnot(proteomics_upload_is_ignored(upload_descriptor, upload_descriptor$datapath))
stopifnot(!proteomics_upload_is_ignored(upload_descriptor, paste0(upload_descriptor$datapath, "_new")))

empty_project <- empty_proteomics_project_payload()
# A new project must never inherit any table or metadata from the active session.
stopifnot(identical(empty_project$settings, list()))
stopifnot(is.null(empty_project$metadata), is.null(empty_project$processed_s2), is.null(empty_project$processed_s3))
stopifnot(identical(empty_project$spqc_metadata_edits, data.frame(
  Sample = character(0), Column = character(0), Value = character(0), stringsAsFactors = FALSE
)))

stopifnot(identical(
  proteomics_abundance_columns(restored_s3),
  c("A_Protein_group_abundance", "B_Protein_group_abundance")
))
stopifnot(identical(
  proteomics_abundance_sample_names(c("[1] RunA.PG.Quantity", "RunB_Protein_group_abundance")),
  c("RunA", "RunB")
))
stopifnot(identical(
  restore_proteomics_sample_names(
    c("9_FA_Protein_group_abundance", "rep1_SPQC_Protein_group_abundance"),
    metadata_samples = c("raw_sample_9", "raw_spqc_1"),
    metadata_header_labels = c("9_FA", "rep1_SPQC")
  ),
  c("raw_sample_9", "raw_spqc_1")
))

measurement_columns <- c(
  "FA_9_quantified_precursors", "FA_9_Protein_group_abundance",
  "O3_9_quantified_precursors", "O3_9_Protein_group_abundance"
)
stopifnot(identical(
  proteomics_measurement_kind(measurement_columns),
  c("precursor", "abundance", "precursor", "abundance")
))
stopifnot(identical(
  order_proteomics_measurement_columns(measurement_columns, c(1L, 1L, 2L, 2L)),
  c("FA_9_quantified_precursors", "O3_9_quantified_precursors", "FA_9_Protein_group_abundance", "O3_9_Protein_group_abundance")
))
stopifnot(identical(
  proteomics_annotation_columns(c(
    "PG.Genes", "SPQC_percent_CV", "Group_O3_vs_FA_log2_fold_change",
    measurement_columns
  )),
  "PG.Genes"
))

app_text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")
stopifnot(!grepl('if (is.null(project_file("protein_no_impute_file")) && !is.null(cache$processed_s2)) return', app_text, fixed = TRUE))
stopifnot(!grepl('if (is.null(project_file("protein_imputed_file")) && !is.null(cache$processed_s3)) return', app_text, fixed = TRUE))

message("Restored project table tests passed.")
