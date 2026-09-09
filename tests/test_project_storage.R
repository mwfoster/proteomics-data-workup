app_dir <- if (file.exists("project_duckdb.R")) getwd() else normalizePath("..", mustWork = TRUE)
source(file.path(app_dir, "project_duckdb.R"))
source(file.path(app_dir, "project_rds.R"))

settings <- list(workflow_tabs = "PCA", cv_conditions = c("Control", "SPQC"))
metadata <- data.frame(Sample = c("S1", "S2"), Condition = c("Control", "Case"), Batch = c("1", "2"), stringsAsFactors = FALSE)
spqc_edits <- data.frame(Sample = "S2", Column = "Batch", Value = "3", stringsAsFactors = FALSE)
processed_s2 <- data.frame(Feature = c("P1", "P2"), S1 = c(10, 20), S2 = c(11, 21), stringsAsFactors = FALSE)
processed_s3 <- data.frame(Feature = c("P1", "P2"), S1 = c(12, 22), S2 = c(13, 23), stringsAsFactors = FALSE)
sample_map <- data.frame(Column = c("S1", "S2"), Sample = c("S1", "S2"), stringsAsFactors = FALSE)
corrected <- data.frame(Feature = c("P1", "P2"), S1 = c(10.5, 20.5), S2 = c(10.5, 20.5), stringsAsFactors = FALSE)
statistics <- data.frame(Feature = c("P1", "P2"), PValue = c(0.01, 0.2), stringsAsFactors = FALSE)
identifications_overview <- data.frame(Condition = c("Case", "Control"), Replicate = c("1", "1"), Precursors = c(100, 90), ProteinGroups = c(50, 45), stringsAsFactors = FALSE)
run_identifications_precursor <- data.frame(SourceLabel = c("Case.1", "Control.1"), Complete = c(80, 70), Shared = c(10, 10), Sparse = c(5, 5), Unique = c(5, 5), stringsAsFactors = FALSE)
run_identifications_protein <- data.frame(SourceLabel = c("Case.1", "Control.1"), Complete = c(40, 35), Shared = c(5, 5), Sparse = c(3, 3), Unique = c(2, 2), stringsAsFactors = FALSE)

assert_round_trip <- function(path) {
  save_proteomics_project_file(
    path = path,
    project_name = "storage-test",
    settings = settings,
    metadata = metadata,
    spqc_metadata_edits = spqc_edits,
    source_file_manifest = data.frame(InputID = "protein_imputed_file", Name = "s3.csv"),
    processed_s2 = processed_s2,
    processed_s3 = processed_s3,
    processed_sample_map = sample_map,
    batch_corrected_s3 = corrected,
    batch_corrected_sample_map = sample_map,
    statistics_tables = statistics,
    identifications_overview = identifications_overview,
    run_identifications_precursor = run_identifications_precursor,
    run_identifications_protein = run_identifications_protein
  )
  loaded <- load_proteomics_project_file(path)
  stopifnot(identical(unlist(loaded$settings$cv_conditions, use.names = FALSE), c("Control", "SPQC")))
  stopifnot(identical(as.data.frame(loaded$metadata, stringsAsFactors = FALSE), metadata))
  stopifnot(identical(as.data.frame(loaded$spqc_metadata_edits, stringsAsFactors = FALSE), spqc_edits))
  stopifnot(identical(as.data.frame(loaded$processed_s2, stringsAsFactors = FALSE), processed_s2))
  stopifnot(identical(as.data.frame(loaded$processed_s3, stringsAsFactors = FALSE), processed_s3))
  stopifnot(identical(as.data.frame(loaded$processed_sample_map, stringsAsFactors = FALSE), sample_map))
  stopifnot(identical(as.data.frame(loaded$batch_corrected_s3, stringsAsFactors = FALSE), corrected))
  stopifnot(identical(as.data.frame(loaded$statistics_tables, stringsAsFactors = FALSE), statistics))
  stopifnot(identical(as.data.frame(loaded$identifications_overview, stringsAsFactors = FALSE), identifications_overview))
  stopifnot(identical(as.data.frame(loaded$run_identifications_precursor, stringsAsFactors = FALSE), run_identifications_precursor))
  stopifnot(identical(as.data.frame(loaded$run_identifications_protein, stringsAsFactors = FALSE), run_identifications_protein))
}

assert_round_trip(tempfile(fileext = ".rds"))

if (proteomics_duckdb_available()) {
  assert_round_trip(tempfile(fileext = ".duckdb"))
} else {
  message("Skipping DuckDB round trip; duckdb, DBI, and jsonlite are required.")
}

message("Project storage round-trip test passed.")
