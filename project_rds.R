proteomics_rds_available <- function() TRUE

save_proteomics_project_rds <- function(
  path, project_name = "proteomics-project", settings = list(), metadata = NULL,
  spqc_metadata_edits = NULL, source_file_manifest = NULL, processed_s2 = NULL,
  processed_s3 = NULL, processed_sample_map = NULL, batch_corrected_s3 = NULL,
  batch_corrected_sample_map = NULL, statistics_tables = NULL,
  identifications_overview = NULL, run_identifications_precursor = NULL,
  run_identifications_protein = NULL
) {
  path <- normalizePath(path, mustWork = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile("proteomics_project_", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(temporary_path, force = TRUE), add = TRUE)
  project <- list(
    manifest = data.frame(
      Key = c("format", "version", "project_name", "saved_at"),
      Value = c("Proteomics Data Workup project", "1", project_name, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      stringsAsFactors = FALSE
    ),
    settings = settings, metadata = metadata, spqc_metadata_edits = spqc_metadata_edits,
    source_file_manifest = source_file_manifest, processed_s2 = processed_s2,
    processed_s3 = processed_s3, processed_sample_map = processed_sample_map,
    batch_corrected_s3 = batch_corrected_s3,
    batch_corrected_sample_map = batch_corrected_sample_map,
    statistics_tables = statistics_tables,
    identifications_overview = identifications_overview,
    run_identifications_precursor = run_identifications_precursor,
    run_identifications_protein = run_identifications_protein
  )
  saveRDS(project, temporary_path, version = 3)
  invisible(readRDS(temporary_path))
  if (file.exists(path)) unlink(path, force = TRUE)
  moved <- file.rename(temporary_path, path)
  if (!moved) moved <- file.copy(temporary_path, path, overwrite = TRUE)
  if (!moved || !file.exists(path)) stop("Could not write RDS project: ", path, call. = FALSE)
  invisible(path)
}

load_proteomics_project_rds <- function(path) {
  project <- readRDS(path)
  if (!is.list(project) || is.null(project$manifest)) stop("Invalid Proteomics Data Workup RDS project.", call. = FALSE)
  project
}
