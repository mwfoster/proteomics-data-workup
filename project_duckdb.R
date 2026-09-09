proteomics_duckdb_available <- function() {
  requireNamespace("DBI", quietly = TRUE) &&
    requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("jsonlite", quietly = TRUE)
}

empty_spqc_metadata_edits <- function() {
  data.frame(Sample = character(0), Column = character(0), Value = character(0), stringsAsFactors = FALSE)
}

proteomics_project_backend_for_path <- function(path) {
  if (identical(tolower(tools::file_ext(path)), "rds")) "rds" else "duckdb"
}

proteomics_project_path <- function(folder, filename) {
  folder <- trimws(as.character(folder)[1L])
  filename <- trimws(as.character(filename)[1L])
  if (is.na(folder) || !nzchar(folder)) stop("Select a project folder.")
  if (is.na(filename) || !nzchar(filename)) stop("Enter a project filename.")
  if (!identical(basename(filename), filename)) stop("Project filename must not include a folder path.")
  if (!grepl("\\.(duckdb|db|rds)$", filename, ignore.case = TRUE)) filename <- paste0(filename, ".duckdb")
  normalizePath(file.path(folder, filename), mustWork = FALSE)
}

proteomics_project_backend_available <- function(path) {
  if (identical(proteomics_project_backend_for_path(path), "rds")) {
    exists("proteomics_rds_available") && proteomics_rds_available()
  } else {
    proteomics_duckdb_available()
  }
}

proteomics_is_no_impute_source <- function(source) {
  source <- as.character(source)
  length(source) > 0L && !is.na(source[[1L]]) && source[[1L]] %in% c("no_impute", "S2")
}

resolve_proteomics_project_table <- function(uploaded_table = NULL, cache = list(), cache_name) {
  if (!is.null(uploaded_table)) {
    return(as.data.frame(uploaded_table, stringsAsFactors = FALSE, check.names = FALSE))
  }
  restored <- cache[[cache_name]]
  if (is.null(restored)) return(NULL)
  as.data.frame(restored, stringsAsFactors = FALSE, check.names = FALSE)
}

proteomics_is_upload_descriptor <- function(source) {
  is.data.frame(source) && nrow(source) == 1L &&
    all(c("name", "size", "type", "datapath") %in% colnames(source))
}

proteomics_upload_is_ignored <- function(file_info, ignored_datapath) {
  if (is.null(file_info) || is.null(ignored_datapath) || !length(ignored_datapath)) return(FALSE)
  datapath <- if (is.data.frame(file_info) && "datapath" %in% colnames(file_info)) {
    as.character(file_info$datapath[1L])
  } else if (!is.null(file_info$datapath)) {
    as.character(file_info$datapath[1L])
  } else ""
  ignored_datapath <- as.character(ignored_datapath[1L])
  !is.na(datapath) && nzchar(datapath) && !is.na(ignored_datapath) && identical(datapath, ignored_datapath)
}

empty_proteomics_project_payload <- function() {
  list(
    settings = list(),
    metadata = NULL,
    spqc_metadata_edits = data.frame(
      Sample = character(0), Column = character(0), Value = character(0),
      stringsAsFactors = FALSE
    ),
    source_file_manifest = NULL,
    processed_s2 = NULL,
    processed_s3 = NULL,
    processed_sample_map = NULL,
    batch_corrected_s3 = NULL,
    batch_corrected_sample_map = NULL,
    statistics_tables = NULL,
    identifications_overview = NULL,
    run_identifications_precursor = NULL,
    run_identifications_protein = NULL
  )
}

use_cached_proteomics_table <- function(source, refresh_revision) {
  revision <- suppressWarnings(as.integer(refresh_revision)[1L])
  is.data.frame(source) && nrow(source) > 0L && !proteomics_is_upload_descriptor(source) &&
    !is.na(revision) && revision <= 0L
}

proteomics_abundance_columns <- function(table) {
  names(table)[grepl("\\.PG\\.Quantity$|_Protein_group_abundance$", names(table))]
}

proteomics_abundance_sample_names <- function(columns) {
  samples <- sub("^\\[[0-9]+\\][[:space:]]*", "", columns)
  sub("\\.PG\\.Quantity$|_Protein_group_abundance$", "", samples)
}

proteomics_measurement_kind <- function(columns) {
  clean <- sub("^\\[[0-9]+\\][[:space:]]*", "", as.character(columns))
  precursor <- grepl("\\.PG\\.NrOfPrecursorsUsedForQuantification$|_quantified_precursors$", clean)
  abundance <- grepl("\\.PG\\.Quantity$|_Protein_group_abundance$", clean)
  ifelse(precursor, "precursor", ifelse(abundance, "abundance", NA_character_))
}

order_proteomics_measurement_columns <- function(columns, sample_rank = seq_along(columns)) {
  kind <- proteomics_measurement_kind(columns)
  sample_rank <- suppressWarnings(as.numeric(sample_rank))
  order_index <- order(match(kind, c("precursor", "abundance")), is.na(sample_rank), sample_rank, seq_along(columns))
  as.character(columns)[order_index]
}

proteomics_annotation_columns <- function(columns) {
  columns <- as.character(columns)
  measurement <- !is.na(proteomics_measurement_kind(columns))
  derived <- grepl("(_percent_CV|_log2_fold_change|_(paired|unpaired)_t_test_p_value|_BH_FDR)$", columns)
  columns[!(measurement | derived)]
}

restore_proteomics_sample_names <- function(columns, metadata_samples, metadata_header_labels) {
  labels <- proteomics_abundance_sample_names(columns)
  if (!any(endsWith(columns, "_Protein_group_abundance"))) return(labels)
  matched <- match(labels, metadata_header_labels)
  labels[!is.na(matched)] <- as.character(metadata_samples[matched[!is.na(matched)]])
  labels
}

invalidate_proteomics_project_cache <- function(cache, changed_input) {
  if (is.null(cache) || !length(cache)) return(cache)
  metadata_inputs <- c(
    "condition_setup_sample_details_file", "condition_setup_template_file",
    "meta_file", "run_order_file", "sample_details_file"
  )
  if (changed_input %in% metadata_inputs) {
    cache$processed_sample_map <- NULL
    cache$batch_corrected_s3 <- NULL
    cache$batch_corrected_sample_map <- NULL
    cache$statistics_tables <- NULL
    return(cache)
  }
  if (identical(changed_input, "protein_no_impute_file")) {
    cache$processed_s2 <- NULL
    cache$processed_sample_map <- NULL
    cache$statistics_tables <- NULL
  }
  if (identical(changed_input, "protein_imputed_file")) {
    cache$processed_s3 <- NULL
    cache$processed_sample_map <- NULL
    cache$batch_corrected_s3 <- NULL
    cache$batch_corrected_sample_map <- NULL
    cache$statistics_tables <- NULL
  }
  cache
}

loaded_proteomics_project_status <- function(status) {
  if (is.null(status) || !nrow(status) || !"Status" %in% colnames(status)) return(status)
  status[!is.na(status$Status) & status$Status != "Not loaded", , drop = FALSE]
}

validate_proteomics_metadata_replacement <- function(replacement, expected_samples) {
  replacement <- as.data.frame(replacement, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"Sample" %in% colnames(replacement)) stop("Modified metadata must contain a Sample column.", call. = FALSE)
  replacement$Sample <- trimws(as.character(replacement$Sample))
  if (any(is.na(replacement$Sample) | !nzchar(replacement$Sample))) stop("The Sample column cannot contain blank values.", call. = FALSE)
  if (anyDuplicated(replacement$Sample)) stop("The Sample column must contain one unique row per sample.", call. = FALSE)
  expected_samples <- trimws(as.character(expected_samples))
  if (!setequal(replacement$Sample, expected_samples)) {
    missing <- setdiff(expected_samples, replacement$Sample)
    unexpected <- setdiff(replacement$Sample, expected_samples)
    details <- c(
      if (length(missing)) paste0("missing: ", paste(missing, collapse = ", ")),
      if (length(unexpected)) paste0("unexpected: ", paste(unexpected, collapse = ", "))
    )
    stop(paste0("Modified metadata samples do not match the active project (", paste(details, collapse = "; "), ")."), call. = FALSE)
  }
  replacement
}

apply_proteomics_metadata_cell_edits <- function(metadata, edits) {
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
  edits <- as.data.frame(edits, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("Sample", "Column", "Value")
  if (!"Sample" %in% colnames(metadata) || !all(required %in% colnames(edits)) || !nrow(edits)) return(metadata)
  for (row_index in seq_len(nrow(edits))) {
    edit_sample <- as.character(edits$Sample[row_index])
    edit_column <- as.character(edits$Column[row_index])
    if (is.na(edit_column) || !nzchar(edit_column) || identical(edit_column, "Sample") || !edit_column %in% colnames(metadata)) next
    matched <- which(as.character(metadata$Sample) == edit_sample)
    if (length(matched)) metadata[[edit_column]][matched] <- as.character(edits$Value[row_index])
  }
  metadata
}

prepare_proteomics_metadata_replacement <- function(replacement) {
  list(
    metadata = as.data.frame(replacement, stringsAsFactors = FALSE, check.names = FALSE),
    spqc_metadata_edits = data.frame(
      Sample = character(0), Column = character(0), Value = character(0),
      stringsAsFactors = FALSE
    )
  )
}

merge_proteomics_metadata_replacement <- function(cached_metadata, replacement, changed_input) {
  cached_metadata <- as.data.frame(cached_metadata, stringsAsFactors = FALSE, check.names = FALSE)
  replacement <- as.data.frame(replacement, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(cached_metadata) || !nrow(replacement)) return(cached_metadata)
  if (identical(changed_input, "run_order_file")) {
    if (!"Sample" %in% colnames(replacement) && "Run Label" %in% colnames(replacement)) replacement$Sample <- as.character(replacement$`Run Label`)
    if ("#" %in% colnames(replacement)) replacement$RunOrder <- replacement$`#`
    cached_key <- replacement_key <- "Sample"
    skip_columns <- c("Sample", "#")
  } else if (identical(changed_input, "sample_details_file") || identical(changed_input, "condition_setup_sample_details_file")) {
    cached_key <- intersect(c("SampleDetailsID", "Sample ID", "SampleID", "ID"), colnames(cached_metadata))[1]
    replacement_key <- intersect(c("SampleDetailsID", "Sample ID", "SampleID", "ID"), colnames(replacement))[1]
    if (!is.na(replacement_key)) {
      replacement[[replacement_key]] <- ifelse(grepl("^ID", as.character(replacement[[replacement_key]]), ignore.case = TRUE), as.character(replacement[[replacement_key]]), paste0("ID", replacement[[replacement_key]]))
    }
    skip_columns <- replacement_key
  } else {
    if (!"Sample" %in% colnames(replacement)) {
      source_key <- intersect(c("Run Label", "File", "File Name", "Filename", "SampleName", "Sample Name"), colnames(replacement))[1]
      if (!is.na(source_key)) replacement$Sample <- as.character(replacement[[source_key]])
    }
    cached_key <- replacement_key <- "Sample"
    skip_columns <- "Sample"
  }
  valid_keys <- length(cached_key) && length(replacement_key) && !is.na(cached_key) && !is.na(replacement_key) && cached_key %in% colnames(cached_metadata) && replacement_key %in% colnames(replacement)
  if (!valid_keys) {
    attr(cached_metadata, "overlay_audit") <- list(matched = 0L, unmatched_cached = nrow(cached_metadata), unmatched_replacement = nrow(replacement))
    return(cached_metadata)
  }
  matched <- match(as.character(cached_metadata[[cached_key]]), as.character(replacement[[replacement_key]]))
  matched_rows <- !is.na(matched)
  for (column in setdiff(colnames(replacement), skip_columns)) {
    if (!column %in% colnames(cached_metadata)) cached_metadata[[column]] <- NA
    cached_metadata[[column]][matched_rows] <- replacement[[column]][matched[matched_rows]]
  }
  attr(cached_metadata, "overlay_audit") <- list(
    matched = as.integer(sum(matched_rows)),
    unmatched_cached = as.integer(sum(!matched_rows)),
    unmatched_replacement = as.integer(sum(!as.character(replacement[[replacement_key]]) %in% as.character(cached_metadata[[cached_key]])))
  )
  cached_metadata
}

proteomics_duckdb_connect <- function(path, read_only = FALSE) {
  if (!proteomics_duckdb_available()) {
    stop("Install R packages 'duckdb', 'DBI', and 'jsonlite' to use DuckDB projects.", call. = FALSE)
  }
  DBI::dbConnect(duckdb::duckdb(), dbdir = path, read_only = read_only)
}

proteomics_duckdb_write_optional <- function(connection, name, value) {
  if (is.null(value)) return(invisible(FALSE))
  DBI::dbWriteTable(connection, name, as.data.frame(value, stringsAsFactors = FALSE, check.names = FALSE), overwrite = TRUE)
  invisible(TRUE)
}

proteomics_duckdb_read_optional <- function(connection, name) {
  if (!DBI::dbExistsTable(connection, name)) return(NULL)
  DBI::dbReadTable(connection, name, check.names = FALSE)
}

save_proteomics_project_duckdb <- function(
  path, project_name = "proteomics-project", settings = list(), metadata = NULL,
  spqc_metadata_edits = NULL, source_file_manifest = NULL, processed_s2 = NULL,
  processed_s3 = NULL, processed_sample_map = NULL, batch_corrected_s3 = NULL,
  batch_corrected_sample_map = NULL, statistics_tables = NULL,
  identifications_overview = NULL, run_identifications_precursor = NULL,
  run_identifications_protein = NULL
) {
  if (!proteomics_duckdb_available()) {
    stop("Install R packages 'duckdb', 'DBI', and 'jsonlite' to use DuckDB projects.", call. = FALSE)
  }
  path <- normalizePath(path, mustWork = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile("proteomics_project_", tmpdir = dirname(path), fileext = ".duckdb")
  on.exit(unlink(temporary_path, force = TRUE), add = TRUE)
  connection <- proteomics_duckdb_connect(temporary_path)
  connected <- TRUE
  on.exit(if (connected) try(DBI::dbDisconnect(connection, shutdown = TRUE), silent = TRUE), add = TRUE)
  manifest <- data.frame(
    Key = c("format", "version", "project_name", "saved_at"),
    Value = c("Proteomics Data Workup project", "1", project_name, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    stringsAsFactors = FALSE
  )
  DBI::dbWriteTable(connection, "project_manifest", manifest, overwrite = TRUE)
  DBI::dbWriteTable(connection, "project_settings", data.frame(
    SettingsJSON = jsonlite::toJSON(settings, auto_unbox = TRUE, null = "null"),
    stringsAsFactors = FALSE
  ), overwrite = TRUE)
  objects <- list(
    metadata = metadata, spqc_metadata_edits = spqc_metadata_edits,
    source_file_manifest = source_file_manifest, processed_s2 = processed_s2,
    processed_s3 = processed_s3, processed_sample_map = processed_sample_map,
    batch_corrected_s3 = batch_corrected_s3,
    batch_corrected_sample_map = batch_corrected_sample_map,
    statistics_tables = statistics_tables,
    identifications_overview = identifications_overview,
    run_identifications_precursor = run_identifications_precursor,
    run_identifications_protein = run_identifications_protein
  )
  for (name in names(objects)) proteomics_duckdb_write_optional(connection, name, objects[[name]])
  DBI::dbDisconnect(connection, shutdown = TRUE)
  connected <- FALSE
  if (file.exists(path)) unlink(path, force = TRUE)
  moved <- file.rename(temporary_path, path)
  if (!moved) moved <- file.copy(temporary_path, path, overwrite = TRUE)
  if (!moved || !file.exists(path)) stop("Could not write project database: ", path, call. = FALSE)
  invisible(path)
}

load_proteomics_project_duckdb <- function(path) {
  connection <- proteomics_duckdb_connect(path, read_only = TRUE)
  on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
  settings_table <- proteomics_duckdb_read_optional(connection, "project_settings")
  settings <- if (!is.null(settings_table) && nrow(settings_table)) {
    jsonlite::fromJSON(settings_table$SettingsJSON[[1]], simplifyVector = FALSE)
  } else list()
  names_to_read <- c(
    "project_manifest", "metadata", "spqc_metadata_edits", "source_file_manifest",
    "processed_s2", "processed_s3", "processed_sample_map", "batch_corrected_s3",
    "batch_corrected_sample_map", "statistics_tables", "identifications_overview",
    "run_identifications_precursor", "run_identifications_protein"
  )
  values <- lapply(names_to_read, function(name) proteomics_duckdb_read_optional(connection, name))
  names(values) <- c("manifest", names_to_read[-1])
  c(list(settings = settings), values)
}

save_proteomics_project_file <- function(path, ...) {
  if (identical(proteomics_project_backend_for_path(path), "rds")) {
    save_proteomics_project_rds(path = path, ...)
  } else {
    save_proteomics_project_duckdb(path = path, ...)
  }
}

load_proteomics_project_file <- function(path) {
  if (identical(proteomics_project_backend_for_path(path), "rds")) {
    load_proteomics_project_rds(path)
  } else {
    load_proteomics_project_duckdb(path)
  }
}
