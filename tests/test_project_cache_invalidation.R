source(if (file.exists("project_duckdb.R")) "project_duckdb.R" else "../project_duckdb.R")

metadata <- data.frame(Sample = c("run_1", "run_2"), SampleDetailsID = c("ID1", "ID2"), Condition = c("Control", "Case"), Batch = c("1", "1"), RunOrder = c(1, 2), stringsAsFactors = FALSE)
cache <- list(
  settings = list(workflow_tabs = "PCA"), metadata = metadata,
  spqc_metadata_edits = data.frame(Sample = "run_2", Column = "Batch", Value = "2"),
  processed_s2 = data.frame(Feature = "P1", run_1 = 10),
  processed_s3 = data.frame(Feature = "P1", run_1 = 11),
  processed_sample_map = data.frame(Sample = c("run_1", "run_2")),
  batch_corrected_s3 = data.frame(Feature = "P1", run_1 = 10.5),
  batch_corrected_sample_map = data.frame(Sample = c("run_1", "run_2")),
  statistics_tables = data.frame(Feature = "P1", PValue = 0.05)
)

metadata_changed <- invalidate_proteomics_project_cache(cache, "run_order_file")
stopifnot(identical(metadata_changed$settings, cache$settings))
stopifnot(identical(metadata_changed$metadata, cache$metadata))
stopifnot(identical(metadata_changed$processed_s2, cache$processed_s2))
stopifnot(identical(metadata_changed$processed_s3, cache$processed_s3))
stopifnot(is.null(metadata_changed$processed_sample_map))
stopifnot(is.null(metadata_changed$batch_corrected_s3))
stopifnot(is.null(metadata_changed$batch_corrected_sample_map))
stopifnot(is.null(metadata_changed$statistics_tables))

protein_changed <- invalidate_proteomics_project_cache(cache, "protein_imputed_file")
stopifnot(identical(protein_changed$metadata, cache$metadata))
stopifnot(identical(protein_changed$processed_s2, cache$processed_s2))
stopifnot(is.null(protein_changed$processed_s3))
stopifnot(is.null(protein_changed$processed_sample_map))
stopifnot(is.null(protein_changed$batch_corrected_s3))
stopifnot(is.null(protein_changed$statistics_tables))

replacement_run_order <- data.frame(`Run Label` = c("run_1", "run_2"), `#` = c(20, 10), stringsAsFactors = FALSE, check.names = FALSE)
run_updated <- merge_proteomics_metadata_replacement(metadata, replacement_run_order, "run_order_file")
stopifnot(identical(run_updated$RunOrder, c(20, 10)))
stopifnot(identical(run_updated$Condition, metadata$Condition))
stopifnot(identical(attr(run_updated, "overlay_audit")$matched, 2L))

replacement_details <- data.frame(ID = c("1", "2"), Batch = c("2", "3"), Cohort = c("A", "B"), stringsAsFactors = FALSE)
details_updated <- merge_proteomics_metadata_replacement(metadata, replacement_details, "sample_details_file")
stopifnot(identical(details_updated$Batch, c("2", "3")))
stopifnot(identical(details_updated$Cohort, c("A", "B")))
stopifnot(identical(details_updated$Condition, metadata$Condition))

replacement_metadata <- data.frame(Sample = c("run_2", "run_1"), Condition = c("Case updated", "Control updated"), stringsAsFactors = FALSE)
validated <- validate_proteomics_metadata_replacement(replacement_metadata, metadata$Sample)
stopifnot(identical(validated$Sample, replacement_metadata$Sample))
stopifnot(inherits(try(validate_proteomics_metadata_replacement(replacement_metadata[, -1, drop = FALSE], metadata$Sample), silent = TRUE), "try-error"))
stopifnot(inherits(try(validate_proteomics_metadata_replacement(rbind(replacement_metadata[1, ], replacement_metadata[1, ]), metadata$Sample), silent = TRUE), "try-error"))
stopifnot(inherits(try(validate_proteomics_metadata_replacement(replacement_metadata[1, ], metadata$Sample), silent = TRUE), "try-error"))

replacement_with_batch <- data.frame(
  Sample = c("run_1", "run_2"),
  Condition = c("Control", "Case"),
  Batch = c("1", "1"),
  stringsAsFactors = FALSE
)
spqc_edits <- data.frame(
  Sample = c("run_2", "run_2", "missing", "run_1"),
  Column = c("Batch", "NewColumn", "Batch", "Sample"),
  Value = c("3", "ignored", "4", "changed"),
  stringsAsFactors = FALSE
)
replacement_with_edits <- apply_proteomics_metadata_cell_edits(replacement_with_batch, spqc_edits)
# Removing the overlay would leave the downloaded replacement metadata at Batch 1.
stopifnot(identical(replacement_with_edits$Batch, c("1", "3")))
stopifnot(identical(replacement_with_edits$Sample, c("run_1", "run_2")))
stopifnot(!"NewColumn" %in% colnames(replacement_with_edits))

replacement_layers <- prepare_proteomics_metadata_replacement(replacement_with_batch)
# Keeping old edits here would silently overwrite the newly uploaded metadata.
stopifnot(identical(replacement_layers$metadata, replacement_with_batch))
stopifnot(identical(replacement_layers$spqc_metadata_edits, data.frame(
  Sample = character(0), Column = character(0), Value = character(0), stringsAsFactors = FALSE
)))

status_rows <- data.frame(FileType = c("Metadata", "Unused"), Status = c("Restored from DuckDB", "Not loaded"), stringsAsFactors = FALSE)
loaded_rows <- loaded_proteomics_project_status(status_rows)
stopifnot(nrow(loaded_rows) == 1L, identical(loaded_rows$FileType, "Metadata"))

message("Project cache invalidation and metadata overlay tests passed.")
