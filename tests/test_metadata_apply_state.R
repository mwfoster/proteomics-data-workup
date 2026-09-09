source(if (file.exists("project_duckdb.R")) "project_duckdb.R" else "../project_duckdb.R")
source(if (file.exists("metadata_state.R")) "metadata_state.R" else "../metadata_state.R")

metadata <- data.frame(
  Sample = c("SPQC_1", "Sample_1"),
  Condition = c("SPQC", "Control"),
  Batch = c("1", "1"),
  stringsAsFactors = FALSE
)
committed <- data.frame(Sample = "SPQC_1", Column = "Batch", Value = "1", stringsAsFactors = FALSE)
draft <- data.frame(Sample = "SPQC_1", Column = "Batch", Value = "2", stringsAsFactors = FALSE)

state <- build_proteomics_applied_state(metadata, c("Sample", "Condition"), committed, draft)
stopifnot(identical(state$metadata$Batch, c("2", "1")))
stopifnot(identical(state$columns, c("Sample", "Condition")))
stopifnot(nrow(state$spqc_metadata_edits) == 1L, state$spqc_metadata_edits$Value[[1L]] == "2")

before <- state
error <- try(build_proteomics_applied_state(metadata, "MissingColumn", committed, draft), silent = TRUE)
stopifnot(inherits(error, "try-error"), identical(state, before))

duplicate_metadata <- rbind(metadata, metadata[1, , drop = FALSE])
error <- try(build_proteomics_applied_state(duplicate_metadata, "Sample", committed, draft), silent = TRUE)
stopifnot(inherits(error, "try-error"))

cat("Metadata apply state tests passed.\n")
