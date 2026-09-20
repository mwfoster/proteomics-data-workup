source(if (file.exists("selection_helpers.R")) "selection_helpers.R" else "../selection_helpers.R")

metadata <- data.frame(
  Sample = c("run-1", "run-2"),
  Name = c("Saos2\u00a0Nuc SEL &amp; BRT 1", "Saos2 Nuc SEL & BRT 2"),
  SampleName = c("Saos2 Nuc SEL & BRT 1-2", "Saos2 Nuc SEL & BRT 2-2"),
  Condition = c("Saos2 Nuc SEL BRT", "Saos2 Nuc SEL BRT"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

resolved <- resolve_proteomics_processed_sample_ids(
  c("Saos2 Nuc SEL & BRT 1", "Saos2 Nuc SEL & BRT 2"),
  metadata,
  metadata$SampleName
)
stopifnot(identical(unname(resolved), c("run-1", "run-2")))

saved_map <- data.frame(
  Sample = c("run-1", "run-2"),
  HeaderLabel = c("old composite label A", "old composite label B"),
  stringsAsFactors = FALSE
)
resolved_saved <- resolve_proteomics_processed_sample_ids(
  c("old composite label A", "not present"),
  metadata,
  metadata$SampleName,
  saved_map
)
stopifnot(identical(unname(resolved_saved), c("run-1", NA_character_)))

resolved_syntactic <- resolve_proteomics_processed_sample_ids(
  c("old.composite.label.A", "old.composite.label.B"),
  metadata,
  metadata$SampleName,
  saved_map
)
stopifnot(identical(unname(resolved_syntactic), c("run-1", "run-2")))

ambiguous_metadata <- rbind(metadata, transform(metadata[1L, , drop = FALSE], Sample = "run-3", SampleName = "Saos2 Nuc SEL & BRT 1-2"))
ambiguous <- resolve_proteomics_processed_sample_ids(
  "Saos2 Nuc SEL & BRT 1",
  ambiguous_metadata,
  ambiguous_metadata$SampleName
)
stopifnot(is.na(ambiguous))

restored_table <- data.frame(
  PG.ProteinGroups = "P1",
  X9_FA_quantified_precursors = 12,
  X9_FA_Protein_group_abundance = 100,
  X9_FA_percent_CV = 7,
  check.names = FALSE
)
historical_map <- data.frame(
  Sample = c("run-1", "run-2"),
  HeaderLabel = c("X9_FA", "X10_O3"),
  stringsAsFactors = FALSE
)
remapped <- remap_proteomics_processed_headers(
  restored_table,
  metadata,
  current_header_labels = metadata$SampleName,
  saved_sample_map = historical_map
)
stopifnot(identical(
  colnames(remapped),
  c(
    "PG.ProteinGroups",
    "Saos2 Nuc SEL & BRT 1-2_quantified_precursors",
    "Saos2 Nuc SEL & BRT 1-2_Protein_group_abundance",
    "X9_FA_percent_CV"
  )
))

app_text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")
stopifnot(grepl("return(remap_restored_protein_headers(source))", app_text, fixed = TRUE))
stopifnot(length(gregexpr("return(remap_restored_protein_headers(source))", app_text, fixed = TRUE)[[1L]]) == 2L)

cat("Restored header remapping tests passed.\n")
