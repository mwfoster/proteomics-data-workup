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

ambiguous_metadata <- rbind(metadata, transform(metadata[1L, , drop = FALSE], Sample = "run-3", SampleName = "Saos2 Nuc SEL & BRT 1-2"))
ambiguous <- resolve_proteomics_processed_sample_ids(
  "Saos2 Nuc SEL & BRT 1",
  ambiguous_metadata,
  ambiguous_metadata$SampleName
)
stopifnot(is.na(ambiguous))

cat("Restored header remapping tests passed.\n")
