source(if (file.exists("selection_helpers.R")) "selection_helpers.R" else "../selection_helpers.R")

metadata <- data.frame(
  Sample = c("run1", "run2", "run3"),
  Subject = c("S01", "S02", ""),
  Group = c("Control", "Case", "Case"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

labels <- compose_proteomics_metadata_labels(
  metadata,
  c("Subject", "Group"),
  fallback = metadata$Sample
)
stopifnot(identical(labels, c("S01_Control", "S02_Case", "Case")))

fallback_labels <- compose_proteomics_metadata_labels(
  metadata,
  character(0),
  fallback = metadata$Sample
)
stopifnot(identical(fallback_labels, c("run1", "run2", "run3")))

stopifnot(identical(feature_facet_label_size(12, 1), 12))
stopifnot(isTRUE(all.equal(feature_facet_label_size(12, 3), 9.6)))

cat("Feature axis-label tests passed.\n")
