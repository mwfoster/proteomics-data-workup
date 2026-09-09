source(if (file.exists("selection_helpers.R")) "selection_helpers.R" else "../selection_helpers.R")

hits <- data.frame(
  Protein = c("GENE2", "GENE1", "fallback"),
  ProteinGroupID = c("PG2", "PG1", "missing"),
  stringsAsFactors = FALSE
)
feature_info <- data.frame(
  Protein = c("GENE1", "GENE2", "fallback"),
  PG.ProteinGroups = c("PG1", "PG2", "PG3"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(identical(
  map_volcano_hits_to_features(hits, feature_info),
  c("GENE2", "GENE1", "fallback")
))
stopifnot(identical(volcano_source_to_feature_source("S2"), "no_impute"))
stopifnot(identical(volcano_source_to_feature_source("S3"), "imputed"))
stopifnot(identical(volcano_source_to_feature_source("S3_batch_corrected"), "S3_batch_corrected"))

expression <- data.frame(
  Feature = c("GENE1", "GENE2"),
  SampleA = c(10, 30),
  SampleB = c(20, 40),
  check.names = FALSE
)
long <- build_multifeature_boxplot_data(expression, c("GENE2", "GENE1"))
stopifnot(identical(as.character(long$Feature), c("GENE2", "GENE2", "GENE1", "GENE1")))
stopifnot(identical(as.character(long$Sample), c("SampleA", "SampleB", "SampleA", "SampleB")))
stopifnot(identical(long$Value, c(30, 40, 10, 20)))

md <- data.frame(Sample = c("A", "A", "B"), Condition = c("X", "X", "Y"), stringsAsFactors = FALSE)
plot_md <- unique_sample_metadata(md, c("A", "A", "B"))
stopifnot(identical(plot_md$Sample, c("A", "B")))
stopifnot(nrow(plot_md) == 2L)

comparison <- make_stats_comparison_id("Condition", "Case", "Control")
stopifnot(identical(comparison, "Condition|||Case|||Control"))
parsed <- parse_stats_comparison_id(comparison)
stopifnot(identical(parsed$group_col, "Condition"), identical(parsed$numerator, "Case"), identical(parsed$denominator, "Control"))
legacy <- parse_stats_comparison_id("Case|||Control")
stopifnot(isTRUE(legacy$legacy), identical(legacy$group_col, "Condition"))
stopifnot(identical(stats_comparison_label(comparison), "Condition: Case vs Control"))
stopifnot(identical(stats_comparison_prefix(comparison), "Condition_Case_vs_Control"))
stopifnot(identical(stats_comparison_prefix("Case|||Control"), "Case_vs_Control"))

app_text <- paste(readLines("app.R", warn = FALSE), collapse = "\n")
stopifnot(!grepl('updateSelectInput(session, "script_box_group_by", selected = "Condition")', app_text, fixed = TRUE))

message("Volcano feature-selection tests passed.")
