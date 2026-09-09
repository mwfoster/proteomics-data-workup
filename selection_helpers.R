volcano_source_to_feature_source <- function(source) {
  source <- as.character(source)[1L]
  if (is.na(source) || !nzchar(source)) return("imputed")
  if (identical(source, "S2")) return("no_impute")
  if (identical(source, "S3")) return("imputed")
  source
}

map_volcano_hits_to_features <- function(hits, feature_info) {
  if (is.null(hits) || !nrow(hits) || is.null(feature_info) || !nrow(feature_info)) return(character(0))
  mapped <- rep(NA_character_, nrow(hits))
  if ("ProteinGroupID" %in% colnames(hits) && "PG.ProteinGroups" %in% colnames(feature_info)) {
    matched <- match(as.character(hits$ProteinGroupID), as.character(feature_info$PG.ProteinGroups))
    mapped[!is.na(matched)] <- as.character(feature_info$Protein[matched[!is.na(matched)]])
  }
  missing <- is.na(mapped) | !nzchar(mapped)
  if (any(missing) && "Protein" %in% colnames(hits) && "Protein" %in% colnames(feature_info)) {
    matched <- match(as.character(hits$Protein[missing]), as.character(feature_info$Protein))
    mapped[which(missing)[!is.na(matched)]] <- as.character(feature_info$Protein[matched[!is.na(matched)]])
  }
  unique(mapped[!is.na(mapped) & nzchar(mapped)])
}

build_multifeature_boxplot_data <- function(expression, selected_features) {
  selected_features <- intersect(as.character(selected_features), as.character(expression$Feature))
  if (!length(selected_features)) return(data.frame(Feature = character(0), Sample = character(0), Value = numeric(0)))
  rows <- lapply(selected_features, function(feature) {
    row <- expression[expression$Feature == feature, , drop = FALSE]
    data.frame(
      Feature = feature,
      Sample = colnames(expression)[-1L],
      Value = suppressWarnings(as.numeric(row[1L, -1L, drop = TRUE])),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

deduplicate_feature_samples <- function(data) {
  if (is.null(data) || !nrow(data)) return(data)
  required <- c("Feature", "Sample")
  if (!all(required %in% colnames(data))) return(data)
  key <- paste(as.character(data$Feature), as.character(data$Sample), sep = "\r")
  data[!duplicated(key), , drop = FALSE]
}

unique_sample_metadata <- function(metadata, samples = NULL) {
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"Sample" %in% colnames(metadata)) return(data.frame(Sample = unique(as.character(samples)), stringsAsFactors = FALSE))
  metadata$Sample <- as.character(metadata$Sample)
  metadata <- metadata[!is.na(metadata$Sample) & nzchar(metadata$Sample), , drop = FALSE]
  metadata <- metadata[!duplicated(metadata$Sample), , drop = FALSE]
  if (!is.null(samples)) {
    wanted <- unique(as.character(samples))
    metadata <- metadata[match(wanted, metadata$Sample), , drop = FALSE]
    metadata$Sample <- wanted
  }
  rownames(metadata) <- NULL
  metadata
}

cv_metadata_group_columns <- function(metadata) {
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
  excluded <- c("Sample", "SampleName", "AnalysisLabel", "Excluded", "ExclusionReason")
  candidates <- setdiff(colnames(metadata), excluded)
  candidates[vapply(candidates, function(column) {
    values <- trimws(as.character(metadata[[column]]))
    values <- values[!is.na(values) & nzchar(values)]
    groups <- table(values)
    length(groups) >= 2L && any(groups >= 2L) && length(groups) <= max(50L, floor(nrow(metadata) / 2L))
  }, logical(1))]
}

table_s1_metadata_choices <- function(available, selected) {
  selected <- unique(as.character(unlist(selected, use.names = FALSE)))
  selected[!is.na(selected) & selected %in% available]
}

retain_metadata_choice <- function(current, choices, default = character(0)) {
  current <- as.character(current)
  if (length(current)) return(if (!is.na(current[1L]) && current[1L] %in% choices) current[1L] else "")
  default <- default[default %in% choices]
  if (length(default)) default[1L] else ""
}

proteomics_multiselect_options <- function(placeholder = NULL, drag = FALSE, ...) {
  plugins <- c(if (isTRUE(drag)) "drag_drop", "remove_button")
  options <- list(plugins = as.list(plugins))
  if (!is.null(placeholder)) options$placeholder <- placeholder
  c(options, list(...))
}

stats_comparison_selectize_options <- function() {
  proteomics_multiselect_options(
    "Type a condition or metadata field to find comparisons",
    drag = TRUE,
    searchField = c("text", "value"),
    closeAfterSelect = FALSE,
    maxOptions = 10000
  )
}

stats_metadata_group_columns <- function(metadata) {
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
  excluded <- c("Sample", "SampleName", "AnalysisLabel", "Excluded", "ExclusionReason")
  candidates <- setdiff(colnames(metadata), excluded)
  candidates[vapply(candidates, function(column) {
    values <- unique(trimws(as.character(metadata[[column]])))
    values <- values[!is.na(values) & nzchar(values)]
    length(values) >= 2L && length(values) <= max(50L, floor(nrow(metadata) / 2L))
  }, logical(1))]
}

facet_text_scale <- function(columns_per_row) {
  columns_per_row <- suppressWarnings(as.integer(columns_per_row)[1L])
  if (is.na(columns_per_row) || columns_per_row <= 1L) return(1)
  if (columns_per_row == 2L) return(0.9)
  if (columns_per_row == 3L) return(0.8)
  0.7
}

cv_sample_groups <- function(metadata, group_col, samples) {
  metadata <- unique_sample_metadata(metadata)
  samples <- as.character(samples)
  if (!group_col %in% colnames(metadata)) return(stats::setNames(rep(NA_character_, length(samples)), samples))
  values <- as.character(metadata[[group_col]][match(samples, metadata$Sample)])
  stats::setNames(values, samples)
}

normalize_proteomics_text <- function(value) {
  value <- as.character(value)
  for (iteration in seq_len(3L)) {
    normalized <- gsub("&#x0*26;", "&", value, ignore.case = TRUE)
    normalized <- gsub("&#0*38;", "&", normalized, ignore.case = TRUE)
    normalized <- gsub("&amp;", "&", normalized, ignore.case = TRUE, fixed = FALSE)
    if (identical(normalized, value)) break
    value <- normalized
  }
  value
}

normalize_proteomics_metadata <- function(metadata) {
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
  text_columns <- vapply(metadata, function(column) is.character(column) || is.factor(column), logical(1))
  metadata[text_columns] <- lapply(metadata[text_columns], normalize_proteomics_text)
  metadata
}

stats_comparison_component <- function(value) gsub("[^A-Za-z0-9]+", "_", trimws(normalize_proteomics_text(value)))

make_stats_comparison_id <- function(group_col, numerator, denominator) {
  paste(normalize_proteomics_text(c(group_col, numerator, denominator)), collapse = "|||")
}

parse_stats_comparison_id <- function(comparison, default_group_col = "Condition") {
  parts <- normalize_proteomics_text(strsplit(as.character(comparison)[1L], "|||", fixed = TRUE)[[1L]])
  if (length(parts) == 2L) {
    return(list(group_col = default_group_col, numerator = parts[[1L]], denominator = parts[[2L]], legacy = TRUE))
  }
  if (length(parts) != 3L) stop("Invalid statistics comparison identifier.", call. = FALSE)
  list(group_col = parts[[1L]], numerator = parts[[2L]], denominator = parts[[3L]], legacy = FALSE)
}

stats_comparison_label <- function(comparison) {
  parsed <- parse_stats_comparison_id(comparison)
  if (isTRUE(parsed$legacy)) paste(parsed$numerator, "vs", parsed$denominator) else paste0(parsed$group_col, ": ", parsed$numerator, " vs ", parsed$denominator)
}

stats_comparison_prefix <- function(comparison) {
  parsed <- parse_stats_comparison_id(comparison)
  pieces <- c(if (!isTRUE(parsed$legacy)) parsed$group_col, parsed$numerator, "vs", parsed$denominator)
  paste(stats_comparison_component(pieces), collapse = "_")
}
