
library(shiny)
library(ggplot2)
library(DT)
library(dplyr)
library(stringr)
library(missMDA)
library(FactoMineR)
library(svglite)

options(shiny.maxRequestSize = 1024 * 1024 * 1024)

project_duckdb_file <- c(file.path(getwd(), "project_duckdb.R"), file.path(dirname(normalizePath("app.R", mustWork = FALSE)), "project_duckdb.R"))
project_duckdb_file <- project_duckdb_file[file.exists(project_duckdb_file)][1]
if (!is.na(project_duckdb_file)) source(project_duckdb_file)

metadata_state_file <- c(
  file.path(getwd(), "metadata_state.R"),
  file.path(dirname(normalizePath("app.R", mustWork = FALSE)), "metadata_state.R")
)
metadata_state_file <- metadata_state_file[file.exists(metadata_state_file)][1]
if (!is.na(metadata_state_file)) source(metadata_state_file)
project_rds_file <- c(file.path(getwd(), "project_rds.R"), file.path(dirname(normalizePath("app.R", mustWork = FALSE)), "project_rds.R"))
project_rds_file <- project_rds_file[file.exists(project_rds_file)][1]
if (!is.na(project_rds_file)) source(project_rds_file)
selection_helpers_file <- c(file.path(getwd(), "selection_helpers.R"), file.path(dirname(normalizePath("app.R", mustWork = FALSE)), "selection_helpers.R"))
selection_helpers_file <- selection_helpers_file[file.exists(selection_helpers_file)][1]
if (!is.na(selection_helpers_file)) source(selection_helpers_file)
imputation_helpers_file <- c(file.path(getwd(), "imputation_helpers.R"), file.path(dirname(normalizePath("app.R", mustWork = FALSE)), "imputation_helpers.R"))
imputation_helpers_file <- imputation_helpers_file[file.exists(imputation_helpers_file)][1]
if (!is.na(imputation_helpers_file)) source(imputation_helpers_file)

condition_replicate_label <- function(condition, replicate, fallback) {
  condition <- trimws(as.character(condition))
  replicate <- trimws(as.character(replicate))
  fallback <- as.character(fallback)
  complete <- !is.na(condition) & nzchar(condition) & !is.na(replicate) & nzchar(replicate)
  out <- fallback
  out[complete] <- paste(condition[complete], replicate[complete], sep = "_")
  out
}

parse_sample_exclusion_text <- function(text) {
  if (is.null(text) || !nzchar(trimws(as.character(text)[1]))) return(character(0))
  values <- unlist(strsplit(as.character(text)[1], "[\r\n]+"))
  values <- trimws(values)
  values <- values[!is.na(values) & nzchar(values)]
  unique(values)
}

normalize_pca_missing_values <- function(matrix_data) {
  matrix_data <- as.matrix(matrix_data)
  matrix_data[!is.finite(matrix_data)] <- NA_real_
  matrix_data
}

subset_pca_samples <- function(expression, metadata, metadata_column = "", included_values = character(0)) {
  expression <- as.matrix(expression)
  sample_ids <- rownames(expression)
  if (is.null(sample_ids)) stop("PCA expression matrix requires sample row names.", call. = FALSE)
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"Sample" %in% colnames(metadata)) stop("PCA metadata requires a Sample column.", call. = FALSE)
  metadata <- metadata[match(sample_ids, as.character(metadata$Sample)), , drop = FALSE]
  metadata$Sample <- sample_ids

  included_values <- trimws(as.character(included_values))
  included_values <- unique(included_values[!is.na(included_values) & nzchar(included_values)])
  use_filter <- length(metadata_column) == 1L && nzchar(metadata_column) &&
    metadata_column %in% colnames(metadata) && length(included_values) > 0L
  keep <- if (use_filter) {
    values <- trimws(as.character(metadata[[metadata_column]]))
    !is.na(values) & values %in% included_values
  } else {
    rep(TRUE, length(sample_ids))
  }

  list(
    expression = expression[keep, , drop = FALSE],
    metadata = metadata[keep, , drop = FALSE],
    excluded_samples = as.integer(sum(!keep))
  )
}

project_input_table <- function(file_info, cached_table = NULL) {
  if (is.null(file_info)) return(cached_table)
  value <- function(name) {
    item <- file_info[[name]]
    if (is.null(item) || !length(item)) "" else as.character(item[[1L]])
  }
  datapath <- value("datapath")
  filename <- value("name")
  if (!nzchar(datapath) || !file.exists(datapath)) stop("The selected table file is not available.", call. = FALSE)
  ext <- tolower(tools::file_ext(filename))
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  read.table(
    datapath,
    header = TRUE,
    sep = sep,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = ""
  )
}

prepare_metadata_editor_view <- function(metadata, selected_columns, locked_columns = character(0)) {
  metadata <- normalize_proteomics_metadata(metadata)
  if (!"Sample" %in% colnames(metadata)) stop("Metadata must contain a Sample column.", call. = FALSE)
  selected_columns <- as.character(selected_columns)
  selected_columns <- selected_columns[selected_columns %in% colnames(metadata)]
  if (!length(selected_columns)) stop("Select at least one column for Table S1. Metadata.", call. = FALSE)
  if ("RunOrder" %in% colnames(metadata)) {
    metadata <- metadata[order(is.na(metadata$RunOrder), metadata$RunOrder, seq_len(nrow(metadata))), , drop = FALSE]
  }
  list(
    display = metadata[, selected_columns, drop = FALSE],
    sample_keys = as.character(metadata$Sample),
    locked_indices = which(selected_columns %in% locked_columns) - 1L
  )
}

sample_exclusion_flags <- function(md, exclusion_terms) {
  if (is.null(md) || nrow(md) == 0) return(logical(0))
  exclusion_terms <- parse_sample_exclusion_text(paste(exclusion_terms, collapse = "\n"))
  if (length(exclusion_terms) == 0) return(rep(FALSE, nrow(md)))
  match_cols <- intersect(
    c("Sample", "SampleName", "AnalysisLabel", "Run Label", "File Name", "Filename", "RunLabel"),
    colnames(md)
  )
  if (length(match_cols) == 0) return(rep(FALSE, nrow(md)))
  row_values <- lapply(seq_len(nrow(md)), function(row_index) {
    values <- unlist(md[row_index, match_cols, drop = FALSE], use.names = FALSE)
    values <- as.character(values)
    values <- values[!is.na(values) & nzchar(trimws(values))]
    stems <- tools::file_path_sans_ext(basename(values))
    unique(tolower(c(values, stems)))
  })
  terms <- tolower(exclusion_terms)
  vapply(row_values, function(values) {
    exact <- any(terms %in% values)
    substring_terms <- terms[nchar(terms) >= 6]
    contained <- length(substring_terms) > 0 && any(vapply(substring_terms, function(term) any(grepl(term, values, fixed = TRUE)), logical(1)))
    exact || contained
  }, logical(1))
}

filter_expression_columns_by_exclusion <- function(expr, exclusion_terms, md = NULL) {
  if (is.null(expr) || ncol(expr) <= 1) return(expr)
  sample_cols <- colnames(expr)[-1]
  if (!is.null(md) && nrow(md) > 0 && "Sample" %in% colnames(md)) {
    flags <- sample_exclusion_flags(md, exclusion_terms)
    excluded_samples <- as.character(md$Sample[flags])
  } else {
    sample_md <- data.frame(Sample = sample_cols, stringsAsFactors = FALSE)
    flags <- sample_exclusion_flags(sample_md, exclusion_terms)
    excluded_samples <- sample_cols[flags]
  }
  keep_cols <- c(TRUE, !sample_cols %in% excluded_samples)
  expr[, keep_cols, drop = FALSE]
}

parse_position_list <- function(text, max_position = 96) {
  if (is.null(text) || !nzchar(trimws(as.character(text)[1]))) return(integer(0))
  parts <- unlist(strsplit(as.character(text)[1], "[,;[:space:]]+"))
  values <- integer(0)
  for (part in parts) {
    part <- trimws(part)
    if (!nzchar(part)) next
    if (grepl("^[0-9]+[-:][0-9]+$", part)) {
      bounds <- as.integer(unlist(strsplit(part, "[-:]")))
      if (length(bounds) == 2 && all(is.finite(bounds))) {
        values <- c(values, seq(bounds[1], bounds[2]))
      }
    } else {
      values <- c(values, suppressWarnings(as.integer(part)))
    }
  }
  values <- values[is.finite(values) & values >= 1 & values <= max_position]
  unique(values)
}

evosep_plate_grid <- function() {
  positions <- seq_len(96)
  data.frame(
    Position = positions,
    Row = LETTERS[((positions - 1) %/% 12) + 1],
    Column = ((positions - 1) %% 12) + 1,
    Well = paste0(LETTERS[((positions - 1) %/% 12) + 1], ((positions - 1) %% 12) + 1),
    stringsAsFactors = FALSE
  )
}

evosep_xml_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

evosep_xml_unescape <- function(x) {
  x <- as.character(x)
  x <- gsub("&quot;", '"', x, fixed = TRUE)
  x <- gsub("&gt;", ">", x, fixed = TRUE)
  x <- gsub("&lt;", "<", x, fixed = TRUE)
  x <- gsub("&amp;", "&", x, fixed = TRUE)
  x
}

evosep_extract_grd_rows <- function(template_text) {
  worksheet_match <- regexpr(
    "(?s)<Worksheet[^>]+ss:Name=\"grdSampleList\"[^>]*>.*?</Worksheet>",
    template_text,
    perl = TRUE
  )
  if (worksheet_match[1] < 0) stop("No grdSampleList worksheet was found in the CSL template.")
  worksheet <- regmatches(template_text, worksheet_match)
  row_matches <- gregexpr("(?s)<Row[^>]*>.*?</Row>", worksheet, perl = TRUE)[[1]]
  if (row_matches[1] < 0) stop("The grdSampleList worksheet does not contain rows.")
  row_blocks <- regmatches(worksheet, list(row_matches))[[1]]
  lapply(row_blocks, function(row_block) {
    data_matches <- gregexpr("(?s)<Data[^>]*>(.*?)</Data>", row_block, perl = TRUE)[[1]]
    if (data_matches[1] < 0) return(character(0))
    data_blocks <- regmatches(row_block, list(data_matches))[[1]]
    values <- sub("^<Data[^>]*>", "", data_blocks)
    values <- sub("</Data>$", "", values)
    evosep_xml_unescape(values)
  })
}

evosep_template_defaults <- function(template_text) {
  rows <- evosep_extract_grd_rows(template_text)
  if (length(rows) < 2 || length(rows[[2]]) < 10) stop("The CSL template must contain at least one sample row with 10 columns.")
  id_rows <- rows[vapply(rows, function(row) length(row) >= 7 && grepl("^ID[0-9]+", row[7]), logical(1))]
  default_row <- if (length(id_rows) > 0) id_rows[[1]] else rows[[2]]
  id_filename <- default_row[7]
  sample_suffix <- sub("^ID[0-9]+", "", id_filename)
  if (identical(sample_suffix, id_filename)) sample_suffix <- ""
  list(
    analysis_method = default_row[2],
    source_tray = default_row[3],
    xcalibur_method = default_row[6],
    output_dir = default_row[9],
    comment = default_row[10],
    sample_filename_suffix = sample_suffix
  )
}

evosep_incremented_suffix <- function(suffix, index) {
  suffix <- as.character(suffix)
  if (!nzchar(suffix)) return(paste0("_", sprintf("%02d", index)))
  if (grepl("^_[0-9]{2}", suffix)) {
    sub("^_[0-9]{2}", paste0("_", sprintf("%02d", index)), suffix)
  } else {
    paste0("_", sprintf("%02d", index), suffix)
  }
}

evosep_normalize_assignment <- function(values) {
  values <- trimws(as.character(values))
  values[is.na(values)] <- ""
  lower <- tolower(values)
  out <- rep("Empty", length(values))
  out[lower %in% c("adh", "standard")] <- "ADH"
  out[lower %in% c("spqc", "qc", "pool", "pooled qc", "pooled_qc")] <- "SPQC"
  out[lower %in% c("sample", "manifest", "id")] <- "Sample"
  out
}

evosep_normalize_plate_table <- function(plate, manifest = NULL, suffix = "", spqc_prefix = "IDSPQC") {
  grid <- evosep_plate_grid()
  required <- c("Position", "Well", "Row", "Column", "Assignment", "ManifestID", "SampleName", "Filename")
  for (column in required) {
    if (!column %in% colnames(plate)) plate[[column]] <- ""
  }
  plate$Position <- suppressWarnings(as.integer(plate$Position))
  plate <- plate[!is.na(plate$Position) & plate$Position %in% grid$Position, , drop = FALSE]
  plate <- plate[!duplicated(plate$Position), , drop = FALSE]
  out <- merge(grid, plate[, setdiff(colnames(plate), c("Row", "Column", "Well")), drop = FALSE], by = "Position", all.x = TRUE, sort = TRUE)
  out <- out[order(out$Position), , drop = FALSE]
  for (column in c("Assignment", "ManifestID", "SampleName", "Filename")) {
    out[[column]][is.na(out[[column]])] <- ""
    out[[column]] <- trimws(as.character(out[[column]]))
  }
  out$Assignment <- evosep_normalize_assignment(out$Assignment)
  out$ManifestID[out$Assignment != "Sample"] <- ""

  manifest_ids <- character(0)
  manifest_names <- character(0)
  if (!is.null(manifest) && nrow(manifest) > 0 && all(c("SampleDetailsID", "SampleName") %in% colnames(manifest))) {
    manifest_ids <- as.character(manifest$SampleDetailsID)
    manifest_names <- as.character(manifest$SampleName)
  }
  name_by_id <- stats::setNames(manifest_names, manifest_ids)
  used_ids <- character(0)
  sample_rows <- which(out$Assignment == "Sample")
  for (row_index in sample_rows) {
    current_id <- trimws(out$ManifestID[row_index])
    if (!nzchar(current_id) || !current_id %in% manifest_ids || current_id %in% used_ids) {
      next_id <- setdiff(manifest_ids, used_ids)[1]
      current_id <- if (is.na(next_id)) "" else next_id
      out$ManifestID[row_index] <- current_id
    }
    if (nzchar(current_id)) used_ids <- c(used_ids, current_id)
    if (nzchar(current_id) && (!nzchar(out$SampleName[row_index]) || out$SampleName[row_index] == current_id)) {
      out$SampleName[row_index] <- unname(name_by_id[current_id])
    }
    if (nzchar(current_id) && !nzchar(out$Filename[row_index])) {
      out$Filename[row_index] <- paste0(current_id, suffix)
    }
  }

  adh_rows <- which(out$Assignment == "ADH")
  for (i in seq_along(adh_rows)) {
    row_index <- adh_rows[i]
    if (!nzchar(out$SampleName[row_index])) out$SampleName[row_index] <- "ADH"
    if (!nzchar(out$Filename[row_index])) out$Filename[row_index] <- paste0("ADH", evosep_incremented_suffix(suffix, i))
  }

  spqc_rows <- which(out$Assignment == "SPQC")
  if (!nzchar(spqc_prefix)) spqc_prefix <- "IDSPQC"
  for (i in seq_along(spqc_rows)) {
    row_index <- spqc_rows[i]
    if (!nzchar(out$SampleName[row_index])) out$SampleName[row_index] <- "SPQC"
    if (!nzchar(out$Filename[row_index])) out$Filename[row_index] <- paste0(spqc_prefix, evosep_incremented_suffix(suffix, i))
  }

  empty_rows <- which(out$Assignment == "Empty")
  out$SampleName[empty_rows] <- ""
  out$Filename[empty_rows] <- ""
  out[, required, drop = FALSE]
}

evosep_make_row_xml <- function(values) {
  cell_types <- c("Number", rep("String", 2), "Number", rep("String", 6))
  cells <- vapply(seq_along(values), function(i) {
    paste0(
      "    <Cell><Data ss:Type=\"",
      cell_types[i],
      "\">",
      evosep_xml_escape(values[i]),
      "</Data></Cell>"
    )
  }, character(1))
  paste(c("   <Row>", cells, "   </Row>"), collapse = "\n")
}

evosep_write_csl_text <- function(template_text, queue_data) {
  worksheet_match <- regexpr(
    "(?s)<Worksheet[^>]+ss:Name=\"grdSampleList\"[^>]*>.*?</Worksheet>",
    template_text,
    perl = TRUE
  )
  if (worksheet_match[1] < 0) stop("No grdSampleList worksheet was found in the CSL template.")
  worksheet <- regmatches(template_text, worksheet_match)
  row_matches <- gregexpr("(?s)<Row[^>]*>.*?</Row>", worksheet, perl = TRUE)[[1]]
  if (row_matches[1] < 0) stop("The grdSampleList worksheet does not contain rows.")
  row_blocks <- regmatches(worksheet, list(row_matches))[[1]]
  header_row <- row_blocks[1]
  data_rows <- vapply(seq_len(nrow(queue_data)), function(i) {
    evosep_make_row_xml(unlist(queue_data[i, c(
      "RunNumber", "AnalysisMethod", "SourceTray", "SourceVial", "SampleName",
      "XcaliburMethod", "XcaliburFilename", "PostAcquisitionProgram",
      "OutputDir", "Comment"
    )], use.names = FALSE))
  }, character(1))
  replacement_rows <- paste(c(header_row, data_rows), collapse = "\n")
  first_row_start <- row_matches[1]
  last_row_index <- length(row_matches)
  last_row_end <- row_matches[last_row_index] + attr(row_matches, "match.length")[last_row_index] - 1
  worksheet <- paste0(
    substr(worksheet, 1, first_row_start - 1),
    replacement_rows,
    substr(worksheet, last_row_end + 1, nchar(worksheet))
  )
  worksheet <- sub(
    "ss:ExpandedRowCount=\"[0-9]+\"",
    paste0("ss:ExpandedRowCount=\"", nrow(queue_data) + 1, "\""),
    worksheet,
    perl = TRUE
  )
  sub(regmatches(template_text, worksheet_match), worksheet, template_text, fixed = TRUE)
}

replicate_pair_plan <- function(md, header_labels, numerator, denominator, group_col = "Condition", pair_col = "Replicate") {
  if (is.null(pair_col) || !length(pair_col) || is.na(pair_col[1]) || !nzchar(pair_col[1])) pair_col <- "Replicate"
  if (!group_col %in% colnames(md)) return(list(balanced = FALSE, numerator_labels = character(0), denominator_labels = character(0), reason = paste0("metadata column not found: ", group_col)))
  if (!pair_col %in% colnames(md)) return(list(balanced = FALSE, numerator_labels = character(0), denominator_labels = character(0), reason = paste0("pairing metadata column not found: ", pair_col)))
  pairing <- data.frame(
    Group = trimws(as.character(md[[group_col]])),
    Replicate = trimws(as.character(md[[pair_col]])),
    HeaderLabel = as.character(header_labels),
    stringsAsFactors = FALSE
  )
  numerator_rows <- pairing[pairing$Group == numerator, , drop = FALSE]
  denominator_rows <- pairing[pairing$Group == denominator, , drop = FALSE]

  blank_replicates <- any(is.na(numerator_rows$Replicate) | !nzchar(numerator_rows$Replicate)) ||
    any(is.na(denominator_rows$Replicate) | !nzchar(denominator_rows$Replicate))
  duplicate_replicates <- anyDuplicated(numerator_rows$Replicate) > 0 ||
    anyDuplicated(denominator_rows$Replicate) > 0
  numerator_replicates <- sort(numerator_rows$Replicate)
  denominator_replicates <- sort(denominator_rows$Replicate)
  enough_replicates <- length(numerator_replicates) >= 2 && length(denominator_replicates) >= 2
  same_replicates <- identical(numerator_replicates, denominator_replicates)

  if (blank_replicates || duplicate_replicates || !enough_replicates || !same_replicates) {
    reason <- if (blank_replicates) {
      "blank replicate values"
    } else if (duplicate_replicates) {
      "duplicate replicate values"
    } else if (!enough_replicates) {
      "fewer than two replicates per group"
    } else {
      "replicate sets differ between groups"
    }
    return(list(
      balanced = FALSE,
      numerator_labels = character(0),
      denominator_labels = character(0),
      reason = paste0(pair_col, ": ", reason)
    ))
  }

  numerator_rows <- numerator_rows[match(numerator_replicates, numerator_rows$Replicate), , drop = FALSE]
  denominator_rows <- denominator_rows[match(numerator_replicates, denominator_rows$Replicate), , drop = FALSE]
  list(
    balanced = TRUE,
    numerator_labels = numerator_rows$HeaderLabel,
    denominator_labels = denominator_rows$HeaderLabel,
    reason = ""
  )
}

metadata_spqc_rows <- function(data) {
  if (is.null(data) || nrow(data) == 0) return(logical(0))
  search_cols <- intersect(c("Sample", "Run Label", "SampleName", "AnalysisLabel"), colnames(data))
  if (!length(search_cols)) return(rep(FALSE, nrow(data)))
  search_text <- do.call(paste, c(lapply(search_cols, function(col) as.character(data[[col]])), sep = " "))
  rows <- grepl("SPQC", search_text, ignore.case = TRUE)
  rows[is.na(rows)] <- FALSE
  rows
}

spqc_metadata_editable_columns <- function(data) {
  setdiff(colnames(data), "Sample")
}

calculate_protein_comparison <- function(report, numerator_cols, denominator_cols, paired = FALSE) {
  numeric_matrix <- function(cols) {
    values <- as.data.frame(report[, cols, drop = FALSE], stringsAsFactors = FALSE)
    values[] <- lapply(values, function(column) suppressWarnings(as.numeric(as.character(column))))
    as.matrix(values)
  }

  numerator_raw <- numeric_matrix(numerator_cols)
  denominator_raw <- numeric_matrix(denominator_cols)
  numerator_raw[!is.finite(numerator_raw) | numerator_raw <= 0] <- NA_real_
  denominator_raw[!is.finite(denominator_raw) | denominator_raw <= 0] <- NA_real_

  if (isTRUE(paired)) {
    valid_pairs <- is.finite(numerator_raw) & is.finite(denominator_raw)
    numerator_raw[!valid_pairs] <- NA_real_
    denominator_raw[!valid_pairs] <- NA_real_
    pair_counts <- rowSums(valid_pairs)
    log2_ratios <- log2(numerator_raw / denominator_raw)
    log2_fc <- rowMeans(log2_ratios, na.rm = TRUE)
    log2_fc[pair_counts < 2 | !is.finite(log2_fc)] <- NaN
    p_value <- vapply(seq_len(nrow(report)), function(i) {
      valid <- valid_pairs[i, ]
      if (sum(valid) < 2) return(NaN)
      test <- tryCatch(
        stats::t.test(log2(numerator_raw[i, valid]), log2(denominator_raw[i, valid]), paired = TRUE),
        error = function(e) NULL
      )
      if (is.null(test)) NaN else unname(test$p.value)
    }, numeric(1))
    return(list(data = data.frame(log2_fc, p_value, check.names = FALSE), method = "paired"))
  }

  numerator_log2 <- log2(numerator_raw)
  denominator_log2 <- log2(denominator_raw)
  log2_fc <- rowMeans(numerator_log2, na.rm = TRUE) - rowMeans(denominator_log2, na.rm = TRUE)
  log2_fc[!is.finite(log2_fc)] <- NaN
  p_value <- vapply(seq_len(nrow(report)), function(i) {
    numerator_values <- numerator_log2[i, ]
    denominator_values <- denominator_log2[i, ]
    numerator_values <- numerator_values[!is.na(numerator_values)]
    denominator_values <- denominator_values[!is.na(denominator_values)]
    if (length(numerator_values) < 2 || length(denominator_values) < 2) return(NaN)
    test <- tryCatch(
      stats::t.test(numerator_values, denominator_values, var.equal = FALSE),
      error = function(e) NULL
    )
    if (is.null(test)) NaN else unname(test$p.value)
  }, numeric(1))
  list(data = data.frame(log2_fc, p_value, check.names = FALSE), method = "unpaired")
}

protein_header_labels_from_metadata <- function(md, header_columns) {
  header_columns <- header_columns[header_columns %in% colnames(md)]
  if (length(header_columns) == 0) {
    header_columns <- if ("SampleName" %in% colnames(md)) {
      "SampleName"
    } else if ("AnalysisLabel" %in% colnames(md)) {
      "AnalysisLabel"
    } else {
      "Sample"
    }
  }
  header_parts <- lapply(header_columns, function(column) {
    values <- as.character(md[[column]])
    values[is.na(values)] <- ""
    if (identical(column, "SampleName") && "Replicate" %in% colnames(md)) {
      missing_sample_name <- !nzchar(values)
      replicate_values <- as.character(md$Replicate)
      replicate_values[is.na(replicate_values)] <- ""
      values[missing_sample_name & nzchar(replicate_values)] <- replicate_values[missing_sample_name & nzchar(replicate_values)]
    }
    values
  })
  header_labels <- vapply(seq_len(nrow(md)), function(row_index) {
    parts <- vapply(header_parts, `[`, character(1), row_index)
    parts <- trimws(parts)
    parts <- parts[!is.na(parts) & nzchar(parts)]
    if (length(parts) == 0) return(as.character(md$Sample[row_index]))
    paste(parts, collapse = "_")
  }, character(1))
  header_labels <- trimws(header_labels)
  missing_header_labels <- is.na(header_labels) | !nzchar(header_labels)
  header_labels[missing_header_labels] <- as.character(md$Sample[missing_header_labels])
  header_labels
}

preferred_metadata_column <- function(columns, preferred_names) {
  columns <- columns[!is.na(columns) & nzchar(columns)]
  if (length(columns) == 0) return(NA_character_)
  normalized_columns <- tolower(gsub("[^a-z0-9]+", "", columns))
  normalized_preferred <- tolower(gsub("[^a-z0-9]+", "", preferred_names))
  for (preferred in normalized_preferred) {
    matched <- which(normalized_columns == preferred)
    if (length(matched) > 0) return(columns[matched[1]])
  }
  NA_character_
}

fill_missing_metadata_values <- function(existing, replacement) {
  replacement <- as.character(replacement)
  replacement[is.na(replacement)] <- ""
  if (is.null(existing)) return(replacement)
  existing <- as.character(existing)
  trimmed_existing <- trimws(existing)
  missing_existing <- is.na(existing) | !nzchar(trimmed_existing) | tolower(trimmed_existing) %in% c("not defined", "na", "n/a")
  existing[missing_existing] <- replacement[missing_existing]
  existing
}

run_label_date_token <- function(values) {
  values <- as.character(values)
  values[is.na(values)] <- ""
  matches <- regmatches(values, gregexpr("(?<![0-9])[0-9]{6}(?![0-9])", values, perl = TRUE))
  vapply(matches, function(tokens) {
    if (length(tokens) == 0) return("")
    tokens[length(tokens)]
  }, character(1))
}

infer_spqc_batches_from_run_dates <- function(samples, batches, spqc_rows) {
  samples <- as.character(samples)
  batches <- as.character(batches)
  batches[is.na(batches)] <- ""
  spqc_rows[is.na(spqc_rows)] <- FALSE
  date_tokens <- run_label_date_token(samples)
  inferred <- batches
  missing_spqc_batch <- spqc_rows & (!nzchar(trimws(inferred)) | tolower(trimws(inferred)) %in% c("not defined", "na", "n/a"))
  if (!any(missing_spqc_batch)) return(inferred)

  for (date_token in unique(date_tokens[missing_spqc_batch])) {
    if (!nzchar(date_token)) next
    same_date_regular <- !spqc_rows & date_tokens == date_token & nzchar(trimws(batches)) &
      !tolower(trimws(batches)) %in% c("not defined", "na", "n/a")
    inherited_batch <- if (any(same_date_regular)) {
      batch_counts <- sort(table(batches[same_date_regular]), decreasing = TRUE)
      names(batch_counts)[1]
    } else {
      date_token
    }
    inferred[missing_spqc_batch & date_tokens == date_token] <- inherited_batch
  }
  inferred
}

parse_spqc_batch_overrides <- function(text) {
  if (is.null(text) || !nzchar(trimws(as.character(text)[1]))) {
    return(data.frame(Pattern = character(0), Batch = character(0), stringsAsFactors = FALSE))
  }
  lines <- unlist(strsplit(as.character(text)[1], "\r?\n"))
  rows <- lapply(lines, function(line) {
    line <- trimws(line)
    if (!nzchar(line) || startsWith(line, "#")) return(NULL)
    parts <- unlist(strsplit(line, "\\s*(=|,|\t)\\s*", perl = TRUE))
    if (length(parts) < 2) return(NULL)
    pattern <- trimws(parts[1])
    batch <- trimws(paste(parts[-1], collapse = " "))
    if (!nzchar(pattern) || !nzchar(batch)) return(NULL)
    data.frame(Pattern = pattern, Batch = batch, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    return(data.frame(Pattern = character(0), Batch = character(0), stringsAsFactors = FALSE))
  }
  dplyr::bind_rows(rows)
}

apply_spqc_batch_overrides <- function(samples, sample_names, batches, spqc_rows, override_text) {
  overrides <- parse_spqc_batch_overrides(override_text)
  if (nrow(overrides) == 0) return(batches)
  search_text <- paste(samples, sample_names)
  search_text[is.na(search_text)] <- ""
  out <- as.character(batches)
  for (i in seq_len(nrow(overrides))) {
    matched <- spqc_rows & grepl(overrides$Pattern[i], search_text, ignore.case = TRUE, fixed = TRUE)
    out[matched] <- overrides$Batch[i]
  }
  out
}

batch_correction_confounded <- function(md, batch_col, group_col) {
  if (is.null(batch_col) || is.null(group_col) || !batch_col %in% colnames(md) || !group_col %in% colnames(md)) {
    return(list(confounded = FALSE, message = "Batch or group column is not available."))
  }
  values <- data.frame(
    Batch = trimws(as.character(md[[batch_col]])),
    Group = trimws(as.character(md[[group_col]])),
    stringsAsFactors = FALSE
  )
  values <- values[!is.na(values$Batch) & nzchar(values$Batch) & !is.na(values$Group) & nzchar(values$Group), , drop = FALSE]
  if (nrow(values) == 0 || length(unique(values$Batch)) < 2 || length(unique(values$Group)) < 2) {
    return(list(confounded = FALSE, message = "Not enough non-empty batch/group values to assess confounding."))
  }
  groups_per_batch <- tapply(values$Group, values$Batch, function(x) length(unique(x)))
  batches_per_group <- tapply(values$Batch, values$Group, function(x) length(unique(x)))
  confounded <- all(groups_per_batch == 1) || all(batches_per_group == 1)
  message <- if (confounded) {
    "Batch and group appear confounded; batch correction may remove biology."
  } else {
    "Batch and group are not perfectly confounded."
  }
  list(confounded = confounded, message = message)
}

batch_correct_prepare_input <- function(report, md, batch_col, group_col, header_label_columns, feature_col_requested = NULL, source_label = "S3", pseudocount = 1e-6, min_batches_for_feature = 2) {
  stopifnot("Sample" %in% colnames(md))
  if (!batch_col %in% colnames(md)) stop("Batch metadata column not found: ", batch_col, call. = FALSE)
  if (!group_col %in% colnames(md)) stop("Group metadata column not found: ", group_col, call. = FALSE)
  if (anyDuplicated(md$Sample)) stop("Metadata Sample values must be unique for batch correction.", call. = FALSE)

  quantity_cols <- proteomics_abundance_columns(report)
  column_labels <- proteomics_abundance_sample_names(quantity_cols)
  restored_processed <- any(endsWith(quantity_cols, "_Protein_group_abundance"))
  metadata_keys <- if (restored_processed) protein_header_labels_from_metadata(md, header_label_columns) else as.character(md$Sample)
  keep <- column_labels %in% metadata_keys
  quantity_cols <- quantity_cols[keep]
  column_labels <- column_labels[keep]
  if (length(quantity_cols) < 3) stop("Need at least 3 matched S3 abundance columns for batch correction.", call. = FALSE)

  md2 <- md[match(column_labels, metadata_keys), , drop = FALSE]
  run_labels <- as.character(md2$Sample)
  valid_metadata <- !is.na(md2[[batch_col]]) & nzchar(trimws(as.character(md2[[batch_col]]))) &
    !is.na(md2[[group_col]]) & nzchar(trimws(as.character(md2[[group_col]])))
  quantity_cols <- quantity_cols[valid_metadata]
  run_labels <- run_labels[valid_metadata]
  md2 <- md2[valid_metadata, , drop = FALSE]
  if (length(quantity_cols) < 3) stop("Need at least 3 samples with non-empty batch and group metadata.", call. = FALSE)
  if (length(unique(as.character(md2[[batch_col]]))) < 2) stop("Need at least 2 batches for batch correction.", call. = FALSE)

  mat <- as.data.frame(report[, quantity_cols, drop = FALSE], stringsAsFactors = FALSE)
  mat[] <- lapply(mat, function(column) suppressWarnings(as.numeric(as.character(column))))
  mat <- as.matrix(mat)
  mat[!is.finite(mat) | mat <= 0] <- NA_real_
  min_batches_for_feature <- suppressWarnings(as.integer(min_batches_for_feature))
  if (!is.finite(min_batches_for_feature) || is.na(min_batches_for_feature) || min_batches_for_feature < 1) min_batches_for_feature <- 1L
  batches <- as.character(md2[[batch_col]])
  batch_levels <- unique(batches[!is.na(batches) & nzchar(trimws(batches))])
  min_batches_for_feature <- min(min_batches_for_feature, length(batch_levels))
  batch_presence <- integer(nrow(mat))
  for (batch in batch_levels) {
    batch_cols <- which(batches == batch)
    if (length(batch_cols) > 0) {
      batch_presence <- batch_presence + as.integer(rowSums(!is.na(mat[, batch_cols, drop = FALSE])) > 0)
    }
  }
  keep_feature_rows <- batch_presence >= min_batches_for_feature
  mat <- mat[keep_feature_rows, , drop = FALSE]
  if (nrow(mat) < 2) stop("Batch correction needs at least 2 protein groups after filtering by batch presence.", call. = FALSE)
  data_as_input <- as.data.frame(log2(mat + pseudocount), check.names = FALSE)
  colnames(data_as_input) <- run_labels
  rownames(data_as_input) <- paste0("F", seq_len(nrow(data_as_input)))

  header_labels <- protein_header_labels_from_metadata(md2, header_label_columns)
  if (anyDuplicated(header_labels)) stop("Batch-corrected sample labels are not unique. Add metadata columns to the header label setting.", call. = FALSE)

  feature_col <- resolve_report_feature_col(report, feature_col_requested)
  feature_values <- if (!is.na(feature_col) && feature_col %in% colnames(report)) protein_feature_labels(report, feature_col) else paste0("Feature_", seq_len(nrow(report)))
  info_cols <- intersect(c("PG.Genes", "PG.ProteinGroups", "PG.ProteinNames", "PG.ProteinDescriptions"), colnames(report))
  feature_info <- report[keep_feature_rows, info_cols, drop = FALSE]
  feature_info$Feature <- make.unique(feature_values[keep_feature_rows])

  sample_map <- data.frame(
    Sample = run_labels,
    HeaderLabel = header_labels,
    Batch = as.character(md2[[batch_col]]),
    Group = as.character(md2[[group_col]]),
    stringsAsFactors = FALSE
  )
  description_as_input <- data.frame(
    ID = run_labels,
    sample = seq_along(run_labels),
    batch = as.character(md2[[batch_col]]),
    group = as.character(md2[[group_col]]),
    stringsAsFactors = FALSE
  )
  list(
    data_as_input = data_as_input,
    description_as_input = description_as_input,
    sample_map = sample_map,
    feature_info = feature_info,
    kept_row_indices = which(keep_feature_rows),
    source_label = source_label,
    confounding = batch_correction_confounded(md2, batch_col, group_col),
    filter_summary = list(
      source_rows = nrow(report),
      kept_rows = nrow(mat),
      removed_rows = nrow(report) - nrow(mat),
      min_batches_for_feature = min_batches_for_feature,
      batch_presence_counts = batch_presence[keep_feature_rows]
    )
  )
}

batch_correct_rebuild_table <- function(report, corrected_values, prepared, suffix = "_batch_corrected_log2") {
  corrected_values <- as.data.frame(corrected_values, check.names = FALSE)
  corrected_values[] <- lapply(corrected_values, function(column) suppressWarnings(as.numeric(as.character(column))))
  corrected_values <- corrected_values[, prepared$sample_map$Sample, drop = FALSE]
  colnames(corrected_values) <- paste0(prepared$sample_map$HeaderLabel, suffix)
  info_cols <- intersect(c("PG.Genes", "PG.ProteinGroups", "PG.ProteinNames", "PG.ProteinDescriptions"), colnames(report))
  report_rows <- if (!is.null(prepared$kept_row_indices)) prepared$kept_row_indices else seq_len(nrow(corrected_values))
  out <- data.frame(report[report_rows, info_cols, drop = FALSE], corrected_values, check.names = FALSE)
  attr(out, "sample_map") <- prepared$sample_map
  attr(out, "confounding") <- prepared$confounding
  out
}

calculate_protein_comparison_log2 <- function(report, numerator_cols, denominator_cols, paired = FALSE) {
  numeric_matrix <- function(cols) {
    values <- as.data.frame(report[, cols, drop = FALSE], stringsAsFactors = FALSE)
    values[] <- lapply(values, function(column) suppressWarnings(as.numeric(as.character(column))))
    as.matrix(values)
  }
  numerator_log2 <- numeric_matrix(numerator_cols)
  denominator_log2 <- numeric_matrix(denominator_cols)

  if (isTRUE(paired)) {
    valid_pairs <- is.finite(numerator_log2) & is.finite(denominator_log2)
    numerator_log2[!valid_pairs] <- NA_real_
    denominator_log2[!valid_pairs] <- NA_real_
    pair_counts <- rowSums(valid_pairs)
    log2_fc <- rowMeans(numerator_log2 - denominator_log2, na.rm = TRUE)
    log2_fc[pair_counts < 2 | !is.finite(log2_fc)] <- NaN
    p_value <- vapply(seq_len(nrow(report)), function(i) {
      valid <- valid_pairs[i, ]
      if (sum(valid) < 2) return(NaN)
      test <- tryCatch(
        stats::t.test(numerator_log2[i, valid], denominator_log2[i, valid], paired = TRUE),
        error = function(e) NULL
      )
      if (is.null(test)) NaN else unname(test$p.value)
    }, numeric(1))
    return(list(data = data.frame(log2_fc, p_value, check.names = FALSE), method = "paired"))
  }

  log2_fc <- rowMeans(numerator_log2, na.rm = TRUE) - rowMeans(denominator_log2, na.rm = TRUE)
  log2_fc[!is.finite(log2_fc)] <- NaN
  p_value <- vapply(seq_len(nrow(report)), function(i) {
    numerator_values <- numerator_log2[i, ]
    denominator_values <- denominator_log2[i, ]
    numerator_values <- numerator_values[is.finite(numerator_values)]
    denominator_values <- denominator_values[is.finite(denominator_values)]
    if (length(numerator_values) < 2 || length(denominator_values) < 2) return(NaN)
    test <- tryCatch(
      stats::t.test(numerator_values, denominator_values, var.equal = FALSE),
      error = function(e) NULL
    )
    if (is.null(test)) NaN else unname(test$p.value)
  }, numeric(1))
  list(data = data.frame(log2_fc, p_value, check.names = FALSE), method = "unpaired")
}

append_log2_stats_to_protein_table <- function(report, md, sample_map, comparisons, paired_comparisons = character(0), include_fdr = TRUE, abundance_suffix = "_batch_corrected_log2", pair_col = "Replicate") {
  if (is.null(comparisons)) comparisons <- character(0)
  if (!length(comparisons)) {
    attr(report, "stats_comparison") <- character(0)
    attr(report, "stats_methods") <- character(0)
    return(report)
  }
  paired_comparisons <- if (is.null(paired_comparisons)) character(0) else paired_comparisons
  md2 <- md[match(sample_map$Sample, as.character(md$Sample)), , drop = FALSE]
  md2$HeaderLabel <- sample_map$HeaderLabel

  stats_blocks <- lapply(comparisons, function(comparison) {
    parsed_comparison <- parse_stats_comparison_id(comparison)
    group_col <- parsed_comparison$group_col
    numerator <- parsed_comparison$numerator
    denominator <- parsed_comparison$denominator
    if (!group_col %in% colnames(md2) || numerator == denominator) stop("Each statistics comparison must use a valid metadata column and different groups.", call. = FALSE)
    requested_paired <- comparison %in% paired_comparisons
    pairing <- replicate_pair_plan(md2, md2$HeaderLabel, numerator, denominator, group_col = group_col, pair_col = pair_col)
    effective_paired <- requested_paired && pairing$balanced
    fallback_reason <- ""

    if (effective_paired) {
      numerator_cols <- paste0(pairing$numerator_labels, abundance_suffix)
      denominator_cols <- paste0(pairing$denominator_labels, abundance_suffix)
      if (!all(numerator_cols %in% colnames(report)) || !all(denominator_cols %in% colnames(report))) {
        effective_paired <- FALSE
        fallback_reason <- "matched abundance columns are missing"
      }
    }
    if (!effective_paired) {
      numerator_cols <- paste0(as.character(md2$HeaderLabel[as.character(md2[[group_col]]) == numerator]), abundance_suffix)
      denominator_cols <- paste0(as.character(md2$HeaderLabel[as.character(md2[[group_col]]) == denominator]), abundance_suffix)
      if (requested_paired && !pairing$balanced) fallback_reason <- pairing$reason
    }
    numerator_cols <- numerator_cols[numerator_cols %in% colnames(report)]
    denominator_cols <- denominator_cols[denominator_cols %in% colnames(report)]
    if (length(numerator_cols) < 2 || length(denominator_cols) < 2) {
      stop("Comparison ", numerator, " vs ", denominator, " requires at least two matched abundance columns per group.", call. = FALSE)
    }

    result <- calculate_protein_comparison_log2(report, numerator_cols, denominator_cols, paired = effective_paired)
    comparison_prefix <- stats_comparison_prefix(comparison)
    comparison_df <- data.frame(result$data$log2_fc, result$data$p_value, check.names = FALSE)
    colnames(comparison_df) <- c(
      paste0(comparison_prefix, "_log2_fold_change"),
      paste0(comparison_prefix, "_", result$method, "_t_test_p_value")
    )
    if (isTRUE(include_fdr)) {
      comparison_df[[paste0(comparison_prefix, "_BH_FDR")]] <- stats::p.adjust(result$data$p_value, method = "BH")
    }
    method_note <- paste0(stats_comparison_label(comparison), ": ", result$method, " on batch-corrected log2 values")
    if (requested_paired && !effective_paired) {
      method_note <- paste0(method_note, " (automatic fallback: ", fallback_reason, ")")
    }
    attr(comparison_df, "method_note") <- method_note
    comparison_df
  })

  stats_df <- do.call(cbind, stats_blocks)
  out <- data.frame(report, stats_df, check.names = FALSE)
  attr(out, "stats_comparison") <- vapply(comparisons, stats_comparison_label, character(1))
  attr(out, "stats_methods") <- vapply(stats_blocks, function(block) attr(block, "method_note"), character(1))
  out
}

comparison_p_value_column <- function(report, prefix) {
  candidates <- paste0(prefix, c("_paired_t_test_p_value", "_unpaired_t_test_p_value"))
  candidates <- candidates[candidates %in% colnames(report)]
  if (length(candidates) == 1) candidates else NA_character_
}

infer_project_number <- function(file_names) {
  file_names <- basename(as.character(file_names))
  file_names <- file_names[!is.na(file_names) & nzchar(file_names)]
  if (!length(file_names)) return(NA_character_)
  matches <- stringr::str_match(file_names, "(?:^|[^0-9])([0-9]{5})(?:[^0-9]|$)")[, 2]
  matches <- matches[!is.na(matches) & nzchar(matches)]
  if (length(matches)) matches[1] else NA_character_
}

project_bundle_filename <- function(file_names, date = Sys.Date(), project_number = NULL) {
  project_number <- clean_project_number(project_number)
  if (is.na(project_number)) project_number <- infer_project_number(file_names)
  date_text <- format(as.Date(date), "%m%d%y")
  if (is.na(project_number)) {
    paste0("proteomics_data_workup_project_", date_text, ".zip")
  } else {
    paste0("proteomics_data_workup_project_", project_number, "_", date_text, ".zip")
  }
}

file_integrity <- function(path) {
  info <- file.info(path)
  if (is.na(info$size) || !file.exists(path)) stop("Cannot inspect project file: ", path)
  list(size_bytes = unname(as.numeric(info$size)), md5 = unname(tools::md5sum(path)))
}

safe_archive_entry <- function(name) {
  if (length(name) != 1 || is.na(name) || !nzchar(name)) return(FALSE)
  normalized <- gsub("\\\\", "/", name)
  !grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))", normalized)
}

clean_project_number <- function(project_number) {
  project_number <- trimws(as.character(project_number)[1])
  if (is.na(project_number) || !nzchar(project_number)) return(NA_character_)
  project_number <- gsub("[^A-Za-z0-9_-]", "_", project_number)
  project_number <- gsub("_+", "_", project_number)
  project_number <- gsub("^_|_$", "", project_number)
  if (nzchar(project_number)) project_number else NA_character_
}

add_project_number_to_filename <- function(filename, project_number) {
  project_number <- clean_project_number(project_number)
  if (is.na(project_number)) return(filename)
  ext <- tools::file_ext(filename)
  stem <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", filename) else filename
  if (startsWith(stem, paste0(project_number, "_"))) return(filename)
  if (nzchar(ext)) paste0(project_number, "_", stem, ".", ext) else paste0(project_number, "_", stem)
}

pca_loadings_matrix <- function(pca_model) {
  if (!is.null(pca_model$rotation)) {
    return(as.matrix(pca_model$rotation))
  }
  if (!is.null(pca_model$var) && !is.null(pca_model$var$coord)) {
    return(as.matrix(pca_model$var$coord))
  }
  matrix(numeric(0), nrow = 0, ncol = 0)
}

rank_pca_loadings <- function(pca_model, protein_info = NULL, pc_mode = "combined", top_n = 50) {
  loadings <- pca_loadings_matrix(pca_model)
  if (nrow(loadings) == 0 || ncol(loadings) < 2) {
    return(data.frame())
  }
  pc1 <- suppressWarnings(as.numeric(loadings[, 1]))
  pc2 <- suppressWarnings(as.numeric(loadings[, 2]))
  proteins <- rownames(loadings)
  if (is.null(proteins) || any(!nzchar(proteins))) proteins <- paste0("Protein_", seq_len(nrow(loadings)))
  out <- data.frame(
    Protein = proteins,
    PCAFeatureID = proteins,
    PC1_Loading = pc1,
    PC2_Loading = pc2,
    Combined_PC1_PC2 = sqrt(pc1^2 + pc2^2),
    stringsAsFactors = FALSE
  )
  if (!is.null(protein_info) && nrow(protein_info) > 0 && "Protein" %in% colnames(protein_info)) {
    out <- merge(out, protein_info, by = "Protein", all.x = TRUE, sort = FALSE)
    annotation_cols <- intersect(c("PG.Genes", "PG.ProteinNames", "PG.ProteinGroups", "PG.ProteinDescriptions"), colnames(out))
    needs_raw_join <- if (length(annotation_cols) == 0) {
      rep(TRUE, nrow(out))
    } else {
      apply(out[, annotation_cols, drop = FALSE], 1, function(row) all(is.na(row) | !nzchar(trimws(as.character(row)))))
    }
    if (any(needs_raw_join)) {
      join_keys <- c("RawFeatureID", "SourceRowIndex", "KeptRowIndex")
      unresolved <- which(needs_raw_join)
      for (join_key in join_keys[join_keys %in% colnames(protein_info)]) {
        if (length(unresolved) == 0) break
        key_values <- trimws(as.character(protein_info[[join_key]]))
        raw_match <- match(trimws(as.character(out$PCAFeatureID[unresolved])), key_values)
        matched <- unresolved[!is.na(raw_match)]
        raw_match <- raw_match[!is.na(raw_match)]
        if (length(matched) > 0) {
          for (column_name in setdiff(colnames(protein_info), "Protein")) {
            if (!column_name %in% colnames(out)) out[[column_name]] <- NA_character_
            out[[column_name]][matched] <- protein_info[[column_name]][raw_match]
          }
          out$Protein[matched] <- protein_info$Protein[raw_match]
          unresolved <- setdiff(unresolved, matched)
        }
      }
    }
  }
  out$Protein <- best_protein_display_label(out, fallback_col = "Protein")
  rank_value <- switch(
    pc_mode,
    "PC1 positive" = out$PC1_Loading,
    "PC1 negative" = -out$PC1_Loading,
    "PC2 positive" = out$PC2_Loading,
    "PC2 negative" = -out$PC2_Loading,
    out$Combined_PC1_PC2
  )
  out <- out[order(rank_value, decreasing = TRUE, na.last = NA), , drop = FALSE]
  top_n <- suppressWarnings(as.integer(top_n))
  if (is.finite(top_n) && top_n > 0) out <- head(out, top_n)
  rownames(out) <- NULL
  out
}

split_gene_symbols <- function(values) {
  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(values)]
  genes <- unlist(strsplit(values, "[;|,[:space:]]+"), use.names = FALSE)
  genes <- trimws(genes)
  unique(genes[!is.na(genes) & nzchar(genes)])
}

protein_feature_labels <- function(df, requested_col = NULL) {
  protein_identifier_cols <- c(
    "PG.Genes", "PG.ProteinNames", "PG.ProteinGroups",
    "Protein.Group", "ProteinGroups", "Genes", "Gene"
  )
  preferred_cols <- c(
    requested_col,
    "PG.Genes",
    "PG.ProteinNames",
    "PG.ProteinGroups",
    "Protein.Group",
    "ProteinGroups",
    "Genes",
    "Gene",
    colnames(df)[1]
  )
  preferred_cols <- unique(preferred_cols[!is.na(preferred_cols) & nzchar(preferred_cols) & preferred_cols %in% colnames(df)])
  if (length(preferred_cols) == 0) return(as.character(seq_len(nrow(df))))

  out <- rep(NA_character_, nrow(df))
  for (column_name in preferred_cols) {
    values <- trimws(as.character(df[[column_name]]))
    numeric_only <- grepl("^[0-9]+$", values)
    if (!column_name %in% protein_identifier_cols) {
      values[numeric_only] <- NA_character_
    }
    usable <- !is.na(values) & nzchar(values) & (is.na(out) | !nzchar(out))
    out[usable] <- values[usable]
  }
  out[is.na(out) | !nzchar(out)] <- as.character(seq_len(sum(is.na(out) | !nzchar(out))))
  out
}

resolve_report_feature_col <- function(df, requested_col = NULL) {
  requested_col <- if (is.null(requested_col) || length(requested_col) == 0) "" else as.character(requested_col)[1]
  if (nzchar(requested_col) && requested_col %in% colnames(df)) return(requested_col)
  candidate_cols <- c("PG.Genes", "PG.ProteinGroups", "PG.ProteinNames", colnames(df)[1])
  candidate <- candidate_cols[candidate_cols %in% colnames(df)][1]
  if (is.na(candidate)) NA_character_ else candidate
}

best_protein_display_label <- function(df, fallback_col = "Protein") {
  preferred_cols <- c("PG.Genes", "PG.ProteinNames", "PG.ProteinGroups", "Protein.Group", "ProteinGroups", fallback_col)
  preferred_cols <- unique(preferred_cols[preferred_cols %in% colnames(df)])
  out <- rep(NA_character_, nrow(df))
  for (column_name in preferred_cols) {
    values <- trimws(as.character(df[[column_name]]))
    usable <- !is.na(values) & nzchar(values) & (is.na(out) | !nzchar(out))
    out[usable] <- values[usable]
  }
  out[is.na(out) | !nzchar(out)] <- as.character(seq_len(sum(is.na(out) | !nzchar(out))))
  out
}

ranked_gene_statistics <- function(proteins, genes, log2_fc, p_value, fdr = NULL, metric = "p_value") {
  genes_list <- lapply(genes, split_gene_symbols)
  significance <- if (identical(metric, "BH_FDR") && !is.null(fdr)) fdr else p_value
  rank_score <- sign(log2_fc) * -log10(pmax(significance, .Machine$double.xmin))
  rows <- lapply(seq_along(genes_list), function(i) {
    if (!is.finite(rank_score[i]) || length(genes_list[[i]]) == 0) return(NULL)
    data.frame(
      Gene = genes_list[[i]],
      Protein = proteins[i],
      RankScore = rank_score[i],
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(stats::setNames(numeric(0), character(0)))
  out <- do.call(rbind, rows)
  if (nrow(out) == 0) return(stats::setNames(numeric(0), character(0)))
  out <- out[is.finite(out$RankScore) & !is.na(out$Gene) & out$Gene != "", , drop = FALSE]
  if (nrow(out) == 0) return(stats::setNames(numeric(0), character(0)))
  out <- out[order(out$Gene, -abs(out$RankScore)), , drop = FALSE]
  out <- out[!duplicated(out$Gene), , drop = FALSE]
  out <- out[order(out$RankScore, decreasing = TRUE), , drop = FALSE]
  ranks <- out$RankScore
  names(ranks) <- out$Gene
  ranks
}

ui <- fluidPage(
  titlePanel("Proteomics Data Workup"),
  tabsetPanel(id = "workflow_tabs",
        tabPanel(
          "Make metadata",
          h4("Open or save project"),
          fluidRow(
            column(width = 6, fileInput("project_open_file", "Open existing project", accept = c(".duckdb", ".db", ".rds"))),
            column(width = 3, actionButton("create_new_project", "Create new project"))
          ),
          textOutput("project_path_preview"),
          fluidRow(
            column(width = 6, downloadButton("download_project_duckdb", "Save As / Download project")),
            column(width = 3, actionButton("clear_active_project", "Clear project", class = "btn-warning"))
          ),
          tags$small("The active project is a temporary session copy. Download it before closing the app to keep your changes."),
          verbatimTextOutput("project_status"),
          tags$small("DuckDB is the default project format. Existing DuckDB and RDS projects can be opened."),
          textInput("download_project_number", "Project number for downloads", value = ""),
          checkboxInput("include_project_number_in_downloads", "Add project number to downloaded filenames", FALSE),
          tags$small("Leave blank to infer the project number from loaded source filenames when possible."),
          checkboxInput("save_download_copy", "Also save downloads to a folder on this computer", FALSE),
          textInput("download_destination_dir", "Download copy destination folder", value = ""),
          verbatimTextOutput("download_destination_note"),
          verbatimTextOutput("project_bundle_note"),
          h4("Active project data"),
          tags$small("Only data currently loaded or restored in this project are shown below."),
          DTOutput("project_files_table"),
          br(),
          downloadButton("download_active_metadata", "Download active metadata CSV"),
          fileInput("metadata_replacement_file", "Upload modified metadata CSV", accept = c(".csv", ".tsv", ".txt")),
          actionButton("apply_metadata_replacement", "Load modified metadata into draft"),
          verbatimTextOutput("metadata_replacement_note"),
          tags$hr(),
          h4("Global sample exclusions"),
          p("Enter one outlier sample per line. Matching samples are removed from analysis views and exports that use sample measurements, while the metadata keeps a record of the exclusion."),
          textAreaInput(
            "sample_exclusions_text",
            "Samples/run labels to exclude",
            value = "",
            rows = 4,
            placeholder = "One per line, for example:\nID12345_01_10751_STAR1_102125.htrms"
          ),
          textInput("sample_exclusion_reason", "Exclusion reason", value = "PCA outlier"),
          verbatimTextOutput("sample_exclusion_note"),
          DTOutput("sample_exclusion_preview"),
          tags$hr(),
          h4("Build analysis metadata"),
          actionButton("apply_metadata_changes", "Apply metadata changes", class = "btn-primary"),
          textOutput("metadata_apply_status"),
          p("Use the condition setup as the metadata basis, the run-order table for run order, and the sample-details workbook to add submitted sample names plus any other workbook columns such as batch and group."),
          fileInput("meta_file", "Select condition setup table", accept = c(".csv", ".tsv", ".txt")),
          fileInput("run_order_file", "Select run-order table", accept = c(".csv", ".tsv", ".txt")),
          fileInput("sample_details_file", "Select order sample details workbook", accept = c(".xlsx", ".xls")),
          tags$small("Choose each file from its location on your computer. The three sources are combined into the analysis metadata table below; sample-details batch/group columns are available for export, PCA, stats, and batch correction."),
          verbatimTextOutput("metadata_loaded_files_note"),
          tags$hr(),
          h4("SPQC assignment"),
          selectInput(
            "spqc_assignment_mode",
            "Assign SPQC samples to",
            choices = c(
              "One shared SPQC group" = "single",
              "Separate SPQC group by inferred batch" = "batch",
              "Separate SPQC group by run-label date" = "date",
              "Keep existing SPQC condition/group values" = "keep"
            ),
            selected = "single"
          ),
          textInput("spqc_group_label", "Shared SPQC group label", value = "SPQC"),
          textInput("spqc_group_prefix", "SPQC subgroup prefix", value = "SPQC"),
          textAreaInput(
            "spqc_batch_overrides",
            "Optional SPQC batch overrides",
            value = "",
            rows = 3,
            placeholder = "One rule per line, for example:\nSPQC_02_OA10222_10751_STAR1_102125 = 1"
          ),
          tags$small("SPQC samples are detected from Run Label/SampleName containing SPQC. Batch can still be inferred from the run-label date even when SPQC samples are not in the manifest."),
          tags$hr(),
          verbatimTextOutput("metadata_build_note"),
          h4("Supplementary table export"),
          selectizeInput(
            "metadata_export_columns",
            "Columns in Table S1. Metadata (drag to reorder)",
            choices = NULL,
            selected = NULL,
            multiple = TRUE,
            options = proteomics_multiselect_options("Select metadata columns for export", drag = TRUE)
          ),
          actionButton("reset_metadata_columns", "Restore default columns"),
          tags$small("The Excel worksheet follows this column order. Rows are ordered by RunOrder when a run-order file is selected."),
          textInput("metadata_workbook_filename", "Excel workbook filename", value = "supplementary_tables.xlsx"),
          h4("Editable metadata"),
          tags$small("Edit any metadata field except Sample, run order, and filename identifiers. Changes remain drafts until Apply metadata changes is clicked."),
          actionButton("discard_metadata_edits", "Discard all metadata edits"),
          DTOutput("metadata_preview"),
          tags$hr(),
          downloadButton("download_built_metadata", "Download Table S1 CSV"),
          downloadButton("download_metadata_workbook", "Download Excel workbook")
        ),
        tabPanel(
          "Protein tables",
          h4("Add protein supplementary tables"),
          p("Start with the non-imputed protein group report. Table S3 can use a Spectronaut-imputed report, app-generated kNN imputation, or Table S2 without imputation."),
          fileInput("protein_no_impute_file", "Select protein report (no imputation)", accept = c(".csv", ".tsv", ".txt")),
          h4("Missing-value imputation for Table S3"),
          radioButtons(
            "s3_imputation_source",
            "Table S3 source",
            choices = c(
              "Upload a Spectronaut-imputed protein report" = "spectronaut",
              "Generate Table S3 with kNN from Table S2" = "knn",
              "Use Table S2 without additional imputation" = "s2"
            ),
            selected = "spectronaut"
          ),
          conditionalPanel(
            condition = "input.s3_imputation_source == 'spectronaut'",
            fileInput("protein_imputed_file", "Select Spectronaut-imputed protein report", accept = c(".csv", ".tsv", ".txt"))
          ),
          conditionalPanel(
            condition = "input.s3_imputation_source == 'knn'",
            numericInput("protein_knn_k", "Number of nearest proteins (k)", value = 10, min = 1, step = 1),
            numericInput("protein_knn_max_missing_percent", "Exclude proteins missing more than (%)", value = 50, min = 0, max = 100, step = 5),
            radioButtons(
              "protein_knn_scope",
              "kNN scope",
              choices = c("All samples together" = "global", "Separately within a metadata field" = "metadata"),
              selected = "global",
              inline = TRUE
            ),
            conditionalPanel(
              condition = "input.protein_knn_scope == 'metadata'",
              selectInput("protein_knn_group_col", "Metadata field defining imputation groups", choices = NULL)
            ),
            actionButton("run_protein_knn", "Generate Table S3 with kNN", class = "btn-primary")
          ),
          verbatimTextOutput("protein_imputation_status"),
          tags$small("App kNN operates on log2 protein-group abundances, changes only missing values, and converts results back to abundance scale. The selected source is exported as 'Table S3. Protein, imputed'."),
          verbatimTextOutput("protein_loaded_files_note"),
          selectizeInput(
            "protein_header_label_columns",
            "Rename sample measurement headers using metadata variables (drag to reorder)",
            choices = NULL,
            selected = c("SampleName"),
            multiple = TRUE,
            options = proteomics_multiselect_options("Select variables such as Condition, SampleName, or Replicate", drag = TRUE)
          ),
          tags$small("Selected values are joined with underscores, for example Condition_SampleName or Condition_Replicate. Use a unique combination so replicate columns do not collide."),
          selectizeInput(
            "protein_quantity_order_columns",
            "Order protein quantity columns by metadata (drag to reorder)",
            choices = NULL,
            selected = NULL,
            multiple = TRUE,
            options = proteomics_multiselect_options("For example: Condition, Batch, Replicate, RunOrder", drag = TRUE)
          ),
          tags$small("Leave blank to keep the source report order. When selected, sample measurement columns are sorted by these metadata fields before export."),
          h4("Protein information columns"),
          selectizeInput(
            "s2_non_data_columns",
            "Table S2 non-data columns (drag to reorder)",
            choices = NULL,
            selected = NULL,
            multiple = TRUE,
            options = proteomics_multiselect_options("Select Table S2 information columns", drag = TRUE)
          ),
          selectizeInput(
            "s3_non_data_columns",
            "Table S3 non-data columns (drag to reorder)",
            choices = NULL,
            selected = NULL,
            multiple = TRUE,
            options = proteomics_multiselect_options("Select Table S3 information columns", drag = TRUE)
          ),
          tags$small("Sample measurement columns remain in the protein tables after the selected non-data columns."),
          selectInput("protein_cv_group_col", "Metadata column for %CV calculation", choices = "Condition", selected = "Condition"),
          selectizeInput(
            "cv_conditions",
            "Groups for %CV calculation",
            choices = NULL,
            selected = "SPQC",
            multiple = TRUE,
            options = proteomics_multiselect_options("Select groups for %CV columns")
          ),
          tags$small("%CV columns use protein group abundance values and are inserted after the protein annotation columns."),
          radioButtons(
            "protein_derived_column_order",
            "Derived column order after protein information columns",
            choices = c("CV columns before statistics" = "cv_before_stats", "Statistics before CV columns" = "stats_before_cv"),
            selected = "cv_before_stats",
            inline = TRUE
          ),
          tags$hr(),
          h4("Protein-group statistics"),
          selectizeInput(
            "stats_group_columns",
            "Metadata columns for comparisons",
            choices = NULL,
            selected = "Condition",
            multiple = TRUE,
            options = proteomics_multiselect_options("Select one or more metadata columns", drag = TRUE)
          ),
          checkboxGroupInput(
            "stats_tables",
            "Add statistics to tables",
            choices = c(
              "Table S2. Protein, no impute" = "S2",
              "Table S3. Protein, imputed" = "S3",
              "Table S3. Protein, imputed, batch-corrected" = "S3_batch_corrected"
            ),
            selected = c("S2", "S3")
          ),
          selectizeInput(
            "stats_comparisons",
            "Comparisons (numerator vs denominator)",
              choices = NULL,
              selected = NULL,
              multiple = TRUE,
              options = stats_comparison_selectize_options()
            ),
            selectizeInput(
              "stats_paired_comparisons",
              "Paired comparisons",
              choices = NULL,
              selected = NULL,
              multiple = TRUE,
              options = proteomics_multiselect_options("Select comparisons to analyze as paired")
            ),
            checkboxInput("stats_bh_fdr", "Include Benjamini-Hochberg FDR", TRUE),
            selectInput("stats_pair_col", "Pair samples by metadata", choices = "Replicate", selected = "Replicate"),
            tags$small("Drag selected comparisons to set their order. Paired comparisons match identical IDs in the selected pairing column, with one sample per ID in each group and at least two pairs. They use the mean paired log2 abundance ratio and a paired two-sided t-test on log2 abundances. Invalid or unbalanced pairing automatically falls back to the unpaired method."),
          actionButton("run_protein_stats", "Recalculate statistics", class = "btn-primary"),
          actionButton("stop_protein_stats", "Stop statistics run"),
          verbatimTextOutput("protein_stats_status"),
          h4("Derived-number export format"),
          numericInput("numeric_sig_figs", "Significant figures", value = 4, min = 1, max = 12, step = 1),
          checkboxInput("scientific_small_values", "Use scientific notation for small derived values", TRUE),
          numericInput("scientific_threshold", "Scientific notation when absolute value is below", value = 0.001, min = 0, step = 0.0001),
          tags$small("Formatting applies to %CV and statistical output columns in Excel. Protein group abundance columns are always displayed at two significant digits."),
          checkboxInput("enable_protein_filters", "Enable Excel filters on S2 and S3 protein tables", FALSE),
          tags$small("Disabled by default because Excel may become unstable when filtering very wide protein-data worksheets."),
          tags$hr(),
          verbatimTextOutput("protein_table_note"),
          h4("Table S2 preview"),
          downloadButton("download_protein_s2_csv", "Download Table S2 CSV"),
          uiOutput("protein_no_impute_preview_ui"),
          h4("Table S3 preview"),
          downloadButton("download_protein_s3_csv", "Download Table S3 CSV"),
          uiOutput("protein_imputed_preview_ui")
        ),
        tabPanel(
          "Batch correction",
          h4("S3 batch correction with ComBat"),
          p("Use the imputed protein table as the default input. Batch is assigned from metadata, and the selected biological group is preserved during correction."),
          fluidRow(
            column(
              width = 4,
              radioButtons(
                "batch_correction_source",
                "Input table",
                choices = c("Table S3. Protein, imputed" = "S3"),
                selected = "S3"
              ),
              selectInput("batch_correction_batch_col", "Batch metadata column", choices = c("Batch")),
              selectInput("batch_correction_group_col", "Biological group to preserve", choices = c("Condition"), selected = "Condition"),
              numericInput("batch_correction_combat_mode", "ComBat mode", value = 1, min = 1, max = 2, step = 1),
              numericInput("batch_correction_pseudocount", "Log2 pseudocount", value = 1e-6, min = 0, step = 1e-6),
              numericInput("batch_correction_min_batches_for_feature", "Keep proteins present in at least this many batches", value = 2, min = 1, step = 1),
              tags$small("Uses direct sva::ComBat when available to avoid HarmonizR parallel-backend errors. Values <= 0 become missing, and a protein is kept if at least one value is present in N batches."),
              actionButton("run_batch_correction", "Run batch correction"),
              downloadButton("download_batch_corrected_s3_csv", "Download corrected S3 CSV")
            ),
            column(
              width = 8,
              verbatimTextOutput("batch_cache_note"),
              verbatimTextOutput("batch_correction_note"),
              h4("Corrected S3 preview"),
              DTOutput("batch_corrected_s3_preview")
            )
          )
        ),
        tabPanel(
          "CV plots",
          fluidRow(
            column(
              width = 4,
              h4("Protein Group-CV Distribution per Condition"),
              radioButtons(
                "cv_plot_source",
                "CV plot data source",
                choices = c(
                  "Upload existing CV distribution table" = "uploaded",
                  "Calculate from non-imputed protein report" = "no_impute",
                  "Calculate from imputed protein report" = "imputed",
                  "Calculate from batch-corrected imputed protein report" = "S3_batch_corrected"
                ),
                selected = "uploaded"
              ),
              conditionalPanel(
                condition = "input.cv_plot_source == 'uploaded'",
                fileInput("cv_distribution_file", "Select CV distribution table", accept = c(".csv", ".tsv", ".txt"))
              ),
              conditionalPanel(
                condition = "input.cv_plot_source == 'no_impute' || input.cv_plot_source == 'imputed' || input.cv_plot_source == 'S3_batch_corrected'",
                selectInput("cv_plot_group_col", "Group CV by metadata field", choices = character(0), selected = "Condition"),
                sliderInput("cv_density_adjust", "Density smoothness", min = 0.1, max = 2, value = 1, step = 0.1),
                tags$small("1.0 uses the default smoothing; lower values produce a less-smoothed curve with more local detail.")
              ),
              selectizeInput("cv_plot_conditions", "Groups to plot", choices = NULL, selected = NULL, multiple = TRUE, options = proteomics_multiselect_options()),
              numericInput("cv_x_cutoff", "Maximum %CV on x-axis", value = 350, min = 1, step = 10),
              textInput("cv_plot_title", "CV plot title", value = "Protein Group-CV Distribution per Condition"),
              checkboxInput("cv_fill_density", "Fill density areas", FALSE),
              h4("CV plot text and line sizes"),
              numericInput("cv_title_size", "Title size", value = 12, min = 6, max = 36, step = 1),
              numericInput("cv_axis_title_size", "Axis title size", value = 11, min = 6, max = 30, step = 1),
              numericInput("cv_axis_text_size", "Axis tick label size", value = 9, min = 5, max = 24, step = 1),
              numericInput("cv_legend_text_size", "Legend text size", value = 9, min = 5, max = 24, step = 1),
              numericInput("cv_median_text_size", "Median label size", value = 3, min = 1, max = 10, step = 0.5),
              numericInput("cv_line_width", "Density line width", value = 0.65, min = 0.1, max = 3, step = 0.05),
              h4("Figure size"),
              numericInput("cv_figure_width", "Export width (inches)", value = 16, min = 2, max = 40, step = 0.5),
              numericInput("cv_figure_height", "Export height (inches)", value = 5, min = 2, max = 30, step = 0.5),
              downloadButton("download_cv_png", "Download CV plot PNG"),
              downloadButton("download_cv_svg", "Download CV plot SVG")
            ),
            column(
              width = 8,
              div(
                id = "cv_plot_panel",
                uiOutput("cv_plot_ui"),
                verbatimTextOutput("cv_plot_note")
              )
            )
          )
        ),
        tabPanel(
          "PCA",
          fluidRow(
            column(
              width = 4,
              h4("ClustVis-like PCA"),
              actionButton("run_clustvis_pca", "Run PCA"),
              tags$hr(),
              radioButtons(
                "clustvis_pca_source",
                "Abundance source",
                choices = c(
                  "Imputed protein report (Table S3 upload)" = "imputed",
                  "Non-imputed protein report (Table S2 upload)" = "no_impute",
                  "Batch-corrected imputed protein report (Table S3)" = "S3_batch_corrected"
                ),
                selected = "imputed"
              ),
              tags$hr(),
              h4("Sample subset"),
              selectInput(
                "clustvis_pca_subset_column",
                "Subset samples by metadata",
                choices = c("All samples" = ""),
                selected = ""
              ),
              selectizeInput(
                "clustvis_pca_subset_values",
                "Values to include",
                choices = NULL,
                selected = NULL,
                multiple = TRUE,
                options = proteomics_multiselect_options("Leave empty to include all values")
              ),
              numericInput("clustvis_pca_min_observed_percent", "Minimum observed samples per protein (%)", value = 70, min = 0, max = 100, step = 5),
              numericInput("clustvis_pca_npcs", "SVD imputation / PCA rank", value = 5, min = 2, max = 20, step = 1),
              checkboxInput("clustvis_pca_scale", "Scale features in prcomp", TRUE),
              checkboxInput("clustvis_pca_ellipses", "Show condition confidence ellipses", TRUE),
              tags$hr(),
              h4("Metadata plotting"),
              selectInput("clustvis_pca_color_by", "Color by metadata", choices = c("Condition"), selected = "Condition"),
              selectInput("clustvis_pca_shape_by", "Shape by metadata", choices = c("Replicate"), selected = "Replicate"),
              selectInput("clustvis_pca_label_by", "Label by metadata", choices = c("None"), selected = "None"),
              numericInput("clustvis_pca_point_size", "Point size", value = 3, min = 0.1, max = 10, step = 0.5),
              numericInput("clustvis_pca_label_size", "Label size", value = 3, min = 0.1, max = 10, step = 0.5),
              checkboxInput("clustvis_pca_no_fill_shapes", "Use no-fill / hollow point shapes", FALSE),
              checkboxInput("clustvis_pca_interactive", "Show interactive PCA with mouseover", FALSE),
              tags$hr(),
              h4("Condition opacity"),
              selectizeInput(
                "clustvis_pca_opacity_override_groups",
                "Apply opacity override to selected color groups",
                choices = NULL,
                selected = NULL,
                multiple = TRUE,
                options = proteomics_multiselect_options("Select one or more groups")
              ),
              sliderInput("clustvis_pca_override_opacity", "Selected-group opacity", min = 0, max = 1, value = 0.95, step = 0.05),
              sliderInput("clustvis_pca_default_opacity", "Other-group opacity", min = 0, max = 1, value = 0.65, step = 0.05),
              tags$small("Groups come from the current 'Color by metadata' field. Colors stay as the standard PCA colors; this only changes opacity."),
              tags$hr(),
              textInput("clustvis_pca_title", "Plot title", value = "ClustVis-like PCA"),
              numericInput("clustvis_pca_width", "Figure width (inches)", value = 10, min = 2, max = 30, step = 0.5),
              numericInput("clustvis_pca_height", "Figure height (inches)", value = 8, min = 2, max = 30, step = 0.5),
              downloadButton("download_clustvis_pca_png", "Download PCA PNG"),
              downloadButton("download_clustvis_pca_svg", "Download PCA SVG")
            ),
            column(
              width = 8,
              uiOutput("clustvis_pca_plot_ui"),
              verbatimTextOutput("clustvis_pca_note"),
              h4("PCA scores"),
              DTOutput("clustvis_pca_scores_table"),
              tags$hr(),
              h4("PCA loadings / protein contributors"),
              tags$small("Proteins have loadings on PC1 and PC2. Larger absolute loadings contribute more to the separation shown in the PCA plot."),
              fluidRow(
                column(
                  width = 5,
                  selectInput(
                    "pca_loading_rank_by",
                    "Rank contributors by",
                    choices = c(
                      "Combined PC1 and PC2" = "combined",
                      "PC1 positive" = "PC1 positive",
                      "PC1 negative" = "PC1 negative",
                      "PC2 positive" = "PC2 positive",
                      "PC2 negative" = "PC2 negative"
                    ),
                    selected = "combined"
                  )
                ),
                column(width = 3, numericInput("pca_loading_top_n", "Top proteins", value = 50, min = 1, max = 1000, step = 5)),
                column(width = 4, downloadButton("download_pca_loadings_csv", "Download PCA loadings CSV"))
              ),
              plotOutput("pca_loadings_plot", height = "320px"),
              DTOutput("pca_loadings_table")
            )
          )
        ),
        tabPanel(
          "Volcano plots",
          fluidRow(
            column(
              width = 4,
              h4("Protein-group volcano plot"),
              radioButtons(
                "volcano_source",
                "Statistics source",
                choices = c(
                  "Table S3. Protein, imputed" = "S3",
                  "Table S2. Protein, no impute" = "S2",
                  "Table S3. Protein, imputed, batch-corrected" = "S3_batch_corrected"
                ),
                selected = "S3"
              ),
              selectInput("volcano_comparison", "Comparison", choices = NULL),
              selectInput(
                "volcano_significance_metric",
                "Significance metric",
                choices = c("Benjamini-Hochberg FDR" = "BH_FDR", "Unadjusted p-value" = "p_value"),
                selected = "BH_FDR"
              ),
              numericInput("volcano_fc_cutoff", "Absolute log2 fold-change cutoff", value = 1, min = 0, step = 0.1),
              numericInput("volcano_sig_cutoff", "Significance cutoff", value = 0.05, min = 0.000000000001, max = 1, step = 0.01),
              selectInput("volcano_label_col", "Label proteins by", choices = NULL),
              radioButtons(
                "volcano_label_mode",
                "Labels to show",
                choices = c(
                  "Highlighted table rows" = "selected",
                  "Top significant proteins" = "automatic",
                  "No labels" = "none"
                ),
                selected = "selected"
              ),
              conditionalPanel(
                condition = "input.volcano_label_mode == 'automatic'",
                numericInput("volcano_max_labels", "Proteins to label per direction", value = 5, min = 0, step = 1)
              ),
              textInput("volcano_title", "Plot title", value = "Protein Group Volcano Plot"),
              numericInput("volcano_point_size", "Point size", value = 2, min = 0.1, step = 0.2),
              numericInput("volcano_label_size", "Label size", value = 3, min = 1, step = 0.5),
              numericInput("volcano_figure_width", "Figure width (inches)", value = 8, min = 2, max = 30, step = 0.5),
              numericInput("volcano_figure_height", "Figure height (inches)", value = 7, min = 2, max = 30, step = 0.5),
              checkboxInput("volcano_interactive", "Show interactive volcano plot", FALSE),
              downloadButton("download_volcano_png", "Download volcano PNG"),
              downloadButton("download_volcano_svg", "Download volcano SVG"),
              downloadButton("download_volcano_html", "Download interactive plot ZIP")
            ),
            column(
              width = 8,
              uiOutput("volcano_plot_ui"),
              verbatimTextOutput("volcano_note"),
              h4("Significant protein counts by comparison"),
              tags$small("Select a row to display that comparison and test metric in the volcano plot. Counts use the fold-change and significance cutoff above and report separate rows for BH FDR and unadjusted p-value. Increased means higher in the numerator condition; decreased means lower in the numerator condition."),
              DTOutput("volcano_counts_table"),
              downloadButton("download_volcano_counts_csv", "Download significant counts CSV"),
              h4("Significant proteins"),
              tags$small("Highlight one or more rows to label those proteins on the volcano plot when 'Highlighted table rows' is selected."),
              DTOutput("volcano_hits_table")
            )
          )
        ),
        tabPanel(
          "Feature plot",
          fluidRow(
            column(
              width = 4,
              h4("Feature abundance plot options"),
              radioButtons(
                "feature_data_source",
                "Feature plot abundance source",
                choices = c(
                  "Imputed protein report (Table S3)" = "imputed",
                  "Non-imputed protein report (Table S2)" = "no_impute",
                  "Batch-corrected imputed protein report (Table S3)" = "S3_batch_corrected"
                ),
                selected = "imputed"
              ),
              selectizeInput(
                "feature_select",
                "Protein / feature(s)",
                choices = NULL,
                selected = NULL,
                multiple = TRUE,
                options = proteomics_multiselect_options("Type to search; select one or more proteins", drag = TRUE, create = FALSE, maxOptions = 5000)
              ),
              numericInput(
                "feature_plot_ncol",
                "Feature plots per row (when multiple selected)",
                value = 3,
                min = 1,
                max = 12,
                step = 1
              ),
              selectizeInput(
                "feature_order_columns",
                "Order samples by metadata (drag to reorder)",
                choices = NULL,
                selected = NULL,
                multiple = TRUE,
                options = proteomics_multiselect_options("For example: Group, Batch, RunOrder", drag = TRUE)
              ),
              selectInput("feature_group_by", "Color/group samples by metadata", choices = c("Condition"), selected = "Condition"),
              selectInput("feature_label_by", "Label samples by metadata", choices = c("AnalysisLabel"), selected = "AnalysisLabel"),
              textInput("feature_plot_title", "Feature plot title", value = ""),
              numericInput("feature_title_size", "Feature title size", value = 18, min = 8, max = 40, step = 1),
              numericInput("feature_bar_width", "Bar width", value = 0.7, min = 0.1, max = 1.0, step = 0.05),
              numericInput("feature_text_size", "Sample label size", value = 3, min = 0.1, step = 0.5),
              radioButtons(
                "feature_value_scale",
                "Feature plot value scale",
                choices = c("Raw abundance" = "raw", "Log2 abundance" = "log2"),
                selected = "raw"
              ),
              selectInput("feature_color_mode", "Bar color values by", choices = c("Z-score", "Centered value", "Raw value"), selected = "Z-score"),
              checkboxInput("feature_symmetric_scale", "Force symmetric color scale around zero", TRUE),
              selectInput("feature_group_style", "Differentiate groups by", choices = c("Outline color", "Colored x labels", "None"), selected = "Outline color"),
              checkboxInput("show_feature_mean", "Show mean line", TRUE),
              checkboxInput("rotate_feature_labels", "Rotate sample labels", TRUE),
              checkboxInput("feature_interactive", "Show interactive/zoomable feature plot", FALSE),
              numericInput("feature_figure_width", "Figure width (inches)", value = 12, min = 2, max = 40, step = 0.5),
              numericInput("feature_figure_height", "Figure height (inches)", value = 6, min = 2, max = 30, step = 0.5)
            ),
            column(
              width = 8,
              uiOutput("feature_plot_ui"),
              verbatimTextOutput("feature_plot_note"),
              tags$hr(),
              downloadButton("download_feature_png", "Download feature PNG"),
              downloadButton("download_feature_svg", "Download feature SVG")
            )
          ),
          tags$hr(),
          h4("Log2 abundance matrix"),
            tags$small("Cells show log2 protein-group abundance and are colored by each protein's abundance z-score across samples (blue = lower, red = higher). Sort any column by clicking its header; select a row to open that protein in the Boxplot tab. The Significant column follows the active Volcano plots comparison and thresholds."),
          DTOutput("feature_matrix_table")
        ),
        tabPanel(
          "Boxplot",
          fluidRow(
            column(
              width = 4,
              h4("Protein boxplot"),
              tags$small("This uses the original script-style protein boxplot. Select a protein row in the Feature plot matrix or choose one below."),
              radioButtons(
                "script_box_source",
                "Abundance source",
                choices = c(
                  "Imputed protein report (Table S3)" = "imputed",
                  "Non-imputed protein report (Table S2)" = "no_impute",
                  "Batch-corrected imputed protein report (Table S3)" = "S3_batch_corrected"
                ),
                selected = "imputed"
              ),
              selectizeInput(
                "script_box_features",
                "Proteins / features",
                choices = NULL,
                selected = NULL,
                multiple = TRUE,
                options = proteomics_multiselect_options("Select one or more proteins", drag = TRUE, create = FALSE, maxOptions = 5000)
              ),
              numericInput(
                "script_box_ncol",
                "Boxplot panels per row (when multiple selected)",
                value = 3,
                min = 1,
                max = 12,
                step = 1
              ),
              selectInput("script_box_group_by", "Group samples by metadata", choices = c("Condition"), selected = "Condition"),
              selectInput("script_box_label_by", "Label points by metadata", choices = c("None"), selected = "None"),
              selectizeInput(
                "script_box_conditions",
                "Conditions/groups to include and order",
                choices = NULL,
                selected = NULL,
                multiple = TRUE,
                options = proteomics_multiselect_options("Leave blank to use all groups", drag = TRUE)
              ),
              radioButtons(
                "script_box_value_scale",
                "Abundance scale",
                choices = c("Raw abundance" = "raw", "Log2 abundance" = "log2"),
                selected = "log2",
                inline = TRUE
              ),
              radioButtons(
                "script_box_plot_style",
                "Plot style",
                choices = c("Boxplot + points" = "boxplot", "Mean +/- SD + points" = "mean_sd"),
                selected = "boxplot",
                inline = TRUE
              ),
              textInput("script_box_title", "Plot title (blank = feature)", value = ""),
              textInput("script_box_y_axis_title", "Y-axis title (blank = automatic)", value = ""),
              numericInput("script_box_point_size", "Point size", value = 3, min = 0.5, max = 10, step = 0.5),
              sliderInput("script_box_point_opacity", "Point opacity", min = 0.05, max = 1, value = 0.85, step = 0.05),
              numericInput("script_box_label_size", "Point label size", value = 3, min = 1, max = 10, step = 0.5),
              numericInput("script_box_text_size", "Text size", value = 18, min = 6, max = 36, step = 1),
              numericInput("script_box_width", "Figure width (inches)", value = 6.5, min = 2, max = 30, step = 0.5),
              numericInput("script_box_height", "Figure height (inches)", value = 5.5, min = 2, max = 30, step = 0.5),
              downloadButton("download_script_box_png", "Download boxplot PNG"),
              downloadButton("download_script_box_svg", "Download boxplot SVG")
            ),
            column(
              width = 8,
              uiOutput("script_box_plot_ui"),
              verbatimTextOutput("script_box_note")
            )
          )
        ),
        tabPanel(
          "Correlation",
          fluidRow(
            column(
              width = 4,
              h4("Protein correlation/regression"),
              radioButtons(
                "correlation_source",
                "Abundance source",
                choices = c(
                  "Imputed protein report (Table S3)" = "imputed",
                  "Non-imputed protein report (Table S2)" = "no_impute",
                  "Batch-corrected imputed protein report (Table S3)" = "S3_batch_corrected"
                ),
                selected = "imputed"
              ),
              selectizeInput(
                "correlation_feature",
                "Reference protein / feature",
                choices = NULL,
                selected = NULL,
                multiple = FALSE,
                options = list(placeholder = "Type to search protein / feature", create = FALSE, maxOptions = 5000)
              ),
              radioButtons(
                "correlation_value_scale",
                "Value scale for regression",
                choices = c("Log2 abundance" = "log2", "Raw abundance" = "raw"),
                selected = "log2",
                inline = TRUE
              ),
              radioButtons(
                "correlation_method",
                "Analysis method",
                choices = c(
                  "Linear regression" = "linear",
                  "Pearson correlation" = "pearson",
                  "Spearman correlation" = "spearman"
                ),
                selected = "linear"
              ),
              selectInput("correlation_group_by", "Within-group metadata column", choices = c("All samples" = "None"), selected = "None"),
              selectizeInput(
                "correlation_groups",
                "Within group(s)",
                choices = NULL,
                selected = NULL,
                multiple = TRUE,
                options = proteomics_multiselect_options("Leave blank to use all groups")
              ),
              selectizeInput(
                "correlation_covariates",
                "Adjust linear regression for metadata covariates",
                choices = NULL,
                selected = NULL,
                multiple = TRUE,
                options = proteomics_multiselect_options("Optional: Batch, Group, age, etc.")
              ),
              selectInput(
                "correlation_rank_by",
                "Rank top proteins by",
                choices = c(
                  "Smallest q-value" = "q_value",
                  "Largest R2" = "r_squared",
                  "Largest absolute beta" = "abs_beta",
                  "Smallest p-value" = "p_value"
                ),
                selected = "q_value"
              ),
              numericInput("correlation_top_n", "Top proteins to plot", value = 25, min = 1, max = 500, step = 1),
              checkboxInput("correlation_exclude_reference", "Exclude reference protein from results", TRUE),
              textInput("correlation_plot_title", "Plot title", value = "Top Protein Correlations"),
              numericInput("correlation_figure_width", "Figure width (inches)", value = 8, min = 2, max = 30, step = 0.5),
              numericInput("correlation_figure_height", "Figure height (inches)", value = 7, min = 2, max = 30, step = 0.5),
              downloadButton("download_correlation_png", "Download lollipop PNG"),
              downloadButton("download_correlation_svg", "Download lollipop SVG"),
              downloadButton("download_correlation_csv", "Download correlation table CSV")
            ),
            column(
              width = 8,
              uiOutput("correlation_lollipop_plot_ui"),
              verbatimTextOutput("correlation_note"),
              h4("All protein regressions"),
              tags$small("Each row compares a target protein against the selected reference protein across matched samples. Linear regression can optionally adjust for metadata covariates. QValue is the Benjamini-Hochberg adjusted p-value."),
              DTOutput("correlation_results_table")
            )
          )
        ),
        tabPanel(
          "Run identifications",
          fluidRow(
            column(
                width = 4,
                h4("Identification inputs"),
                fileInput("identifications_overview_file", "Select IdentificationsOverview table", accept = c(".csv", ".tsv", ".txt")),
                fileInput("run_identifications_precursor_file", "Select stacked precursor identifications table", accept = c(".csv", ".tsv", ".txt")),
                fileInput("run_identifications_protein_file", "Select stacked protein-group identifications table", accept = c(".csv", ".tsv", ".txt")),
                selectInput(
                  "identification_metric",
                  "Plot metric",
                  choices = c("Precursors" = "Precursors", "Protein groups" = "ProteinGroups"),
                  selected = "Precursors"
                ),
              textInput("identification_overview_title", "Overview plot title", value = "Identifications per Run"),
              textInput("run_identifications_title", "Stacked plot title", value = "Run Identifications"),
              numericInput("identification_title_size", "Plot title size", value = 13, min = 6, max = 36, step = 1),
              numericInput("identification_axis_text_size", "Axis label size", value = 9, min = 5, max = 24, step = 1),
              numericInput("identification_legend_size", "Legend text size", value = 9, min = 5, max = 24, step = 1),
              numericInput("identification_overview_width", "Overview width (inches)", value = 14, min = 2, max = 40, step = 0.5),
              numericInput("identification_overview_height", "Overview height (inches)", value = 5, min = 2, max = 30, step = 0.5),
              numericInput("run_identifications_width", "Stacked width (inches)", value = 16, min = 2, max = 40, step = 0.5),
              numericInput("run_identifications_height", "Stacked height (inches)", value = 5, min = 2, max = 30, step = 0.5),
              downloadButton("download_identification_overview_png", "Download overview PNG"),
              downloadButton("download_identification_overview_svg", "Download overview SVG"),
              downloadButton("download_run_identifications_png", "Download stacked PNG"),
              downloadButton("download_run_identifications_svg", "Download stacked SVG")
            ),
            column(
              width = 8,
              h4("Per-run identification count"),
              uiOutput("identification_overview_plot_ui"),
              h4("Identification coverage per run"),
              uiOutput("run_identifications_plot_ui"),
              verbatimTextOutput("run_identifications_note")
            )
          )
        )
  )
)

server <- function(input, output, session) {

  normalize_project_file_info <- function(file_info) {
    if (is.null(file_info)) return(NULL)
    value <- function(name, default = "") {
      if (is.data.frame(file_info) && name %in% colnames(file_info)) {
        return(as.character(file_info[[name]][1]))
      }
      if (!is.null(file_info[[name]])) return(as.character(file_info[[name]][1]))
      default
    }
    size_value <- if (is.data.frame(file_info) && "size" %in% colnames(file_info)) {
      suppressWarnings(as.numeric(file_info$size[1]))
    } else if (!is.null(file_info$size)) {
      suppressWarnings(as.numeric(file_info$size[1]))
    } else {
      NA_real_
    }
    data.frame(
      name = value("name"),
      size = size_value,
      type = value("type"),
      datapath = value("datapath"),
      stringsAsFactors = FALSE
    )
  }

  read_uploaded_table <- function(file_info) {
    file_info <- normalize_project_file_info(file_info)
    validate(need(!is.null(file_info) && file.exists(file_info$datapath), "The selected table file is not available. Reopen the project ZIP or upload the file again."))
    ext <- tolower(tools::file_ext(file_info$name))
    sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
    read.table(
      file_info$datapath,
      header = TRUE,
      sep = sep,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "\"",
      comment.char = ""
    )
  }

  read_sample_details_workbook <- function(file_info) {
    validate(need(requireNamespace("readxl", quietly = TRUE), "Package 'readxl' is required to read the sample-details workbook."))
    file_info <- normalize_project_file_info(file_info)
    validate(need(!is.null(file_info) && file.exists(file_info$datapath), "The selected sample-details workbook is not available. Reopen the project ZIP or upload the file again."))

    workbook_reader <- if (tolower(tools::file_ext(file_info$name)) == "xls") readxl::read_xls else readxl::read_xlsx
    unheaded_details <- as.data.frame(
      workbook_reader(file_info$datapath, col_names = FALSE),
      stringsAsFactors = FALSE
    )
    header_row <- which(
      as.character(unheaded_details[[1]]) == "ID" &
        as.character(unheaded_details[[2]]) == "Name"
    )
    header_row <- header_row[length(header_row)]
    validate(need(!is.na(header_row), "Could not find the ID / Name sample table in the sample-details workbook."))

    details <- as.data.frame(
      workbook_reader(file_info$datapath, skip = header_row - 1),
      stringsAsFactors = FALSE
    )
    validate(need(all(c("ID", "Name") %in% colnames(details)), "Sample-details workbook must contain ID and Name columns."))

    details <- details[!is.na(details$ID) & details$ID != "[end samples]", , drop = FALSE]
    details$ID <- as.character(details$ID)
    details$SampleDetailsID <- paste0("ID", details$ID)
    details$SampleName <- as.character(details$Name)
    details
  }

  build_export_filename <- function(filename) {
    if (is.null(filename) || !nzchar(filename)) filename <- "export"
    ext <- tools::file_ext(filename)
    stem <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", filename) else filename
    out <- if (isTRUE(input$append_timestamp)) {
      stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      if (nzchar(ext)) paste0(stem, "_", stamp, ".", ext) else paste0(stem, "_", stamp)
    } else {
      filename
    }
    if (isTRUE(input$include_project_number_in_downloads)) {
      out <- add_project_number_to_filename(out, project_number_for_downloads())
    }
    out
  }

  download_destination_path <- function() {
    path <- input$download_destination_dir
    if (is.null(path) || length(path) == 0) path <- ""
    path <- trimws(as.character(path)[1])
    if (!nzchar(path)) return(NA_character_)
    normalizePath(path, winslash = "/", mustWork = FALSE)
  }

  unique_destination_file <- function(destination_dir, filename) {
    candidate <- file.path(destination_dir, basename(filename))
    if (!file.exists(candidate)) return(candidate)
    ext <- tools::file_ext(candidate)
    stem <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", candidate) else candidate
    stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    if (nzchar(ext)) paste0(stem, "_", stamp, ".", ext) else paste0(stem, "_", stamp)
  }

  save_download_copy <- function(file, filename) {
    destination <- proteomics_download_copy_destination(input$save_download_copy, input$download_destination_dir)
    if (!isTRUE(destination$copy)) {
      if (nzchar(destination$message)) download_destination_message(destination$message)
      return(invisible(FALSE))
    }
    destination_dir <- destination$path
    destination_file <- unique_destination_file(destination_dir, filename)
    copied <- file.copy(file, destination_file, overwrite = FALSE)
    if (!isTRUE(copied) || !file.exists(destination_file)) {
      download_destination_message(paste0("Browser download completed; additional folder copy could not be written to: ", destination_file))
      return(invisible(FALSE))
    }
    download_destination_message(paste0("Saved copy: ", destination_file))
    invisible(TRUE)
  }

  project_file_ids <- c(
    "condition_setup_sample_details_file", "condition_setup_template_file",
    "evosep_manifest_file", "evosep_template_file",
    "meta_file", "run_order_file", "sample_details_file", "protein_no_impute_file",
    "protein_imputed_file", "cv_distribution_file", "data_file",
    "identifications_overview_file", "run_identifications_precursor_file",
    "run_identifications_protein_file"
  )
  project_file_labels <- c(
    condition_setup_sample_details_file = "Condition setup sample-details workbook",
    condition_setup_template_file = "Condition setup template",
    evosep_manifest_file = "Evosep queue manifest workbook",
    evosep_template_file = "Evosep queue CSL template",
    meta_file = "Table S1 metadata basis / condition setup",
    run_order_file = "Run-order table",
    sample_details_file = "Order sample-details workbook",
    protein_no_impute_file = "Table S2 protein report, no impute",
    protein_imputed_file = "Table S3 protein report, imputed",
    cv_distribution_file = "Uploaded CV distribution table",
    data_file = "Separate PCA abundance upload",
    identifications_overview_file = "Identifications overview table",
    run_identifications_precursor_file = "Run identifications, precursors",
    run_identifications_protein_file = "Run identifications, protein groups"
  )
  imported_project_files <- reactiveVal(list())
  restored_project_settings <- reactiveVal(list())
  project_restore_token <- reactiveVal(0L)
  project_bundle_message <- reactiveVal("No saved project is open. Upload source files directly or open a project ZIP.")
  download_destination_message <- reactiveVal("Browser downloads use your browser's download settings. Enable folder copies to save an additional copy from the app.")
  project_db_cache <- reactiveVal(list())
  ignored_project_upload_paths <- reactiveVal(list())
  pending_new_project_path <- reactiveVal("")
  generated_s3_result <- reactiveVal(NULL)
  active_metadata_override <- reactiveVal(NULL)
  metadata_replacement_message <- reactiveVal("No replacement metadata applied this session.")
  active_project_path <- reactiveVal("")
  active_project_display_name <- reactiveVal("")
  active_project_last_saved <- reactiveVal(NULL)
  active_project_last_error <- reactiveVal("")
  active_project_autosave_requested_at <- reactiveVal(NULL)
  active_project_autosave_pending_reason <- reactiveVal("")
  active_project_autosave_include_derived <- reactiveVal(FALSE)
  project_database_message <- reactiveVal("No active DuckDB/RDS project is open.")
  cached_batch_corrected_s3_result <- reactiveVal(NULL)
  batch_cache_message <- reactiveVal("No saved batch correction is loaded.")
spqc_metadata_edits <- reactiveVal(empty_spqc_metadata_edits())
spqc_metadata_draft_edits <- reactiveVal(empty_spqc_metadata_edits())
spqc_clear_pending <- reactiveVal(FALSE)
applied_metadata_state <- reactiveVal(empty_proteomics_metadata_state())
metadata_apply_revision <- reactiveVal(0L)
metadata_editor_revision <- reactiveVal(0L)
metadata_apply_message <- reactiveVal("Metadata changes have not been applied.")
metadata_draft_tracking_enabled <- reactiveVal(TRUE)
protein_stats_paused <- reactiveVal(FALSE)
  protein_stats_refresh_revision <- reactiveVal(0L)
  stopped_stats_tables <- reactiveVal(character(0))

  observeEvent(input$run_protein_stats, {
    protein_stats_paused(FALSE)
    protein_stats_refresh_revision(protein_stats_refresh_revision() + 1L)
    showNotification("Statistics recalculation requested.", type = "message", duration = 3)
  })

  project_file <- function(id) {
    imported <- imported_project_files()[[id]]
    if (!is.null(imported)) return(normalize_project_file_info(imported))
    direct <- normalize_project_file_info(input[[id]])
    if (proteomics_upload_is_ignored(direct, ignored_project_upload_paths()[[id]])) return(NULL)
    direct
  }

  current_project_upload_paths <- function() {
    paths <- lapply(project_file_ids, function(id) {
      file_info <- normalize_project_file_info(input[[id]])
      if (is.null(file_info)) "" else as.character(file_info$datapath[1L])
    })
    stats::setNames(paths, project_file_ids)
  }

  protein_source_table <- function(source) {
    no_impute <- proteomics_is_no_impute_source(source)
    if (no_impute) {
      file_info <- project_file("protein_no_impute_file")
      uploaded <- if (is.null(file_info)) NULL else read_uploaded_table(file_info)
      table <- resolve_proteomics_project_table(uploaded, project_db_cache(), "processed_s2")
    } else {
      s3_file <- project_file("protein_imputed_file")
      uploaded_s3 <- if (is.null(s3_file)) NULL else read_uploaded_table(s3_file)
      s2_file <- project_file("protein_no_impute_file")
      uploaded_s2 <- if (is.null(s2_file)) NULL else read_uploaded_table(s2_file)
      s2 <- resolve_proteomics_project_table(uploaded_s2, project_db_cache(), "processed_s2")
      generated <- generated_s3_result()
      generated <- if (is.null(generated)) NULL else generated$data
      method <- input$s3_imputation_source
      if (is.null(method) || !length(method)) method <- "spectronaut"
      cached_s3 <- project_db_cache()$processed_s3
      if (identical(method, "knn") && !is.null(s2_file) && is.null(generated)) cached_s3 <- NULL
      table <- resolve_protein_s3_source(method, uploaded_s3, cached_s3, generated, s2)
    }
    validate(need(!is.null(table), paste0("Upload ", if (no_impute) "Table S2" else "Table S3", " or open a DuckDB project containing it.")))
    table
  }

  protein_source_available <- function(source) {
    no_impute <- proteomics_is_no_impute_source(source)
    if (no_impute) return(!is.null(project_file("protein_no_impute_file")) || !is.null(project_db_cache()$processed_s2))
    method <- input$s3_imputation_source
    if (is.null(method) || !length(method)) method <- "spectronaut"
    if (identical(method, "s2")) return(protein_source_available("S2"))
    if (identical(method, "knn")) {
      return(!is.null(generated_s3_result()) || (is.null(project_file("protein_no_impute_file")) && !is.null(project_db_cache()$processed_s3)))
    }
    !is.null(project_file("protein_imputed_file")) || !is.null(project_db_cache()$processed_s3)
  }

  project_number_for_downloads <- function() {
    manual_project_number <- clean_project_number(input$download_project_number)
    if (!is.na(manual_project_number)) return(manual_project_number)
    priority_ids <- c(
      "meta_file", "run_order_file", "protein_no_impute_file", "protein_imputed_file",
      setdiff(project_file_ids, c("meta_file", "run_order_file", "protein_no_impute_file", "protein_imputed_file"))
    )
    source_names <- vapply(priority_ids, function(id) {
      source <- project_file(id)
      if (is.null(source) || is.null(source$name)) "" else as.character(source$name)
    }, character(1))
    infer_project_number(source_names)
  }

  project_file_status <- reactive({
    imported <- imported_project_files()
    cache <- project_db_cache()
    restored_cache_name <- c(
      meta_file = "metadata",
      protein_no_impute_file = "processed_s2",
      protein_imputed_file = "processed_s3"
    )
    data.frame(
      FileType = unname(project_file_labels[project_file_ids]),
      InputID = project_file_ids,
      Status = vapply(project_file_ids, function(id) {
        if (!is.null(imported[[id]])) {
          "Restored from project ZIP"
        } else if (!is.null(input[[id]])) {
          "Uploaded directly"
        } else if (id %in% names(restored_cache_name) && !is.null(cache[[restored_cache_name[[id]]]])) {
          paste0("Restored from DuckDB (", nrow(cache[[restored_cache_name[[id]]]]), " rows)")
        } else {
          "Not loaded"
        }
      }, character(1)),
      FileName = vapply(project_file_ids, function(id) {
        source <- project_file(id)
        if (!is.null(source) && !is.null(source$name) && nzchar(source$name)) {
          source$name
        } else if (id %in% names(restored_cache_name) && !is.null(cache[[restored_cache_name[[id]]]])) {
          paste0(restored_cache_name[[id]], " (database table)")
        } else ""
      }, character(1)),
      stringsAsFactors = FALSE
    )
  })

  lapply(project_file_ids, function(id) {
    observeEvent(input[[id]], {
    if (id %in% c("meta_file", "run_order_file", "sample_details_file")) {
      active_metadata_override(NULL)
      spqc_metadata_draft_edits(empty_spqc_metadata_edits())
      metadata_replacement_message("Source metadata changed; no replacement metadata is active.")
      }
      imported <- imported_project_files()
      if (!is.null(imported[[id]])) {
        imported[[id]] <- NULL
        imported_project_files(imported)
        project_bundle_message(paste0("Opened project is active; replaced project source file: ", input[[id]]$name, "."))
      }
      cache <- project_db_cache()
      if (length(cache)) {
        if (id %in% c("meta_file", "run_order_file", "sample_details_file")) {
          metadata_apply_message("Unapplied metadata changes.")
          project_database_message(paste0("Metadata source loaded into draft: ", input[[id]]$name, ". Click Apply metadata changes to update subsequent tabs."))
        } else {
          project_db_cache(invalidate_proteomics_project_cache(cache, id))
          cached_batch_corrected_s3_result(NULL)
          if (id %in% c("condition_setup_sample_details_file", "condition_setup_template_file")) {
            project_database_message(paste0("Condition-setup source replaced: ", input[[id]]$name, "."))
          } else if (id %in% c("protein_no_impute_file", "protein_imputed_file")) {
            project_database_message(paste0("Protein source replaced: ", input[[id]]$name, ". Dependent processed results were invalidated."))
          }
        }
      }
    }, ignoreInit = TRUE)
  })

  project_settings <- function() {
    ids <- c(
        "project_include_readme", "project_include_session_info",
        "download_project_number", "include_project_number_in_downloads",
        "save_download_copy", "download_destination_dir",
        "condition_setup_run_label_suffix", "condition_setup_condition_col",
        "condition_setup_replicate_order_col",
        "evosep_source_tray", "evosep_output_dir", "evosep_filename_suffix",
        "evosep_comment", "evosep_adh_positions", "evosep_spqc_positions",
        "evosep_sample_positions", "evosep_spqc_prefix", "evosep_analysis_method",
        "evosep_xcalibur_method",
        "metadata_export_columns", "metadata_workbook_filename", "s2_non_data_columns",
      "s3_non_data_columns", "protein_header_label_mode", "protein_header_label_columns",
      "protein_quantity_order_columns", "s3_imputation_source", "protein_knn_k",
      "protein_knn_max_missing_percent", "protein_knn_scope", "protein_knn_group_col",
      "spqc_assignment_mode", "spqc_group_label", "spqc_group_prefix", "spqc_batch_overrides",
      "sample_exclusions_text", "sample_exclusion_reason", "protein_cv_group_col", "cv_conditions",
      "protein_derived_column_order", "stats_tables", "stats_group_columns", "stats_comparisons", "stats_paired_comparisons", "stats_pair_col",
      "stats_bh_fdr", "numeric_sig_figs", "scientific_small_values",
      "scientific_threshold", "enable_protein_filters", "cv_plot_source",
      "batch_correction_source", "batch_correction_batch_col", "batch_correction_group_col",
      "batch_correction_combat_mode", "batch_correction_pseudocount", "batch_correction_min_batches_for_feature",
      "cv_plot_group_col", "cv_plot_conditions", "cv_density_adjust", "cv_x_cutoff", "cv_plot_title",
      "cv_fill_density", "cv_title_size", "cv_axis_title_size", "cv_axis_text_size",
      "cv_legend_text_size", "cv_median_text_size", "cv_line_width",
      "cv_figure_width", "cv_figure_height", "pca_data_source",
      "data_layout", "auto_detect_annotations", "annotation_rows", "report_feature_col",
      "pca_min_observed_percent", "zscore", "max_ncp", "use_estimated_ncp",
      "manual_ncp", "sample_id_col", "label_mode", "color_by", "shape_by",
      "label_points", "label_missing_points", "point_size", "label_size", "plot_title", "plot_subtitle",
      "pca_figure_width", "pca_figure_height",
      "clustvis_pca_source", "clustvis_pca_subset_column", "clustvis_pca_subset_values",
      "clustvis_pca_min_observed_percent", "clustvis_pca_npcs",
      "clustvis_pca_scale", "clustvis_pca_ellipses", "clustvis_pca_color_by",
      "clustvis_pca_shape_by", "clustvis_pca_label_by", "clustvis_pca_point_size",
      "clustvis_pca_label_size", "clustvis_pca_title", "clustvis_pca_width",
      "clustvis_pca_height", "clustvis_pca_opacity_override_groups",
      "clustvis_pca_no_fill_shapes", "clustvis_pca_interactive", "clustvis_pca_override_opacity",
      "clustvis_pca_default_opacity", "pca_loading_rank_by", "pca_loading_top_n",
      "append_timestamp", "png_filename", "svg_filename", "volcano_source",
      "volcano_comparison", "volcano_significance_metric", "volcano_fc_cutoff",
      "volcano_sig_cutoff", "volcano_label_col", "volcano_label_mode",
      "volcano_max_labels", "volcano_title", "volcano_point_size",
      "volcano_label_size", "volcano_figure_width", "volcano_figure_height",
      "volcano_interactive", "gsea_source", "gsea_comparison", "gsea_rank_metric",
      "gsea_gene_col", "gsea_species", "gsea_collection", "gsea_min_size",
      "gsea_max_size", "gsea_top_n", "gsea_plot_width", "gsea_plot_height",
      "feature_data_source", "feature_select", "feature_order_columns", "feature_group_by", "feature_label_by",
      "feature_plot_title", "feature_title_size", "feature_bar_width",
      "feature_text_size", "feature_value_scale", "feature_color_mode", "feature_symmetric_scale",
      "feature_group_style", "show_feature_mean", "rotate_feature_labels",
      "feature_interactive", "feature_figure_width", "feature_figure_height", "feature_plot_ncol",
      "box_group_by", "box_plot_style", "box_value_scale", "box_plot_title",
      "box_y_axis_title", "box_point_size", "box_text_size",
      "box_figure_width", "box_figure_height",
      "script_box_source", "script_box_features", "script_box_group_by", "script_box_label_by",
      "script_box_conditions", "script_box_value_scale", "script_box_plot_style", "script_box_title",
      "script_box_y_axis_title", "script_box_point_size", "script_box_point_opacity", "script_box_label_size", "script_box_text_size",
      "script_box_width", "script_box_height", "script_box_ncol",
      "correlation_source", "correlation_feature", "correlation_value_scale",
      "correlation_method", "correlation_group_by", "correlation_groups",
      "correlation_covariates", "correlation_rank_by", "correlation_top_n", "correlation_exclude_reference",
      "correlation_plot_title", "correlation_figure_width", "correlation_figure_height",
      "identification_metric", "identification_overview_title",
      "run_identifications_title", "identification_title_size",
      "identification_axis_text_size", "identification_legend_size",
      "identification_overview_width", "identification_overview_height",
      "run_identifications_width", "run_identifications_height"
    )
    settings <- stats::setNames(lapply(ids, function(id) input[[id]]), ids)
    settings$workflow_tabs <- input$workflow_tabs
    settings
  }

  batch_correction_cache_signature <- function() {
    validate(need(requireNamespace("jsonlite", quietly = TRUE), "Package 'jsonlite' is required to cache batch correction results."))
    file_signatures <- lapply(c("meta_file", "run_order_file", "sample_details_file", "protein_imputed_file"), function(id) {
      source <- project_file(id)
      if (is.null(source) || !file.exists(source$datapath)) return(NULL)
      integrity <- file_integrity(source$datapath)
      list(
        input_id = id,
        name = source$name,
        size_bytes = integrity$size_bytes,
        md5 = integrity$md5
      )
    })
    names(file_signatures) <- c("meta_file", "run_order_file", "sample_details_file", "protein_imputed_file")
    jsonlite::toJSON(
      list(
        artifact_type = "S3_batch_corrected_result",
        app_name = "Proteomics Data Workup",
        files = file_signatures,
        settings = project_settings()
      ),
      auto_unbox = TRUE,
      null = "null"
    )
  }

  restore_project_settings <- function(settings) {
    pair_col <- if (is.null(settings$stats_pair_col)) "Replicate" else as.character(settings$stats_pair_col)[1L]
    updateSelectInput(session, "stats_pair_col", choices = unique(c("Replicate", pair_col)), selected = pair_col)
    if (is.null(settings$script_box_features) && !is.null(settings$script_box_feature)) {
      settings$script_box_features <- settings$script_box_feature
    }
    update_text <- c(
        "metadata_workbook_filename", "report_feature_col", "sample_id_col", "plot_title",
        "plot_subtitle", "png_filename", "svg_filename", "volcano_title",
        "feature_plot_title", "identification_overview_title", "run_identifications_title",
        "cv_plot_title", "box_plot_title", "box_y_axis_title",
        "clustvis_pca_title", "script_box_title", "script_box_y_axis_title",
        "correlation_plot_title",
        "download_project_number", "download_destination_dir",
        "spqc_group_label", "spqc_group_prefix", "spqc_batch_overrides",
        "sample_exclusions_text", "sample_exclusion_reason",
        "condition_setup_run_label_suffix", "evosep_source_tray", "evosep_output_dir",
        "evosep_filename_suffix", "evosep_comment", "evosep_adh_positions",
        "evosep_spqc_positions", "evosep_sample_positions", "evosep_spqc_prefix",
        "evosep_analysis_method", "evosep_xcalibur_method"
    )
    update_numeric <- c(
      "numeric_sig_figs", "scientific_threshold", "cv_density_adjust", "cv_x_cutoff",
      "cv_title_size", "cv_axis_title_size", "cv_axis_text_size", "cv_legend_text_size",
      "cv_median_text_size", "cv_line_width", "annotation_rows", "pca_min_observed_percent",
      "max_ncp", "manual_ncp", "point_size", "label_size", "volcano_fc_cutoff",
      "volcano_sig_cutoff", "volcano_max_labels", "volcano_point_size",
      "volcano_label_size", "feature_title_size", "feature_bar_width",
      "feature_text_size", "identification_title_size", "identification_axis_text_size",
      "identification_legend_size", "cv_figure_width", "cv_figure_height",
      "pca_figure_width", "pca_figure_height", "volcano_figure_width",
      "volcano_figure_height", "feature_figure_width", "feature_figure_height", "feature_plot_ncol",
      "box_point_size", "box_text_size", "box_figure_width", "box_figure_height",
      "clustvis_pca_min_observed_percent", "clustvis_pca_npcs", "clustvis_pca_point_size",
      "clustvis_pca_label_size", "clustvis_pca_width", "clustvis_pca_height",
      "script_box_point_size", "script_box_point_opacity", "script_box_label_size", "script_box_text_size", "script_box_width", "script_box_height", "script_box_ncol",
      "correlation_top_n", "correlation_figure_width", "correlation_figure_height",
      "identification_overview_width", "identification_overview_height",
      "run_identifications_width", "run_identifications_height",
      "pca_loading_top_n", "gsea_min_size", "gsea_max_size", "gsea_top_n",
      "gsea_plot_width", "gsea_plot_height", "clustvis_pca_override_opacity",
      "clustvis_pca_default_opacity", "batch_correction_combat_mode",
      "batch_correction_pseudocount", "batch_correction_min_batches_for_feature",
      "protein_knn_k", "protein_knn_max_missing_percent"
    )
    update_checkbox <- c(
      "project_include_readme", "project_include_session_info",
      "include_project_number_in_downloads", "save_download_copy",
      "stats_bh_fdr", "scientific_small_values", "enable_protein_filters",
      "cv_fill_density", "auto_detect_annotations", "zscore", "use_estimated_ncp",
      "label_points", "label_missing_points", "append_timestamp", "volcano_interactive", "feature_symmetric_scale",
      "show_feature_mean", "rotate_feature_labels", "clustvis_pca_scale", "clustvis_pca_ellipses",
      "clustvis_pca_no_fill_shapes", "clustvis_pca_interactive", "evosep_fill_manifest_order",
      "correlation_exclude_reference", "feature_interactive"
    )
    update_radio <- c(
      "cv_plot_source", "pca_data_source", "volcano_source", "volcano_label_mode",
      "feature_data_source", "feature_value_scale", "box_value_scale", "clustvis_pca_source", "script_box_source",
      "script_box_value_scale", "script_box_plot_style", "protein_header_label_mode",
      "protein_derived_column_order", "spqc_assignment_mode", "gsea_source", "gsea_rank_metric",
      "gsea_species", "gsea_collection", "batch_correction_source",
      "correlation_source", "correlation_value_scale", "correlation_method",
      "s3_imputation_source", "protein_knn_scope"
    )
    update_select <- c(
        "data_layout", "label_mode", "color_by", "shape_by", "volcano_comparison",
        "volcano_significance_metric", "volcano_label_col", "feature_color_mode",
        "feature_group_style", "feature_group_by", "feature_label_by", "box_group_by", "box_plot_style", "identification_metric",
        "clustvis_pca_color_by", "clustvis_pca_shape_by", "clustvis_pca_label_by", "clustvis_pca_subset_column",
        "pca_loading_rank_by", "gsea_comparison", "gsea_gene_col",
        "correlation_rank_by", "correlation_group_by",
        "script_box_label_by", "batch_correction_batch_col", "batch_correction_group_col",
        "condition_setup_condition_col", "condition_setup_replicate_order_col", "protein_cv_group_col", "cv_plot_group_col",
        "protein_knn_group_col"
    )
    update_selectize <- c(
      "s2_non_data_columns", "s3_non_data_columns",
      "cv_plot_conditions", "clustvis_pca_subset_values",
      "feature_select", "script_box_features",
      "clustvis_pca_opacity_override_groups", "correlation_feature",
      "correlation_groups", "correlation_covariates"
    )
    for (id in update_text) if (!is.null(settings[[id]])) updateTextInput(session, id, value = settings[[id]])
    for (id in update_numeric) if (!is.null(settings[[id]])) updateNumericInput(session, id, value = settings[[id]])
    for (id in update_checkbox) if (!is.null(settings[[id]])) updateCheckboxInput(session, id, value = settings[[id]])
    if (!is.null(settings$stats_tables)) updateCheckboxGroupInput(session, "stats_tables", selected = unlist(settings$stats_tables))
    for (id in update_radio) if (!is.null(settings[[id]])) updateRadioButtons(session, id, selected = settings[[id]])
    for (id in update_select) if (!is.null(settings[[id]])) updateSelectInput(session, id, selected = settings[[id]])
    for (id in update_selectize) if (!is.null(settings[[id]])) updateSelectizeInput(session, id, selected = unlist(settings[[id]]), server = TRUE)
    if (!is.null(settings$workflow_tabs)) {
      session$onFlushed(function() updateTabsetPanel(session, "workflow_tabs", selected = settings$workflow_tabs), once = TRUE)
    }
  }

  write_project_bundle <- function(file) {
    validate(need(requireNamespace("jsonlite", quietly = TRUE) && requireNamespace("zip", quietly = TRUE), "Packages 'jsonlite' and 'zip' are required to export a project ZIP."))
    bundle_dir <- tempfile("proteomics_project_")
    dir.create(bundle_dir)
    on.exit(unlink(bundle_dir, recursive = TRUE, force = TRUE), add = TRUE)
    progress <- shiny::Progress$new(session, min = 0, max = length(project_file_ids) + 2)
    progress$set(message = "Preparing project ZIP", value = 0)
    on.exit(progress$close(), add = TRUE)
    file_manifest <- list()
    for (index in seq_along(project_file_ids)) {
      id <- project_file_ids[index]
      source <- project_file(id)
      if (!is.null(source)) {
        safe_basename <- gsub("[^A-Za-z0-9._-]", "_", basename(source$name))
        stored_name <- paste0(id, "_", safe_basename)
        destination <- file.path(bundle_dir, stored_name)
        copied <- file.copy(source$datapath, destination, overwrite = TRUE, copy.date = TRUE)
        validate(need(isTRUE(copied) && file.exists(destination), paste0("Could not copy project source file: ", source$name)))
        integrity <- file_integrity(destination)
        file_manifest[[id]] <- list(
          original_name = source$name,
          stored_name = stored_name,
          size_bytes = integrity$size_bytes,
          md5 = integrity$md5
        )
      }
      progress$set(value = index, detail = paste("Checked", index, "of", length(project_file_ids), "file slots"))
    }
    manifest <- list(
      format = "Proteomics Data Workup project",
      version = 2,
      app_schema_version = 2,
      exported_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      included_file_count = length(file_manifest),
      files = file_manifest,
      settings = project_settings()
    )
    jsonlite::write_json(manifest, file.path(bundle_dir, "project_manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
    format_value <- function(value) {
      if (is.null(value) || length(value) == 0) return("Not set")
      value <- unlist(value, use.names = FALSE)
      value <- value[!is.na(value) & nzchar(as.character(value))]
      if (length(value) == 0) return("Not set")
      paste(as.character(value), collapse = ", ")
    }
    comparison_labels <- function(value) {
      if (is.null(value) || length(value) == 0) return("Not set")
      labels <- vapply(unlist(value, use.names = FALSE), stats_comparison_label, character(1))
      labels <- labels[!is.na(labels) & nzchar(labels)]
      if (length(labels) == 0) return("Not set")
      paste(labels, collapse = ", ")
    }
    bool_label <- function(value) {
      if (isTRUE(value)) "Yes" else "No"
    }
    file_label <- function(id) {
      if (!is.null(file_manifest[[id]])) file_manifest[[id]]$original_name else "Not included"
    }
    source_lines <- unlist(lapply(names(project_file_labels), function(id) {
      paste0("- ", project_file_labels[[id]], ": ", file_label(id))
    }), use.names = FALSE)
    settings <- manifest$settings
    summary_lines <- c(
      "Proteomics Data Workup project summary",
      paste0("Exported: ", manifest$exported_at),
      "",
      "Source files",
      source_lines,
      "",
      "Table S1. Metadata",
      paste0("- Columns and order: ", format_value(settings$metadata_export_columns)),
      paste0("- Workbook filename: ", format_value(settings$metadata_workbook_filename)),
      paste0("- Metadata sample ID column: ", format_value(settings$sample_id_col)),
      "",
      "Protein Tables",
      paste0("- Table S3 source: ", format_value(settings$s3_imputation_source)),
      paste0("- kNN nearest proteins: ", format_value(settings$protein_knn_k)),
      paste0("- kNN maximum missingness (%): ", format_value(settings$protein_knn_max_missing_percent)),
      paste0("- kNN scope: ", format_value(settings$protein_knn_scope)),
      paste0("- kNN metadata field: ", format_value(settings$protein_knn_group_col)),
      paste0("- Table S2 non-data columns and order: ", format_value(settings$s2_non_data_columns)),
      paste0("- Table S3 non-data columns and order: ", format_value(settings$s3_non_data_columns)),
      paste0("- Measurement header variables: ", format_value(settings$protein_header_label_columns)),
      paste0("- Protein quantity column order: ", format_value(settings$protein_quantity_order_columns)),
      paste0("- CV conditions: ", format_value(settings$cv_conditions)),
      paste0("- Derived column order: ", if (identical(settings$protein_derived_column_order, "stats_before_cv")) "Statistics before CV columns" else "CV columns before statistics"),
      paste0("- Statistics tables: ", format_value(settings$stats_tables)),
      paste0("- Statistics comparisons: ", comparison_labels(settings$stats_comparisons)),
      paste0("- Paired statistics comparisons: ", comparison_labels(settings$stats_paired_comparisons)),
      paste0("- Include Benjamini-Hochberg FDR: ", bool_label(settings$stats_bh_fdr)),
      paste0("- Significant figures for derived numbers: ", format_value(settings$numeric_sig_figs)),
      paste0("- Scientific notation for small derived values: ", bool_label(settings$scientific_small_values)),
      paste0("- Scientific notation threshold: ", format_value(settings$scientific_threshold)),
      paste0("- Excel filters on S2/S3: ", bool_label(settings$enable_protein_filters)),
      "",
      "Batch Correction",
      paste0("- Source: ", format_value(settings$batch_correction_source)),
      paste0("- Batch metadata column: ", format_value(settings$batch_correction_batch_col)),
      paste0("- Biological group preserved: ", format_value(settings$batch_correction_group_col)),
      paste0("- ComBat mode: ", format_value(settings$batch_correction_combat_mode)),
      paste0("- Log2 pseudocount: ", format_value(settings$batch_correction_pseudocount)),
      paste0("- Minimum batches for protein retention: ", format_value(settings$batch_correction_min_batches_for_feature)),
      "",
      "CV Plot",
      paste0("- Data source: ", format_value(settings$cv_plot_source)),
      paste0("- Conditions plotted: ", format_value(settings$cv_plot_conditions)),
      paste0("- Density smoothness: ", format_value(settings$cv_density_adjust)),
      paste0("- X-axis cutoff: ", format_value(settings$cv_x_cutoff)),
      paste0("- Fill density areas: ", bool_label(settings$cv_fill_density)),
      "",
      "PCA",
      paste0("- Abundance source: ", format_value(settings$pca_data_source)),
      paste0("- Feature label column: ", format_value(settings$report_feature_col)),
      paste0("- Minimum observed samples per protein (%): ", format_value(settings$pca_min_observed_percent)),
      paste0("- Z-score features before PCA: ", bool_label(settings$zscore)),
      paste0("- Estimate optimal ncp: ", bool_label(settings$use_estimated_ncp)),
      paste0("- Max ncp to test: ", format_value(settings$max_ncp)),
      paste0("- Manual ncp: ", format_value(settings$manual_ncp)),
      paste0("- Label points by: ", format_value(settings$label_mode)),
      paste0("- Color by: ", format_value(settings$color_by)),
      paste0("- Shape by: ", format_value(settings$shape_by)),
      "",
      "Volcano Plot",
      paste0("- Statistics source: ", format_value(settings$volcano_source)),
      paste0("- Comparison: ", comparison_labels(settings$volcano_comparison)),
      paste0("- Significance metric: ", if (identical(settings$volcano_significance_metric, "BH_FDR")) "BH FDR" else format_value(settings$volcano_significance_metric)),
      paste0("- Absolute log2 fold-change cutoff: ", format_value(settings$volcano_fc_cutoff)),
      paste0("- Significance cutoff: ", format_value(settings$volcano_sig_cutoff)),
      paste0("- Label proteins by: ", format_value(settings$volcano_label_col)),
      paste0("- Label mode: ", format_value(settings$volcano_label_mode)),
      paste0("- Top labels per direction: ", format_value(settings$volcano_max_labels)),
      paste0("- Interactive volcano enabled: ", bool_label(settings$volcano_interactive)),
      "",
      "Feature Plot and Box Plot",
      paste0("- Feature abundance source: ", format_value(settings$feature_data_source)),
      paste0("- Selected feature: ", format_value(settings$feature_select)),
      paste0("- Feature plot value scale: ", format_value(settings$feature_value_scale)),
      paste0("- Feature color mode: ", format_value(settings$feature_color_mode)),
      paste0("- Box plot group by: ", format_value(settings$box_group_by)),
      paste0("- Box plot value scale: ", format_value(settings$box_value_scale)),
      "",
      "Run Identifications",
      paste0("- Identification metric: ", format_value(settings$identification_metric)),
      "",
      "Machine-readable settings are also saved in project_manifest.json."
    )
    writeLines(summary_lines, file.path(bundle_dir, "project_summary_README.txt"))
    if (isTRUE(input$project_include_readme)) {
      required_packages <- c("shiny", "ggplot2", "plotly", "htmlwidgets", "DT", "dplyr", "stringr", "missMDA", "FactoMineR", "svglite", "readxl", "openxlsx", "jsonlite", "DBI", "duckdb", "zip", "msigdbr", "BiocManager")
      optional_packages <- c("fgsea", "HarmonizR")
      install_line <- paste0("install.packages(c(", paste(sprintf("\"%s\"", required_packages), collapse = ", "), "))")
      optional_line <- paste0("install.packages(c(", paste(sprintf("\"%s\"", optional_packages), collapse = ", "), "))")
      writeLines(
        c(
          "Proteomics Data Workup project bundle",
          "",
          "Contents",
          "- Uploaded source data files used by the app",
          "- project_manifest.json with saved app settings and original file names",
          "- Optional sessionInfo.txt with the exporting user's R/package versions",
          "",
          "Install required R packages",
          install_line,
          "",
          "Install optional packages for advanced tabs",
          optional_line,
          "BiocManager::install(\"fgsea\")",
          "",
          "Batch correction note",
          "The Batch correction tab uses HarmonizR / ComBat and defaults to Table S3 imputed protein abundance data.",
          "",
          "Open the app",
          "shiny::runApp(\"path/to/app_pca_svd_with_feature_plot_adjustable_title.R\")",
          "",
          "Restore this project",
          "1. Open the Make metadata tab.",
          "2. Select this ZIP under 'Select saved Proteomics Data Workup project'.",
          "3. Click 'Open project ZIP'.",
          "",
          "Sharing note",
          "This ZIP includes uploaded study data. Confirm that sharing the data is permitted before sending it to another user.",
          "",
          "Why no .RData environment file?",
          "A saved R workspace is not needed to reopen this project. The source files and settings are portable and are less dependent on a particular R session or computer."
        ),
        file.path(bundle_dir, "README.txt")
      )
    }
    if (isTRUE(input$project_include_session_info)) {
      capture.output(utils::sessionInfo(), file = file.path(bundle_dir, "sessionInfo.txt"))
    }
    if (proteomics_duckdb_available()) {
      cache_path <- file.path(bundle_dir, "project_cache.duckdb")
      validate(need(save_active_project_state("ZIP project cache", include_derived = TRUE, destination = cache_path), "Could not create the embedded DuckDB project cache."))
    }
    temporary_zip <- tempfile("proteomics_project_validated_", fileext = ".zip")
    on.exit(unlink(temporary_zip, force = TRUE), add = TRUE)
    zip::zipr(temporary_zip, list.files(bundle_dir, recursive = TRUE, full.names = TRUE), root = bundle_dir)
    progress$set(value = length(project_file_ids) + 1, detail = "Validating project ZIP")
    listing <- utils::unzip(temporary_zip, list = TRUE)
    validate(need(!anyDuplicated(listing$Name), "Created project ZIP contains duplicate entries."))
    validate(need(all(vapply(listing$Name, safe_archive_entry, logical(1))), "Created project ZIP contains an unsafe path."))
    validate(need("project_manifest.json" %in% listing$Name, "Created project ZIP is missing its manifest."))
    for (id in names(file_manifest)) {
      entry <- file_manifest[[id]]
      row <- listing[listing$Name == entry$stored_name, , drop = FALSE]
      validate(need(nrow(row) == 1, paste0("Created project ZIP is missing: ", entry$stored_name)))
      validate(need(as.numeric(row$Length) == as.numeric(entry$size_bytes), paste0("Created project ZIP has the wrong size for: ", entry$stored_name)))
    }
    copied_zip <- file.copy(temporary_zip, file, overwrite = TRUE)
    validate(need(isTRUE(copied_zip) && file.exists(file), "Could not finalize the validated project ZIP."))
    progress$set(value = length(project_file_ids) + 2, detail = "Project ZIP ready")
  }

  load_project_bundle_file <- function(file_info) {
    validate(need(requireNamespace("jsonlite", quietly = TRUE), "Package 'jsonlite' is required to open a project ZIP."))
    import_dir <- tempfile("proteomics_project_import_")
    dir.create(import_dir)
    progress <- shiny::Progress$new(session, min = 0, max = 4)
    progress$set(message = "Opening project ZIP", value = 0, detail = "Inspecting archive")
    on.exit(progress$close(), add = TRUE)
    metadata_draft_tracking_enabled(FALSE)
    project_db_cache(list())
    applied_metadata_state(empty_proteomics_metadata_state())
    spqc_metadata_draft_edits(empty_spqc_metadata_edits())
    spqc_clear_pending(FALSE)
    listing <- tryCatch(utils::unzip(file_info$datapath, list = TRUE), error = function(error) NULL)
    validate(need(!is.null(listing) && nrow(listing) > 0, "Project ZIP is unreadable or empty."))
    archive_entries <- as.character(listing$Name)
    validate(need(!anyDuplicated(archive_entries), "Project ZIP contains duplicate archive entries and was not opened."))
    validate(need(all(vapply(archive_entries, safe_archive_entry, logical(1))), "Project ZIP contains an unsafe file path and was not opened."))
    validate(need("project_manifest.json" %in% archive_entries, "Project ZIP does not contain project_manifest.json."))
    utils::unzip(file_info$datapath, files = "project_manifest.json", exdir = import_dir)
    manifest_path <- file.path(import_dir, "project_manifest.json")
    validate(need(file.exists(manifest_path), "Project ZIP does not contain project_manifest.json."))
    manifest <- tryCatch(jsonlite::read_json(manifest_path, simplifyVector = FALSE), error = function(error) NULL)
    validate(need(!is.null(manifest), "Project manifest is malformed JSON."))
    validate(need(identical(manifest$format, "Proteomics Data Workup project"), "This ZIP is not a Proteomics Data Workup project."))
    manifest_version <- suppressWarnings(as.integer(unlist(manifest$version, use.names = FALSE)[1]))
    validate(need(!is.na(manifest_version) && manifest_version %in% c(1L, 2L), paste0("Unsupported project manifest version: ", manifest_version, ".")))
    manifest_files <- manifest$files
    if (is.null(manifest_files)) manifest_files <- list()
    unknown_ids <- setdiff(names(manifest_files), project_file_ids)
    validate(need(!length(unknown_ids), paste0("Project manifest contains unknown file types: ", paste(unknown_ids, collapse = ", "), ".")))

    candidate_entries <- list()
    stored_names <- character(0)
    for (id in names(manifest_files)) {
      entry <- manifest_files[[id]]
      stored_name <- as.character(unlist(entry$stored_name, use.names = FALSE)[1])
      original_name <- as.character(unlist(entry$original_name, use.names = FALSE)[1])
      validate(need(safe_archive_entry(stored_name), paste0("Project manifest contains an unsafe stored filename for ", id, ".")))
      validate(need(!stored_name %in% stored_names, paste0("Project manifest repeats stored filename: ", stored_name, ".")))
      archive_row <- listing[archive_entries == stored_name, , drop = FALSE]
      validate(need(nrow(archive_row) == 1, paste0("Project data file is missing: ", stored_name)))
      if (manifest_version >= 2L) {
        expected_size <- suppressWarnings(as.numeric(unlist(entry$size_bytes, use.names = FALSE)[1]))
        expected_md5 <- as.character(unlist(entry$md5, use.names = FALSE)[1])
        validate(need(is.finite(expected_size) && expected_size >= 0, paste0("Project manifest has an invalid size for: ", stored_name)))
        validate(need(nzchar(expected_md5), paste0("Project manifest has no checksum for: ", stored_name)))
        validate(need(as.numeric(archive_row$Length) == expected_size, paste0("Project archive size does not match its manifest for: ", stored_name)))
      } else {
        expected_size <- as.numeric(archive_row$Length)
        expected_md5 <- NA_character_
      }
      candidate_entries[[id]] <- list(
        original_name = original_name,
        stored_name = stored_name,
        expected_size = expected_size,
        expected_md5 = expected_md5
      )
      stored_names <- c(stored_names, stored_name)
    }

    progress$set(value = 1, detail = paste("Extracting", length(stored_names), "project files"))
    if (length(stored_names)) utils::unzip(file_info$datapath, files = stored_names, exdir = import_dir)
    restored_files <- list()
    verified_count <- 0L
    for (id in names(candidate_entries)) {
      entry <- candidate_entries[[id]]
      source_path <- file.path(import_dir, entry$stored_name)
      validate(need(file.exists(source_path), paste0("Project data file could not be extracted: ", entry$stored_name)))
      integrity <- file_integrity(source_path)
      validate(need(as.numeric(integrity$size_bytes) == as.numeric(entry$expected_size), paste0("Extracted project file has the wrong size: ", entry$stored_name)))
      if (manifest_version >= 2L) {
        validate(need(identical(tolower(integrity$md5), tolower(entry$expected_md5)), paste0("Project file checksum failed: ", entry$stored_name)))
        verified_count <- verified_count + 1L
      }
      restored_files[[id]] <- normalize_project_file_info(list(
        name = entry$original_name,
        datapath = source_path,
        size = file.info(source_path)$size,
        type = ""
      ))
    }
    progress$set(value = 2, detail = "Committing verified project files")
    imported_project_files(restored_files)
    settings <- manifest$settings
    if (is.null(settings)) settings <- list()
    if ("project_cache.duckdb" %in% archive_entries && proteomics_duckdb_available()) {
      utils::unzip(file_info$datapath, files = "project_cache.duckdb", exdir = import_dir)
      embedded <- load_proteomics_project_file(file.path(import_dir, "project_cache.duckdb"))
      project_db_cache(embedded)
      if (!is.null(embedded$settings)) settings <- embedded$settings
      spqc_metadata_draft_edits(empty_spqc_metadata_edits())
      if (!is.null(embedded$spqc_metadata_edits)) spqc_metadata_edits(embedded$spqc_metadata_edits)
    }
    restored_project_settings(settings)
    embedded_metadata <- project_db_cache()$metadata
    if (!is.null(embedded_metadata)) {
      initialize_applied_metadata(embedded_metadata, settings, spqc_metadata_edits())
    } else {
      session$onFlushed(function() {
        tryCatch(
          initialize_applied_metadata(
            isolate(draft_metadata()),
            settings,
            isolate(spqc_metadata_edits())
          ),
          error = function(e) metadata_apply_message(paste0("Restored metadata requires Apply: ", conditionMessage(e)))
        )
      }, once = TRUE)
    }
    project_restore_token(isolate(project_restore_token()) + 1L)
    restore_project_settings(settings)
    session$onFlushed(function() {
      session$onFlushed(function() {
        metadata_draft_tracking_enabled(TRUE)
        applied_state <- isolate(applied_metadata_state())
        if (!is.null(applied_state$metadata)) metadata_apply_message("All metadata changes applied.")
      }, once = TRUE)
    }, once = TRUE)
    restore_stage <- function(stage) {
      session$onFlushed(function() {
        restore_project_settings(settings)
        if (stage < 2L) {
          restore_stage(stage + 1L)
        } else {
          isolate(sync_restored_project_inputs(settings))
        }
      }, once = TRUE)
    }
    restore_stage(1L)
    progress$set(value = 3, detail = "Restoring project settings")
    restored_labels <- unname(project_file_labels[names(restored_files)])
    missing_ids <- setdiff(project_file_ids, names(restored_files))
    missing_labels <- unname(project_file_labels[missing_ids])
    project_bundle_message(paste0(
      "Opened project ZIP (manifest version ", manifest_version, ") with ", length(restored_files), " source data files restored.\n",
      if (manifest_version >= 2L) paste0("Checksum verified: ", verified_count, " files.\n") else "Legacy version-1 bundle: files were checked for presence and size; checksums were not available.\n",
      "Restored: ", if (length(restored_labels)) paste(restored_labels, collapse = "; ") else "none", ".\n",
      "Not included in this ZIP: ", if (length(missing_labels)) paste(missing_labels, collapse = "; ") else "none", ".\n",
      "Note: restored files are active in the app, but Shiny does not repopulate file-upload boxes. New individual uploads will override restored project files."
    ))
    progress$set(value = 4, detail = "Project ready")
    invisible(restored_files)
  }

  observeEvent(input$load_project_bundle, {
    req(input$project_bundle_file)
    load_project_bundle_file(input$project_bundle_file)
  })

  active_project_path_from_input <- function() {
    path <- active_project_path()
    if (is.null(path) || !length(path) || !nzchar(path)) stop("Open or create a project first.", call. = FALSE)
    path
  }

  initialize_applied_metadata <- function(metadata, settings = list(), edits = empty_spqc_metadata_edits()) {
    if (is.null(metadata) || !nrow(metadata)) {
      applied_metadata_state(empty_proteomics_metadata_state())
      metadata_apply_message("Metadata changes have not been applied.")
      return(invisible(FALSE))
    }
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
    columns <- unlist(settings$metadata_export_columns, use.names = FALSE)
    columns <- columns[columns %in% colnames(metadata)]
    if (!length(columns)) columns <- metadata_default_columns(colnames(metadata))
    if (!length(columns)) columns <- colnames(metadata)
    state <- build_proteomics_applied_state(metadata, columns, edits, empty_spqc_metadata_edits())
    applied_metadata_state(state)
    spqc_metadata_edits(state$spqc_metadata_edits)
    spqc_metadata_draft_edits(empty_spqc_metadata_edits())
    spqc_clear_pending(FALSE)
    metadata_apply_revision(isolate(metadata_apply_revision()) + 1L)
    metadata_apply_message("All metadata changes applied.")
    invisible(TRUE)
  }

  output$project_path_preview <- renderText({
    name <- active_project_display_name()
    if (is.null(name) || !length(name) || !nzchar(name)) "Active project: none" else paste0("Active project: ", name, " (session copy)")
  })

  load_project_into_session <- function(path, display_name = basename(path), activate_path = FALSE) {
    metadata_draft_tracking_enabled(FALSE)
    loaded <- load_proteomics_project_file(path)
    active_project_display_name(display_name)
    ignored_project_upload_paths(current_project_upload_paths())
    active_metadata_override(NULL)
    metadata_replacement_message("No replacement metadata applied this session.")
    validate(need(is.list(loaded) && !is.null(loaded$manifest), "The selected file is not a valid Proteomics Data Workup project."))
    project_db_cache(loaded)
    imported_project_files(list())
    spqc_metadata_draft_edits(empty_spqc_metadata_edits())
    spqc_metadata_edits(if (is.null(loaded$spqc_metadata_edits)) empty_spqc_metadata_edits() else loaded$spqc_metadata_edits)
    settings <- if (is.null(loaded$settings)) list() else loaded$settings
    initialize_applied_metadata(loaded$metadata, settings, spqc_metadata_edits())
    protein_stats_refresh_revision(0L)
    protein_stats_paused(FALSE)
    restored_project_settings(settings)
    project_restore_token(isolate(project_restore_token()) + 1L)
    restore_project_settings(settings)
    session$onFlushed(function() isolate(sync_restored_project_inputs(settings)), once = TRUE)
    session$onFlushed(function() {
      session$onFlushed(function() {
        metadata_draft_tracking_enabled(TRUE)
        applied_state <- isolate(applied_metadata_state())
        if (!is.null(applied_state$metadata)) metadata_apply_message("All metadata changes applied.")
      }, once = TRUE)
    }, once = TRUE)
    cached_batch_corrected_s3_result(NULL)
    if (!is.null(loaded$batch_corrected_s3)) {
      corrected_table <- as.data.frame(loaded$batch_corrected_s3, stringsAsFactors = FALSE, check.names = FALSE)
      abundance_columns <- grep("_batch_corrected_Protein_group_abundance$", colnames(corrected_table), value = TRUE)
      if (length(abundance_columns)) {
        corrected_abundance <- as.data.frame(lapply(corrected_table[abundance_columns], function(x) suppressWarnings(as.numeric(as.character(x)))), check.names = FALSE)
        sample_names <- sub("_batch_corrected_Protein_group_abundance$", "", abundance_columns)
        colnames(corrected_abundance) <- sample_names
        corrected_log2 <- as.data.frame(lapply(corrected_abundance, function(x) log2(x)), check.names = FALSE)
        feature_column <- resolve_report_feature_col(corrected_table, settings$report_feature_col)
        features <- if (!is.na(feature_column)) protein_feature_labels(corrected_table, feature_column) else paste0("Feature_", seq_len(nrow(corrected_table)))
        cached_batch_corrected_s3_result(list(
          table = corrected_table,
          expression = data.frame(Feature = make.unique(features), corrected_log2, check.names = FALSE),
          corrected_log2 = corrected_log2,
          corrected_abundance = corrected_abundance,
          prepared = list(sample_map = if (is.null(loaded$batch_corrected_sample_map)) data.frame(Sample = sample_names) else loaded$batch_corrected_sample_map),
          combat_mode = NA_integer_, correction_method = "Restored project cache"
        ))
      }
    }
    if (isTRUE(activate_path)) {
      normalized_path <- normalizePath(path, mustWork = FALSE)
      active_project_path(normalized_path)
    }
    restored_counts <- c(
      if (!is.null(loaded$metadata)) paste0(nrow(loaded$metadata), " metadata rows") else NULL,
      if (!is.null(loaded$processed_s2)) paste0(nrow(loaded$processed_s2), " S2 proteins") else NULL,
      if (!is.null(loaded$processed_s3)) paste0(nrow(loaded$processed_s3), " S3 proteins") else NULL
    )
    project_database_message(paste0(
      "Opened ", toupper(proteomics_project_backend_for_path(path)), " project: ", display_name, ".",
      if (length(restored_counts)) paste0(" Restored ", paste(restored_counts, collapse = "; "), ".") else " No cached analysis tables were found."
    ))
    loaded
  }

  save_active_project_state <- function(reason = "manual save", include_derived = TRUE, destination = NULL) {
    path <- if (is.null(destination)) active_project_path() else destination
    if (is.null(path) || !nzchar(path)) path <- active_project_path_from_input()
    existing <- if (file.exists(path)) tryCatch(load_proteomics_project_file(path), error = function(e) list()) else list()
    keep <- function(value, name) {
      if (is.null(value) || (is.data.frame(value) && nrow(value) == 0)) existing[[name]] else value
    }
    metadata <- tryCatch(built_metadata(), error = function(e) NULL)
    settings_to_save <- project_settings()
    applied_state <- applied_metadata_state()
    if (!is.null(applied_state$metadata)) settings_to_save$metadata_export_columns <- applied_state$columns
    s2 <- if (isTRUE(include_derived)) tryCatch(protein_no_impute_table(), error = function(e) NULL) else NULL
    s3 <- if (isTRUE(include_derived)) tryCatch(protein_imputed_table(), error = function(e) NULL) else NULL
    corrected_result <- if (isTRUE(include_derived)) tryCatch(batch_corrected_s3_result(), error = function(e) NULL) else NULL
    corrected <- if (!is.null(corrected_result)) corrected_result$table else NULL
    sample_map <- if (!is.null(s3)) attr(s3, "sample_map") else NULL
    identifications_overview <- tryCatch(
      project_input_table(project_file("identifications_overview_file"), project_db_cache()$identifications_overview),
      error = function(e) NULL
    )
    run_identifications_precursor <- tryCatch(
      project_input_table(project_file("run_identifications_precursor_file"), project_db_cache()$run_identifications_precursor),
      error = function(e) NULL
    )
    run_identifications_protein <- tryCatch(
      project_input_table(project_file("run_identifications_protein_file"), project_db_cache()$run_identifications_protein),
      error = function(e) NULL
    )
    project_name <- project_number_for_downloads()
    if (is.na(project_name)) project_name <- tools::file_path_sans_ext(basename(path))
    success <- tryCatch({
      save_proteomics_project_file(
        path = path, project_name = project_name, settings = settings_to_save,
        metadata = keep(metadata, "metadata"), spqc_metadata_edits = spqc_metadata_edits(),
        source_file_manifest = tryCatch(project_file_status(), error = function(e) existing$source_file_manifest),
        processed_s2 = keep(s2, "processed_s2"), processed_s3 = keep(s3, "processed_s3"),
        processed_sample_map = keep(sample_map, "processed_sample_map"),
        batch_corrected_s3 = keep(corrected, "batch_corrected_s3"),
        batch_corrected_sample_map = if (!is.null(corrected_result) && !is.null(corrected_result$prepared$sample_map)) corrected_result$prepared$sample_map else existing$batch_corrected_sample_map,
        statistics_tables = existing$statistics_tables,
        identifications_overview = keep(identifications_overview, "identifications_overview"),
        run_identifications_precursor = keep(run_identifications_precursor, "run_identifications_precursor"),
        run_identifications_protein = keep(run_identifications_protein, "run_identifications_protein")
      )
      TRUE
    }, error = function(e) {
      active_project_last_error(conditionMessage(e))
      project_database_message(paste0("Project save failed (", reason, "): ", conditionMessage(e)))
      showNotification(conditionMessage(e), type = "error")
      FALSE
    })
    if (success) {
      if (is.null(destination)) active_project_path(normalizePath(path, mustWork = FALSE))
      active_project_last_saved(Sys.time())
      active_project_last_error("")
      project_database_message(paste0("Saved active project (", reason, "): ", normalizePath(path, mustWork = FALSE)))
    }
    success
  }

  autosave_active_project <- function(reason, include_derived = FALSE) {
    path <- active_project_path()
    if (!isTRUE(input$autosave_project) || is.null(path) || !nzchar(path)) return(invisible(FALSE))
    active_project_autosave_pending_reason(reason)
    active_project_autosave_include_derived(isTRUE(include_derived))
    active_project_autosave_requested_at(Sys.time())
    invisible(TRUE)
  }

  autosave_project_signal <- shiny::debounce(reactive(list(
    requested_at = active_project_autosave_requested_at(),
    reason = active_project_autosave_pending_reason(),
    include_derived = active_project_autosave_include_derived()
  )), 3000)

  observeEvent(autosave_project_signal(), {
    signal <- autosave_project_signal()
    if (is.null(signal$requested_at) || !nzchar(signal$reason)) return()
    active_project_autosave_pending_reason("")
    active_project_autosave_include_derived(FALSE)
    save_active_project_state(signal$reason, include_derived = isTRUE(signal$include_derived))
  }, ignoreInit = TRUE)

observeEvent(input$project_open_file, {
  req(input$project_open_file)
  tryCatch({
    ext <- tolower(tools::file_ext(input$project_open_file$name))
    if (!ext %in% c("duckdb", "db", "rds")) stop("Select a DuckDB, DB, or RDS project.", call. = FALSE)
    extension <- paste0(".", ext)
    session_path <- tempfile("proteomics_project_", fileext = extension)
    if (!file.copy(input$project_open_file$datapath, session_path, overwrite = TRUE)) {
      stop("Could not copy the selected project into this app session.", call. = FALSE)
    }
    load_project_into_session(session_path, input$project_open_file$name, TRUE)
  }, error = function(e) showNotification(conditionMessage(e), type = "error", duration = NULL))
})

observeEvent(input$create_new_project, {
  showModal(modalDialog(
    title = "Create new project",
    textInput("new_project_filename", "Project filename", value = "proteomics_project.duckdb"),
    tags$p("The project will be kept in this app session. Use Download copy to save it to your computer."),
    easyClose = TRUE,
    footer = tagList(
      modalButton("Cancel"),
      actionButton("confirm_create_new_project", "Create project", class = "btn-primary")
    )
  ))
})

observeEvent(input$confirm_create_new_project, {
  tryCatch({
    filename <- trimws(as.character(input$new_project_filename)[1L])
    if (is.na(filename) || !nzchar(filename)) stop("Enter a project filename.", call. = FALSE)
    filename <- basename(filename)
    if (!grepl("\\.duckdb$", filename, ignore.case = TRUE)) filename <- paste0(filename, ".duckdb")
    extension <- ".duckdb"
    path <- tempfile("proteomics_project_", fileext = extension)
    removeModal()
      empty_project <- empty_proteomics_project_payload()
      do.call(
        save_proteomics_project_file,
        c(
          list(path = path, project_name = tools::file_path_sans_ext(basename(path))),
          empty_project
        )
      )
      ignored_project_upload_paths(current_project_upload_paths())
      load_project_into_session(path, filename, TRUE)
      active_project_last_saved(Sys.time())
      project_database_message(paste0("Created new session project: ", filename, ". Use Download copy to save it to your computer."))
      showNotification("Created a new empty project. Previous session data were not copied.", type = "message")
  }, error = function(e) {
    active_project_last_error(conditionMessage(e))
    showNotification(conditionMessage(e), type = "error", duration = NULL)
  })
})

observeEvent(input$clear_active_project, {
  showModal(modalDialog(
    title = "Clear the active project?",
    p("This will clear all data, metadata, cached analyses, selections, and the active project from this app session."),
    p(tags$strong("Unapplied or undownloaded changes will be lost.")),
    p("No file on your computer will be deleted or modified."),
    easyClose = TRUE,
    footer = tagList(
      modalButton("Cancel"),
      actionButton("confirm_clear_active_project", "Clear project", class = "btn-danger")
    )
  ))
})

observeEvent(input$confirm_clear_active_project, {
  removeModal()
  session$reload()
})

observeEvent(input$save_active_project, {
  path <- active_project_path()
  if (is.null(path) || !length(path) || !nzchar(path)) {
    showNotification("Open or create a project first.", type = "warning")
    return()
  }
  if (!identical(metadata_apply_message(), "All metadata changes applied.")) {
    showNotification("Saving the last applied metadata; draft metadata changes are not included.", type = "warning")
  }
  save_active_project_state("manual save", include_derived = TRUE)
})

  observeEvent(input$apply_metadata_replacement, {
    req(input$metadata_replacement_file)
    tryCatch({
      current <- draft_metadata()
      replacement <- read_uploaded_table(input$metadata_replacement_file)
      replacement <- validate_proteomics_metadata_replacement(replacement, current$Sample)
      replacement_layers <- prepare_proteomics_metadata_replacement(replacement)
      replacement <- replacement_layers$metadata
      spqc_metadata_draft_edits(replacement_layers$spqc_metadata_edits)
      active_metadata_override(replacement)
      metadata_replacement_message(paste0(
        "Loaded ", nrow(replacement), " metadata rows from ", input$metadata_replacement_file$name,
        " into the draft. Click Apply metadata changes to update subsequent tabs."
      ))
      metadata_apply_message("Unapplied metadata changes.")
      showNotification("Modified metadata loaded into the draft.", type = "message")
    }, error = function(e) {
      metadata_replacement_message(paste0("Metadata replacement failed: ", conditionMessage(e)))
      showNotification(conditionMessage(e), type = "error")
    })
  })

  observeEvent(input$run_batch_correction, {
    autosave_active_project("batch correction", include_derived = TRUE)
  }, ignoreInit = TRUE)

  output$project_status <- renderText({
    cache <- project_db_cache()
    paste(
      project_database_message(),
      paste0("Active project: ", if (nzchar(active_project_display_name())) active_project_display_name() else "none"),
      paste0("Last saved: ", if (is.null(active_project_last_saved())) "not saved this session" else format(active_project_last_saved(), "%Y-%m-%d %H:%M:%S")),
      paste0("Restored metadata rows: ", if (is.null(cache$metadata)) 0 else nrow(cache$metadata)),
      paste0("Restored S2 proteins: ", if (is.null(cache$processed_s2)) 0 else nrow(cache$processed_s2)),
      paste0("Restored S3 proteins: ", if (is.null(cache$processed_s3)) 0 else nrow(cache$processed_s3)),
      if (nzchar(active_project_last_error())) paste0("Last error: ", active_project_last_error()) else NULL,
      sep = "\n"
    )
  })

  output$download_project_duckdb <- downloadHandler(
    filename = function() {
      name <- active_project_display_name()
      if (is.null(name) || !length(name) || !nzchar(name)) {
        paste0("proteomics_project_", format(Sys.Date(), "%m%d%y"), ".duckdb")
      } else {
        paste0(tools::file_path_sans_ext(basename(name)), ".duckdb")
      }
    },
    content = function(file) validate(need(save_active_project_state("download", include_derived = TRUE, destination = file), "Could not create DuckDB project download."))
  )

  condition_setup_template <- reactive({
    if (is.null(project_file("condition_setup_template_file"))) return(NULL)
    read_uploaded_table(project_file("condition_setup_template_file"))
  })

  condition_setup_headers <- reactive({
    template <- condition_setup_template()
    if (!is.null(template)) return(colnames(template))
    c("#", "Reference", "Run Label", "Condition", "Fraction", "Replicate", "Quantity Correction Factor", "Label", "Color", "File Name")
  })

  condition_setup_sample_details <- reactive({
    req(project_file("condition_setup_sample_details_file"))
    read_sample_details_workbook(project_file("condition_setup_sample_details_file"))
  })

  condition_setup_inferred_suffix <- reactive({
    template <- condition_setup_template()
    if (is.null(template) || !"Run Label" %in% colnames(template) || nrow(template) == 0) return("")
    run_label <- as.character(template$`Run Label`[which(!is.na(template$`Run Label`) & nzchar(template$`Run Label`))[1]])
    if (is.na(run_label) || !nzchar(run_label)) return("")
    sub("^ID[0-9]+", "", run_label)
  })

  parse_condition_setup_replicate <- function(values) {
    values <- trimws(as.character(values))
    labels <- values
    replicates <- rep(NA_character_, length(values))
    pattern <- "^(.*?)[ _.-]+(?:rep(?:licate)?|r)?[ _.-]*0*([0-9]+)$"
    matched <- grepl(pattern, values, ignore.case = TRUE)
    if (any(matched)) {
      parsed_labels <- trimws(sub(pattern, "\\1", values[matched], ignore.case = TRUE))
      parsed_labels <- sub("[ _.-]+$", "", parsed_labels)
      parsed_replicates <- sub(pattern, "\\2", values[matched], ignore.case = TRUE)
      usable <- nzchar(parsed_labels) & grepl("[A-Za-z]", parsed_labels) & grepl("^[0-9]+$", parsed_replicates)
      matched_idx <- which(matched)
      labels[matched_idx[usable]] <- parsed_labels[usable]
      replicates[matched_idx[usable]] <- as.character(as.integer(parsed_replicates[usable]))
    }
    data.frame(Label = labels, Replicate = replicates, stringsAsFactors = FALSE)
  }

  infer_condition_setup_label_col <- function(details) {
    choices <- setdiff(colnames(details), c("SampleDetailsID", "SampleName"))
    if (!length(choices)) return(NA_character_)

    score_column <- function(column_name) {
      values <- trimws(as.character(details[[column_name]]))
      values <- values[!is.na(values) & nzchar(values)]
      if (!length(values)) return(-Inf)

      unique_n <- length(unique(values))
      repeated_n <- sum(table(values) > 1)
      numeric_ratio <- mean(!is.na(suppressWarnings(as.numeric(values))))
      lower_name <- tolower(column_name)
      parsed <- parse_condition_setup_replicate(values)
      parsed_detected <- !is.na(parsed$Replicate) & nzchar(parsed$Label)
      if (lower_name %in% c("name", "sample name") && sum(parsed_detected) >= 2) {
        parsed_labels <- parsed$Label[parsed_detected]
        parsed_unique_n <- length(unique(parsed_labels))
        parsed_repeated_n <- sum(table(parsed_labels) > 1)
        if (parsed_unique_n >= 1 && parsed_unique_n <= max(20, ceiling(nrow(details) * 0.6))) {
          return(80 + min(10, parsed_repeated_n) + max(0, 10 - parsed_unique_n))
        }
      }
      excluded_names <- c("id", "name", "solvent", "solvent/buffer", "volume", "ph", "quantity")
      if (lower_name %in% excluded_names) return(-Inf)
      if (grepl("id|sample|volume|quantity|amount|mass|ph|solvent", lower_name)) return(-Inf)
      if (unique_n < 2 || unique_n > max(20, ceiling(nrow(details) * 0.6))) return(-Inf)
      if (repeated_n < 1) return(-Inf)

      name_bonus <- if (grepl("condition|group|label|class|type|status|disease|notes", lower_name)) 30 else 0
      notes_bonus <- if (identical(lower_name, "notes")) 20 else 0
      text_bonus <- if (numeric_ratio < 0.2) 10 else 0
      repeat_bonus <- min(10, repeated_n)
      compact_bonus <- max(0, 10 - unique_n)
      name_bonus + notes_bonus + text_bonus + repeat_bonus + compact_bonus
    }

    scores <- vapply(choices, score_column, numeric(1))
    if (all(!is.finite(scores))) {
      if ("Notes" %in% choices) return("Notes")
      return(choices[1])
    }
    names(which.max(scores))
  }

  observeEvent(condition_setup_template(), {
    suffix <- condition_setup_inferred_suffix()
    if (nzchar(suffix) && !nzchar(isolate(input$condition_setup_run_label_suffix))) {
      updateTextInput(session, "condition_setup_run_label_suffix", value = suffix)
    }
  }, ignoreInit = TRUE)

  observeEvent(condition_setup_sample_details(), {
    details <- condition_setup_sample_details()
    choices <- setdiff(colnames(details), c("SampleDetailsID", "SampleName"))
    inferred_condition <- infer_condition_setup_label_col(details)
    condition_choices <- c("Infer from sample details" = "__infer__", stats::setNames(choices, choices))
    selected_condition <- if (!is.na(inferred_condition)) "__infer__" else choices[1]
    selected_replicate <- if ("ID" %in% choices) "ID" else choices[1]
    replicate_choices <- stats::setNames(choices, choices)
    if ("ID" %in% choices) names(replicate_choices)[replicate_choices == "ID"] <- "Sample ID order"
    updateSelectInput(session, "condition_setup_condition_col", choices = condition_choices, selected = selected_condition)
    updateSelectInput(session, "condition_setup_replicate_order_col", choices = replicate_choices, selected = selected_replicate)
  }, ignoreInit = FALSE)

  condition_setup_table <- reactive({
    details <- condition_setup_sample_details()
    condition_col <- input$condition_setup_condition_col
    inferred_condition_col <- infer_condition_setup_label_col(details)
    if (is.null(condition_col) || identical(condition_col, "__infer__") || !condition_col %in% colnames(details)) {
      condition_col <- if (!is.na(inferred_condition_col)) inferred_condition_col else colnames(details)[1]
    }
    replicate_order_col <- input$condition_setup_replicate_order_col
    if (is.null(replicate_order_col) || !replicate_order_col %in% colnames(details)) {
      replicate_order_col <- if ("ID" %in% colnames(details)) "ID" else condition_col
    }
    suffix <- input$condition_setup_run_label_suffix
    if (is.null(suffix) || !nzchar(suffix)) suffix <- condition_setup_inferred_suffix()
    run_label <- paste0(details$SampleDetailsID, suffix)
    parsed_conditions <- parse_condition_setup_replicate(details[[condition_col]])
    parsed_condition_has_replicate <- !is.na(parsed_conditions$Replicate) & nzchar(parsed_conditions$Label)
    condition_values <- ifelse(parsed_condition_has_replicate, parsed_conditions$Label, as.character(details[[condition_col]]))
    condition_values[is.na(condition_values) | !nzchar(condition_values)] <- "Not Defined"

    order_values <- details[[replicate_order_col]]
    order_key <- suppressWarnings(as.numeric(as.character(order_values)))
    if (all(is.na(order_key))) order_key <- tolower(as.character(order_values))
    replicate_values <- integer(nrow(details))
    for (condition in unique(condition_values)) {
      idx <- which(condition_values == condition)
      replicate_values[idx[order(order_key[idx], seq_along(idx), na.last = TRUE)]] <- seq_along(idx)
    }
    replicate_values[parsed_condition_has_replicate] <- as.integer(parsed_conditions$Replicate[parsed_condition_has_replicate])

    condition_levels <- unique(condition_values)
    known_colors <- c(
      "CKD" = "#00A651",
      "CKDu" = "#597DFF",
      "Control" = "#ED3438",
      "SPQC" = "#FF00CB",
      "Not Defined" = "#505050"
    )
    extra_levels <- setdiff(condition_levels, names(known_colors))
    if (length(extra_levels)) {
      known_colors <- c(known_colors, stats::setNames(scales::hue_pal()(length(extra_levels)), extra_levels))
    }
    color_values <- unname(known_colors[condition_values])

    out <- data.frame(
      check.names = FALSE,
      "#" = seq_len(nrow(details)),
      "Reference" = "False",
      "Run Label" = run_label,
      "Condition" = condition_values,
      "Fraction" = "NA",
      "Replicate" = replicate_values,
      "Quantity Correction Factor" = 1,
      "Label" = condition_values,
      "Color" = color_values,
      "File Name" = sub("\\.[Hh][Tt][Rr][Mm][Ss]$", "", run_label)
    )
    headers <- condition_setup_headers()
    missing_headers <- setdiff(headers, colnames(out))
    for (header in missing_headers) out[[header]] <- NA_character_
    out <- out[, headers, drop = FALSE]

    template <- condition_setup_template()
    if (!is.null(template) && "Run Label" %in% colnames(template)) {
      template_missing_headers <- setdiff(headers, colnames(template))
      for (header in template_missing_headers) template[[header]] <- NA_character_
      template <- template[, headers, drop = FALSE]
      template_run_labels <- as.character(template$`Run Label`)
      generated_run_labels <- as.character(out$`Run Label`)
      template_only <- template[!(template_run_labels %in% generated_run_labels), , drop = FALSE]
      if (nrow(template_only) > 0) {
        if ("Run Label" %in% colnames(template_only)) {
          qc_rows <- grepl("^IDSPQC_", as.character(template_only$`Run Label`), ignore.case = TRUE)
          if (any(qc_rows)) {
            if ("Condition" %in% colnames(template_only)) template_only$Condition[qc_rows] <- "SPQC"
            if ("Label" %in% colnames(template_only)) template_only$Label[qc_rows] <- "SPQC"
            if ("Color" %in% colnames(template_only)) template_only$Color[qc_rows] <- "#FF00CB"
            if ("Replicate" %in% colnames(template_only)) {
              qc_replicates <- sub("^IDSPQC_0*([0-9]+).*", "\\1", as.character(template_only$`Run Label`[qc_rows]), ignore.case = TRUE)
              qc_replicates[!grepl("^[0-9]+$", qc_replicates)] <- as.character(seq_len(sum(qc_rows)))[!grepl("^[0-9]+$", qc_replicates)]
              template_only$Replicate[qc_rows] <- qc_replicates
            }
          }
        }
        out[] <- lapply(out, as.character)
        template_only[] <- lapply(template_only, as.character)
        out <- dplyr::bind_rows(out, template_only)
        if ("#" %in% colnames(out)) {
          row_order <- suppressWarnings(as.numeric(as.character(out[["#"]])))
          out <- out[order(is.na(row_order), row_order, seq_len(nrow(out))), , drop = FALSE]
        }
      }
    }

    out
  })

  evosep_manifest_details <- reactive({
    req(project_file("evosep_manifest_file"))
    details <- read_sample_details_workbook(project_file("evosep_manifest_file"))
    details[, c("SampleDetailsID", "SampleName", setdiff(colnames(details), c("SampleDetailsID", "SampleName"))), drop = FALSE]
  })

  evosep_template_text <- reactive({
    req(project_file("evosep_template_file"))
    paste(readLines(project_file("evosep_template_file")$datapath, warn = FALSE), collapse = "\n")
  })

  observeEvent(evosep_template_text(), {
    defaults <- tryCatch(evosep_template_defaults(evosep_template_text()), error = function(error) NULL)
    if (is.null(defaults)) return()
    updateTextInput(session, "evosep_source_tray", value = defaults$source_tray)
    updateTextInput(session, "evosep_output_dir", value = defaults$output_dir)
    updateTextInput(session, "evosep_filename_suffix", value = defaults$sample_filename_suffix)
    updateTextInput(session, "evosep_comment", value = defaults$comment)
    updateTextAreaInput(session, "evosep_analysis_method", value = defaults$analysis_method)
    updateTextAreaInput(session, "evosep_xcalibur_method", value = defaults$xcalibur_method)
  }, ignoreInit = TRUE)

  evosep_assignment_state <- reactiveVal(NULL)

  evosep_default_assignment_table <- reactive({
    grid <- evosep_plate_grid()
    grid$Assignment <- ""
    grid$ManifestID <- ""
    grid$SampleName <- ""
    grid$Filename <- ""

    suffix <- input$evosep_filename_suffix
    if (is.null(suffix)) suffix <- ""

    adh_positions <- parse_position_list(input$evosep_adh_positions)
    spqc_positions <- parse_position_list(input$evosep_spqc_positions)
    sample_positions <- parse_position_list(input$evosep_sample_positions)

    if (length(adh_positions)) {
      idx <- match(adh_positions, grid$Position)
      grid$Assignment[idx] <- "ADH"
      grid$SampleName[idx] <- "ADH"
      grid$Filename[idx] <- paste0("ADH", vapply(seq_along(idx), function(i) evosep_incremented_suffix(suffix, i), character(1)))
    }

    if (length(spqc_positions)) {
      idx <- match(spqc_positions, grid$Position)
      spqc_prefix <- input$evosep_spqc_prefix
      if (is.null(spqc_prefix) || !nzchar(spqc_prefix)) spqc_prefix <- "IDSPQC"
      grid$Assignment[idx] <- "SPQC"
      grid$SampleName[idx] <- "SPQC"
      grid$Filename[idx] <- paste0(spqc_prefix, vapply(seq_along(idx), function(i) evosep_incremented_suffix(suffix, i), character(1)))
    }

    if (isTRUE(input$evosep_fill_manifest_order) && !is.null(project_file("evosep_manifest_file"))) {
      manifest <- evosep_manifest_details()
      sample_positions <- setdiff(sample_positions, c(adh_positions, spqc_positions))
      usable_n <- min(nrow(manifest), length(sample_positions))
      if (usable_n > 0) {
        idx <- match(sample_positions[seq_len(usable_n)], grid$Position)
        grid$Assignment[idx] <- "Sample"
        grid$ManifestID[idx] <- manifest$SampleDetailsID[seq_len(usable_n)]
        grid$SampleName[idx] <- manifest$SampleName[seq_len(usable_n)]
        grid$Filename[idx] <- paste0(manifest$SampleDetailsID[seq_len(usable_n)], suffix)
      }
    }

    manifest <- if (!is.null(project_file("evosep_manifest_file"))) evosep_manifest_details() else NULL
    evosep_normalize_plate_table(
      grid,
      manifest = manifest,
      suffix = suffix,
      spqc_prefix = input$evosep_spqc_prefix
    )
  })

  observeEvent(input$evosep_reset_assignment_table, {
    evosep_assignment_state(evosep_default_assignment_table())
  }, ignoreInit = TRUE)

  observeEvent(input$evosep_assignment_table_cell_edit, {
    current <- evosep_assignment_state()
    if (is.null(current)) current <- evosep_default_assignment_table()
    edited <- DT::editData(current, input$evosep_assignment_table_cell_edit, rownames = FALSE)
    manifest <- if (!is.null(project_file("evosep_manifest_file"))) evosep_manifest_details() else NULL
    evosep_assignment_state(evosep_normalize_plate_table(
      edited,
      manifest = manifest,
      suffix = input$evosep_filename_suffix,
      spqc_prefix = input$evosep_spqc_prefix
    ))
  })

  observeEvent(project_file("evosep_manifest_file"), {
    evosep_assignment_state(NULL)
  }, ignoreInit = TRUE)

  observeEvent(project_file("evosep_template_file"), {
    evosep_assignment_state(NULL)
  }, ignoreInit = TRUE)

  evosep_assigned_plate <- reactive({
    state <- evosep_assignment_state()
    if (!is.null(state)) {
      manifest <- if (!is.null(project_file("evosep_manifest_file"))) evosep_manifest_details() else NULL
      return(evosep_normalize_plate_table(
        state,
        manifest = manifest,
        suffix = input$evosep_filename_suffix,
        spqc_prefix = input$evosep_spqc_prefix
      ))
    }
    evosep_default_assignment_table()
  })

  evosep_queue_table <- reactive({
    plate <- evosep_assigned_plate()
    assigned <- plate[plate$Assignment != "Empty", , drop = FALSE]
    assigned <- assigned[order(assigned$Position), , drop = FALSE]
    validate(need(nrow(assigned) > 0, "Assign at least one ADH, SPQC, or manifest sample position."))
    data.frame(
      RunNumber = seq_len(nrow(assigned)),
      AnalysisMethod = input$evosep_analysis_method,
      SourceTray = input$evosep_source_tray,
      SourceVial = assigned$Position,
      SampleName = ifelse(assigned$Assignment %in% c("ADH", "SPQC"), assigned$SampleName, ""),
      XcaliburMethod = input$evosep_xcalibur_method,
      XcaliburFilename = assigned$Filename,
      PostAcquisitionProgram = "",
      OutputDir = input$evosep_output_dir,
      Comment = input$evosep_comment,
      Assignment = assigned$Assignment,
      ManifestID = assigned$ManifestID,
      ManifestSampleName = assigned$SampleName,
      PlateRow = assigned$Row,
      PlateColumn = assigned$Column,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  raw_data <- reactive({
    pca_source <- input$pca_data_source
    if (is.null(pca_source) || length(pca_source) != 1 || !nzchar(pca_source)) pca_source <- "imputed"
    switch(
      pca_source,
      "imputed" = {
        protein_source_table("imputed")
      },
      "no_impute" = {
        protein_source_table("no_impute")
      },
      "custom" = {
        req(project_file("data_file"))
        read_uploaded_table(project_file("data_file"))
      }
    )
  })

  report_quantity_columns <- reactive({
    proteomics_abundance_columns(raw_data())
  })

  use_report_layout <- reactive({
    pca_source <- input$pca_data_source
    if (is.null(pca_source) || length(pca_source) != 1 || !nzchar(pca_source)) pca_source <- "imputed"
    data_layout <- input$data_layout
    if (is.null(data_layout) || length(data_layout) != 1 || !nzchar(data_layout)) data_layout <- "Auto detect"
    if (pca_source %in% c("imputed", "no_impute")) return(TRUE)
    if (data_layout == "DIA-NN protein report (.PG.Quantity columns)") return(TRUE)
    if (data_layout == "Feature-by-sample matrix") return(FALSE)
    length(report_quantity_columns()) > 0
  })

  detected_annotation_count <- reactive({
    if (use_report_layout()) return(0)
    df <- raw_data()
    first_col <- trimws(as.character(df[[1]]))
    known <- c("group", "subject", "timepoint", "batch", "class")
    auto_n <- sum(tolower(first_col[seq_len(min(10, length(first_col)))]) %in% known)
    if (isTRUE(input$auto_detect_annotations) && auto_n > 0) auto_n else input$annotation_rows
  })

  detected_annotations <- reactive({
    if (use_report_layout()) return(NULL)
    df <- raw_data()
    ann_n <- min(detected_annotation_count(), nrow(df))
    if (ann_n < 1) return(NULL)

    ann_df <- df[seq_len(ann_n), , drop = FALSE]
    sample_names <- colnames(df)[-1]
    out <- data.frame(Sample = sample_names, stringsAsFactors = FALSE)

    for (i in seq_len(ann_n)) {
      ann_name <- as.character(ann_df[i, 1])
      if (is.na(ann_name) || ann_name == "") ann_name <- paste0("Annotation_", i)
      out[[ann_name]] <- as.character(unlist(ann_df[i, -1, drop = TRUE]))
    }

    out
  })

  resolve_report_feature_col <- function(df, requested_col = NULL) {
    requested_col <- if (is.null(requested_col) || length(requested_col) == 0) "" else as.character(requested_col)[1]
    if (nzchar(requested_col) && requested_col %in% colnames(df)) return(requested_col)
    candidate_cols <- c("PG.Genes", "PG.ProteinGroups", "PG.ProteinNames", colnames(df)[1])
    candidate <- candidate_cols[candidate_cols %in% colnames(df)][1]
    if (is.na(candidate)) NA_character_ else candidate
  }

  metadata_for_exclusion_matching <- function() {
    if (!is.null(project_db_cache()$metadata) || !is.null(project_file("meta_file")) || !is.null(project_file("run_order_file")) || !is.null(project_file("sample_details_file"))) {
      built_metadata()
    } else {
      NULL
    }
  }

  protein_sample_names_for_table <- function(table, abundance_columns) {
    md <- tryCatch(built_metadata(), error = function(e) NULL)
    if (is.null(md) || !"Sample" %in% colnames(md)) return(proteomics_abundance_sample_names(abundance_columns))
    header_labels <- protein_header_labels_from_metadata(md, input$protein_header_label_columns)
    labels <- proteomics_abundance_sample_names(abundance_columns)
    if (!any(endsWith(abundance_columns, "_Protein_group_abundance"))) return(labels)
    resolved <- resolve_proteomics_processed_sample_ids(
      labels,
      md,
      header_labels,
      project_db_cache()$processed_sample_map
    )
    labels[!is.na(resolved)] <- resolved[!is.na(resolved)]
    labels
  }

  expression_data <- reactive({
    df <- raw_data()

    if (use_report_layout()) {
      quantity_cols <- report_quantity_columns()
      validate(
        need(length(quantity_cols) >= 3, "DIA-NN report layout requires at least 3 .PG.Quantity sample columns.")
      )

      feature_col <- resolve_report_feature_col(df, input$report_feature_col)
      validate(need(!is.na(feature_col), "No feature identifier column found in the report."))

      mat_df <- df[, quantity_cols, drop = FALSE]
      mat_df[] <- lapply(mat_df, function(x) suppressWarnings(as.numeric(as.character(x))))
      sample_names <- protein_sample_names_for_table(df, quantity_cols)
      colnames(mat_df) <- make.unique(sample_names)
      mat_df <- filter_expression_columns_by_exclusion(
        data.frame(Feature = seq_len(nrow(mat_df)), mat_df, check.names = FALSE),
        input$sample_exclusions_text,
        metadata_for_exclusion_matching()
      )[, -1, drop = FALSE]

      keep_rows <- rowSums(!is.na(mat_df)) > 0
      feature_names <- protein_feature_labels(df, feature_col)
      return(data.frame(
        Feature = make.unique(feature_names[keep_rows]),
        mat_df[keep_rows, , drop = FALSE],
        check.names = FALSE
      ))
    }

    ann_n <- min(detected_annotation_count(), nrow(df))

    validate(
      need(ncol(df) >= 3, "Data file must have at least 3 columns: feature column + sample columns."),
      need(nrow(df) > ann_n + 1, "Not enough non-annotation rows for plotting.")
    )

    feature_col <- df[(ann_n + 1):nrow(df), 1, drop = TRUE]
    mat_df <- df[(ann_n + 1):nrow(df), -1, drop = FALSE]
    mat_df[] <- lapply(mat_df, function(x) suppressWarnings(as.numeric(as.character(x))))
    mat_df <- filter_expression_columns_by_exclusion(
      data.frame(Feature = seq_len(nrow(mat_df)), mat_df, check.names = FALSE),
      input$sample_exclusions_text
    )[, -1, drop = FALSE]

    keep_rows <- rowSums(!is.na(mat_df)) > 0
    feature_col <- feature_col[keep_rows]
    mat_df <- mat_df[keep_rows, , drop = FALSE]

    data.frame(Feature = make.unique(as.character(feature_col)), mat_df, check.names = FALSE)
  })

  expression_data_for_protein_source <- function(source) {
    if (identical(source, "S3_batch_corrected")) {
      return(filter_expression_columns_by_exclusion(batch_corrected_s3_expression_data(), input$sample_exclusions_text, metadata_for_exclusion_matching()))
    }
    df <- protein_source_table(source)
    quantity_cols <- proteomics_abundance_columns(df)
    validate(need(length(quantity_cols) >= 3, "Protein report requires at least 3 abundance sample columns."))

    feature_col <- resolve_report_feature_col(df, input$report_feature_col)
    validate(need(!is.na(feature_col), "No feature identifier column found in the report."))

    mat_df <- df[, quantity_cols, drop = FALSE]
    mat_df[] <- lapply(mat_df, function(x) suppressWarnings(as.numeric(as.character(x))))
    sample_names <- protein_sample_names_for_table(df, quantity_cols)
    colnames(mat_df) <- make.unique(sample_names)
    mat_df <- filter_expression_columns_by_exclusion(
      data.frame(Feature = seq_len(nrow(mat_df)), mat_df, check.names = FALSE),
      input$sample_exclusions_text,
      metadata_for_exclusion_matching()
    )[, -1, drop = FALSE]

    keep_rows <- rowSums(!is.na(mat_df)) > 0
    feature_names <- protein_feature_labels(df, feature_col)
    data.frame(
      Feature = make.unique(feature_names[keep_rows]),
      mat_df[keep_rows, , drop = FALSE],
      check.names = FALSE
    )
  }

  protein_info_for_source <- function(source, feature_col_requested = NULL) {
    if (identical(source, "S3_batch_corrected")) {
      result <- batch_corrected_s3_result()
      info <- result$prepared$feature_info
      out <- data.frame(
        Protein = info$Feature,
        RawFeatureID = info$Feature,
        SourceRowIndex = seq_len(nrow(info)),
        KeptRowIndex = seq_len(nrow(info)),
        stringsAsFactors = FALSE
      )
      for (column_name in intersect(c("PG.Genes", "PG.ProteinGroups", "PG.ProteinNames", "PG.ProteinDescriptions"), colnames(info))) {
        out[[column_name]] <- as.character(info[[column_name]])
      }
      return(out)
    }
    df <- protein_source_table(source)
    feature_col <- resolve_report_feature_col(df, feature_col_requested)
    validate(need(!is.na(feature_col), "No feature identifier column found in the protein report."))
    quantity_cols <- proteomics_abundance_columns(df)
    mat_df <- df[, quantity_cols, drop = FALSE]
    mat_df[] <- lapply(mat_df, function(x) suppressWarnings(as.numeric(as.character(x))))
    keep_rows <- rowSums(!is.na(mat_df)) > 0
    feature_names <- protein_feature_labels(df, feature_col)
    info_cols <- intersect(c("PG.Genes", "PG.ProteinGroups", "PG.ProteinNames", "PG.ProteinDescriptions"), colnames(df))
    out <- data.frame(
      Protein = make.unique(feature_names[keep_rows]),
      RawFeatureID = trimws(as.character(df[[1]][keep_rows])),
      SourceRowIndex = which(keep_rows),
      KeptRowIndex = seq_len(sum(keep_rows)),
      stringsAsFactors = FALSE
    )
    for (column_name in info_cols) out[[column_name]] <- as.character(df[[column_name]][keep_rows])
    out
  }

  compact_choice_text <- function(values, max_chars = 90) {
    values <- unique(trimws(as.character(values)))
    values <- values[!is.na(values) & nzchar(values)]
    if (!length(values)) return(character(0))
    text <- paste(values, collapse = "; ")
    if (nchar(text) > max_chars) paste0(substr(text, 1, max_chars - 3), "...") else text
  }

  protein_feature_choices <- function(source, expr = NULL) {
    features <- if (!is.null(expr) && "Feature" %in% colnames(expr)) {
      as.character(expr$Feature)
    } else {
      as.character(protein_info_for_source(source, input$report_feature_col)$Protein)
    }
    features <- features[!is.na(features) & nzchar(features)]
    if (!length(features)) return(stats::setNames(character(0), character(0)))

    info <- tryCatch(protein_info_for_source(source, input$report_feature_col), error = function(e) NULL)
    if (is.null(info) || !"Protein" %in% colnames(info)) {
      return(stats::setNames(features, features))
    }
    info <- info[match(features, as.character(info$Protein)), , drop = FALSE]
    label_columns <- intersect(
      c("PG.Genes", "PG.ProteinGroups", "PG.ProteinNames", "PG.ProteinDescriptions", "RawFeatureID"),
      colnames(info)
    )
    labels <- vapply(seq_along(features), function(index) {
      parts <- c(features[index])
      for (column_name in label_columns) {
        value <- compact_choice_text(info[[column_name]][index])
        if (length(value) && !identical(value, features[index])) parts <- c(parts, value)
      }
      paste(unique(parts), collapse = " | ")
    }, character(1))
    stats::setNames(features, labels)
  }

  sample_metadata_for_plotting <- function(samples) {
    md <- if (!is.null(project_db_cache()$metadata) || !is.null(project_file("meta_file")) || !is.null(project_file("run_order_file")) || !is.null(project_file("sample_details_file"))) {
      active_metadata()
    } else {
      data.frame(Sample = samples, stringsAsFactors = FALSE)
    }
    out <- unique_sample_metadata(md, samples)
    if (!"Condition" %in% colnames(out)) out$Condition <- "All"
    if (!"Replicate" %in% colnames(out)) out$Replicate <- seq_len(nrow(out))
    if (!"AnalysisLabel" %in% colnames(out)) out$AnalysisLabel <- out$Sample
    out
  }

  active_metadata <- reactive({
    md <- normalize_proteomics_metadata(built_metadata())
    if (!"Excluded" %in% colnames(md)) return(md)
    excluded <- as.logical(md$Excluded)
    excluded[is.na(excluded)] <- FALSE
    md[!excluded, , drop = FALSE]
  })

  parsed_data <- reactive({
    expr <- expression_data()
    mat <- as.matrix(expr[, -1, drop = FALSE])
    rownames(mat) <- expr$Feature
    X <- t(mat)
    keep <- colSums(!is.na(X)) > 0
    X <- X[, keep, drop = FALSE]

    validate(
      need(nrow(X) >= 3, "Need at least 3 samples for PCA."),
      need(ncol(X) >= 2, "Need at least 2 non-empty features after removing all-NA features.")
    )

    data.frame(Sample = rownames(X), X, check.names = FALSE)
  })

  metadata_df <- reactive({
    req(project_file("meta_file"))
    read_uploaded_table(project_file("meta_file"))
  })

  run_order_df <- reactive({
    req(project_file("run_order_file"))
    read_uploaded_table(project_file("run_order_file"))
  })

  sample_details_df <- reactive({
    req(project_file("sample_details_file"))
    details <- read_sample_details_workbook(project_file("sample_details_file"))
    details[, c("SampleDetailsID", "SampleName", setdiff(colnames(details), c("SampleDetailsID", "SampleName"))), drop = FALSE]
  })

  resolve_metadata_sample_id_col <- function(md, requested_col = NULL) {
    requested_col <- if (is.null(requested_col) || length(requested_col) == 0) "" else as.character(requested_col)[1]
    if (nzchar(requested_col) && requested_col %in% colnames(md)) return(requested_col)
    fallback_cols <- c("Run Label", "Sample", "#")
    fallback <- fallback_cols[fallback_cols %in% colnames(md)][1]
    if (!is.na(fallback) && nzchar(fallback)) return(fallback)
    colnames(md)[1]
  }

  draft_metadata <- reactive({
    sample_details_used_as_basis <- FALSE
    replacement_override <- active_metadata_override()
    cache <- project_db_cache()
    metadata_sources <- c("meta_file", "run_order_file", "sample_details_file")
    uploaded_metadata <- vapply(metadata_sources, function(id) !is.null(project_file(id)), logical(1))
    used_cached_overlay <- !is.null(cache$metadata) || !is.null(replacement_override)
    if (!is.null(replacement_override)) {
      built <- as.data.frame(replacement_override, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (used_cached_overlay) {
      built <- as.data.frame(cache$metadata, stringsAsFactors = FALSE, check.names = FALSE)
      if (uploaded_metadata[["meta_file"]]) built <- merge_proteomics_metadata_replacement(built, metadata_df(), "meta_file")
      if (uploaded_metadata[["run_order_file"]]) built <- merge_proteomics_metadata_replacement(built, run_order_df(), "run_order_file")
      if (uploaded_metadata[["sample_details_file"]]) built <- merge_proteomics_metadata_replacement(built, sample_details_df(), "sample_details_file")
    } else if (!is.null(project_file("meta_file"))) {
      md <- metadata_df()
      id_col <- resolve_metadata_sample_id_col(md, input$sample_id_col)
      validate(need(id_col %in% colnames(md), paste0("Metadata sample ID column not found: ", id_col)))
      validate(need(!anyDuplicated(md[[id_col]]), paste0("Metadata sample ID column contains duplicate values: ", id_col)))
      built <- md
      built$Sample <- as.character(built[[id_col]])
    } else if (!is.null(project_file("sample_details_file"))) {
      built <- sample_details_df()
      built$Sample <- as.character(built$SampleDetailsID)
      sample_details_used_as_basis <- TRUE
    } else {
      built <- data.frame(Sample = parsed_data()$Sample, stringsAsFactors = FALSE)
    }
    detected_sample_details_id <- stringr::str_extract(built$Sample, "ID[0-9]+")
    if ("SampleDetailsID" %in% colnames(built)) {
      built$SampleDetailsID <- fill_missing_metadata_values(built$SampleDetailsID, detected_sample_details_id)
    } else {
      built$SampleDetailsID <- detected_sample_details_id
    }

    if (!used_cached_overlay && !is.null(project_file("run_order_file"))) {
      run_order <- run_order_df()
      validate(need(all(c("Run Label", "#") %in% colnames(run_order)), "Run-order table must contain 'Run Label' and '#' columns."))
      validate(need(!anyDuplicated(run_order$`Run Label`), "Run-order table contains duplicate Run Label values."))
      run_order$Sample <- as.character(run_order$`Run Label`)
      run_order$RunOrder <- run_order$`#`
      run_order <- run_order[, c("Sample", "RunOrder", setdiff(colnames(run_order), c("Sample", "RunOrder", "Run Label", "#"))), drop = FALSE]
      built <- built %>% left_join(run_order, by = "Sample")
      built <- built %>% arrange(is.na(RunOrder), RunOrder)

      batch_col <- preferred_metadata_column(colnames(built), c("Batch", "batch"))
      if (!is.na(batch_col)) {
        built$Batch <- fill_missing_metadata_values(if ("Batch" %in% colnames(built)) built$Batch else NULL, built[[batch_col]])
      }

      group_col <- preferred_metadata_column(colnames(built), c("Group", "group", "Condition", "condition"))
      if (!is.na(group_col)) {
        built$Group <- fill_missing_metadata_values(if ("Group" %in% colnames(built)) built$Group else NULL, built[[group_col]])
        built$Condition <- fill_missing_metadata_values(if ("Condition" %in% colnames(built)) built$Condition else NULL, built[[group_col]])
      }
    }

    if (!used_cached_overlay && !is.null(project_file("sample_details_file")) && !isTRUE(sample_details_used_as_basis)) {
      built <- built %>% left_join(sample_details_df(), by = "SampleDetailsID")
      for (column in colnames(built)[grepl("\\.x$", colnames(built))]) {
        base_column <- sub("\\.x$", "", column)
        paired_column <- paste0(base_column, ".y")
        if (paired_column %in% colnames(built)) {
          built[[base_column]] <- fill_missing_metadata_values(built[[column]], built[[paired_column]])
          built[[column]] <- NULL
          built[[paired_column]] <- NULL
        }
      }
      for (column in colnames(built)[grepl("\\.y$", colnames(built))]) {
        base_column <- sub("\\.y$", "", column)
        if (!base_column %in% colnames(built)) {
          names(built)[names(built) == column] <- base_column
        }
      }
      batch_col <- preferred_metadata_column(colnames(built), c("Batch", "batch"))
      if (!is.na(batch_col)) {
        built$Batch <- fill_missing_metadata_values(if ("Batch" %in% colnames(built)) built$Batch else NULL, built[[batch_col]])
      }
      group_col <- preferred_metadata_column(colnames(built), c("Group", "group", "Condition", "condition"))
      if (!is.na(group_col)) {
        built$Group <- fill_missing_metadata_values(if ("Group" %in% colnames(built)) built$Group else NULL, built[[group_col]])
        built$Condition <- fill_missing_metadata_values(if ("Condition" %in% colnames(built)) built$Condition else NULL, built[[group_col]])
      }
    } else if (!used_cached_overlay) {
      built$SampleName <- NA_character_
    }

    batch_col <- preferred_metadata_column(colnames(built), c("Batch", "batch"))
    if (!is.na(batch_col)) {
      built$Batch <- fill_missing_metadata_values(if ("Batch" %in% colnames(built)) built$Batch else NULL, built[[batch_col]])
    }
    group_col <- preferred_metadata_column(colnames(built), c("Group", "group", "Condition", "condition"))
    if (!is.na(group_col)) {
      built$Group <- fill_missing_metadata_values(if ("Group" %in% colnames(built)) built$Group else NULL, built[[group_col]])
      built$Condition <- fill_missing_metadata_values(if ("Condition" %in% colnames(built)) built$Condition else NULL, built[[group_col]])
    }

    spqc_rows <- grepl("SPQC", built$Sample, ignore.case = TRUE) |
      ("SampleName" %in% colnames(built) & grepl("SPQC", built$SampleName, ignore.case = TRUE))
    spqc_rows[is.na(spqc_rows)] <- FALSE
    if (any(spqc_rows)) {
      built$Group <- fill_missing_metadata_values(if ("Group" %in% colnames(built)) built$Group else NULL, rep("", nrow(built)))
      built$Condition <- fill_missing_metadata_values(if ("Condition" %in% colnames(built)) built$Condition else NULL, rep("", nrow(built)))
      built$Batch <- infer_spqc_batches_from_run_dates(built$Sample, if ("Batch" %in% colnames(built)) built$Batch else rep("", nrow(built)), spqc_rows)
      built$Batch <- apply_spqc_batch_overrides(
        samples = built$Sample,
        sample_names = if ("SampleName" %in% colnames(built)) built$SampleName else rep("", nrow(built)),
        batches = built$Batch,
        spqc_rows = spqc_rows,
        override_text = input$spqc_batch_overrides
      )

      spqc_mode <- input$spqc_assignment_mode
      if (is.null(spqc_mode) || !spqc_mode %in% c("single", "batch", "date", "keep")) spqc_mode <- "single"
      spqc_label <- trimws(as.character(input$spqc_group_label)[1])
      if (is.na(spqc_label) || !nzchar(spqc_label)) spqc_label <- "SPQC"
      spqc_prefix <- trimws(as.character(input$spqc_group_prefix)[1])
      if (is.na(spqc_prefix) || !nzchar(spqc_prefix)) spqc_prefix <- "SPQC"

      spqc_values <- rep(spqc_label, nrow(built))
      if (identical(spqc_mode, "batch")) {
        batch_values <- trimws(as.character(built$Batch))
        batch_values[is.na(batch_values) | !nzchar(batch_values)] <- run_label_date_token(built$Sample)[is.na(batch_values) | !nzchar(batch_values)]
        batch_values[is.na(batch_values) | !nzchar(batch_values)] <- "UnknownBatch"
        spqc_values <- paste(spqc_prefix, batch_values, sep = "_")
      } else if (identical(spqc_mode, "date")) {
        date_values <- run_label_date_token(built$Sample)
        date_values[is.na(date_values) | !nzchar(date_values)] <- "UnknownDate"
        spqc_values <- paste(spqc_prefix, date_values, sep = "_")
      }

      if (!identical(spqc_mode, "keep")) {
        built$Group[spqc_rows] <- spqc_values[spqc_rows]
        built$Condition[spqc_rows] <- spqc_values[spqc_rows]
        if ("Label" %in% colnames(built)) built$Label[spqc_rows] <- spqc_values[spqc_rows]
      }
    }

    committed_spqc_edits <- if (isTRUE(spqc_clear_pending())) empty_spqc_metadata_edits() else spqc_metadata_edits()
    built <- apply_proteomics_metadata_cell_edits(built, committed_spqc_edits)

    condition <- if ("Condition" %in% colnames(built)) as.character(built$Condition) else rep("", nrow(built))
    replicate <- if ("Replicate" %in% colnames(built)) as.character(built$Replicate) else rep("", nrow(built))
    label_value <- ifelse(!is.na(built$SampleName) & built$SampleName != "", built$SampleName, replicate)
    label_value[grepl("SPQC", built$Sample, ignore.case = TRUE) & (is.na(label_value) | label_value == "")] <- "QC"
    built$AnalysisLabel <- ifelse(condition != "" & !is.na(condition), paste(condition, label_value, sep = "_"), label_value)
    built$AnalysisLabel[is.na(built$AnalysisLabel) | built$AnalysisLabel == ""] <- built$Sample[is.na(built$AnalysisLabel) | built$AnalysisLabel == ""]

    excluded <- sample_exclusion_flags(built, input$sample_exclusions_text)
    built$Excluded <- excluded
    reason <- input$sample_exclusion_reason
    if (is.null(reason) || !nzchar(trimws(reason))) reason <- "Excluded sample"
    built$ExclusionReason <- ifelse(excluded, reason, "")

    built
  })

  built_metadata <- reactive({
    metadata_apply_revision()
    state <- applied_metadata_state()
    validate(need(!is.null(state$metadata), "Apply metadata changes on the Make metadata tab first."))
    state$metadata
  })

  metadata_default_columns <- function(columns) {
    preferred <- c(
      "RunOrder", "SampleName", "Condition", "Group", "Batch", "AnalysisLabel", "Excluded", "ExclusionReason", "Run Label",
      "File Name", "Fraction", "Quantity Correction Factor", "Reference", "Color"
    )
    preferred[preferred %in% columns]
  }

  update_protein_header_label_columns <- function(columns, saved_selected = NULL) {
    columns <- table_s1_metadata_choices(columns, downstream_metadata_columns())
    preferred <- c("Condition", "SampleName", "Replicate", "AnalysisLabel", "Run Label", "Sample")
    choices <- unique(c(preferred[preferred %in% columns], columns))
    current <- isolate(input$protein_header_label_columns)
    selected <- current[current %in% choices]
    saved_selected <- unlist(saved_selected, use.names = FALSE)
    saved_selected <- saved_selected[saved_selected %in% choices]
    if (length(selected) == 0 && length(saved_selected) > 0) selected <- saved_selected
    if (length(selected) == 0) {
      selected <- if ("SampleName" %in% choices) {
        "SampleName"
      } else if ("AnalysisLabel" %in% choices) {
        "AnalysisLabel"
      } else {
        character(0)
      }
    }
    updateSelectizeInput(
      session,
      "protein_header_label_columns",
      choices = choices,
      selected = selected,
      server = FALSE
    )
  }

  update_protein_quantity_order_columns <- function(columns, saved_selected = NULL, prefer_saved = FALSE) {
    columns <- table_s1_metadata_choices(columns, downstream_metadata_columns())
    choices <- setdiff(columns, c("SampleDetailsID", "Excluded", "ExclusionReason"))
    current <- isolate(input$protein_quantity_order_columns)
    selected <- current[current %in% choices]
    saved_selected <- unlist(saved_selected, use.names = FALSE)
    saved_selected <- saved_selected[saved_selected %in% choices]
    if (isTRUE(prefer_saved) && length(saved_selected) > 0) selected <- saved_selected
    if (length(selected) == 0 && length(saved_selected) > 0) selected <- saved_selected
    if (length(selected) == 0) {
      defaults <- c("Condition", "Batch", "Replicate", "RunOrder")
      selected <- defaults[defaults %in% choices]
    }
    updateSelectizeInput(
      session,
      "protein_quantity_order_columns",
      choices = choices,
      selected = selected,
      server = TRUE
    )
  }

observeEvent(draft_metadata(), {
  columns <- colnames(draft_metadata())
    current <- isolate(input$metadata_export_columns)
    selected <- current[current %in% columns]
    saved <- unlist(restored_project_settings()$metadata_export_columns, use.names = FALSE)
    saved <- saved[saved %in% columns]
    if (is.null(current) && length(saved) > 0) selected <- saved
    if (is.null(current) && length(selected) == 0) selected <- metadata_default_columns(columns)
    updateSelectizeInput(
      session,
      "metadata_export_columns",
      choices = columns,
      selected = selected,
      server = FALSE
    )
}, ignoreInit = FALSE)

  downstream_metadata_columns <- reactive({
    state <- applied_metadata_state()
    validate(need(!is.null(state$metadata), "Apply metadata changes on the Make metadata tab first."))
    table_s1_metadata_choices(colnames(state$metadata), state$columns)
  })

  scoped_metadata_columns <- function(md) {
    table_s1_metadata_choices(colnames(md), downstream_metadata_columns())
  }

  observeEvent(list(built_metadata(), downstream_metadata_columns()), {
    columns <- downstream_metadata_columns()
    update_protein_header_label_columns(columns, restored_project_settings()$protein_header_label_columns)
    update_protein_quantity_order_columns(columns, restored_project_settings()$protein_quantity_order_columns)
  }, ignoreInit = FALSE)

  observeEvent(input$reset_metadata_columns, {
    columns <- colnames(draft_metadata())
    updateSelectizeInput(
      session,
      "metadata_export_columns",
      choices = columns,
      selected = metadata_default_columns(columns),
      server = FALSE
    )
  })

  exported_metadata <- reactive({
    built <- built_metadata()
    selected <- applied_metadata_state()$columns
    validate(need(length(selected) > 0, "Select at least one column for Table S1. Metadata."))
    selected <- selected[selected %in% colnames(built)]
    validate(need(length(selected) > 0, "No selected columns are available in the built metadata."))

    if ("RunOrder" %in% colnames(built)) {
      built <- built %>% arrange(is.na(RunOrder), RunOrder)
    }
    built[, selected, drop = FALSE]
  })

  draft_exported_metadata <- reactive({
    built <- draft_metadata()
    built <- apply_proteomics_metadata_cell_edits(built, spqc_metadata_draft_edits())
    selected <- input$metadata_export_columns
    validate(need(length(selected) > 0, "Select at least one column for Table S1. Metadata."))
    selected <- selected[selected %in% colnames(built)]
    validate(need(length(selected) > 0, "No selected columns are available in the draft metadata."))
    if ("RunOrder" %in% colnames(built)) built <- built %>% arrange(is.na(RunOrder), RunOrder)
    built[, selected, drop = FALSE]
  })

  update_protein_condition_controls <- function(saved_settings = list(), update_group_selector = TRUE, prefer_saved_pair = FALSE) {
    md <- active_metadata()
    pair_choices <- scoped_metadata_columns(md)
    current_pair <- isolate(input$stats_pair_col)
    saved_pair <- unlist(saved_settings$stats_pair_col, use.names = FALSE)
    selected_pair <- retain_metadata_choice(if (isTRUE(prefer_saved_pair) && length(saved_pair)) saved_pair else current_pair, pair_choices, c(saved_pair, "Replicate"))
    updateSelectInput(session, "stats_pair_col", choices = c("Select pairing field..." = "", stats::setNames(pair_choices, pair_choices)), selected = selected_pair)
    cv_group_columns <- table_s1_metadata_choices(cv_metadata_group_columns(md), downstream_metadata_columns())
    current_cv_group_col <- isolate(input$protein_cv_group_col)
    saved_cv_group_col <- as.character(saved_settings$protein_cv_group_col)[1L]
    selected_cv_group_col <- retain_metadata_choice(current_cv_group_col, cv_group_columns, c(saved_cv_group_col, "Condition"))
    updateSelectInput(session, "protein_cv_group_col", choices = stats::setNames(cv_group_columns, cv_group_columns), selected = selected_cv_group_col)
    conditions <- if (length(selected_cv_group_col) && nzchar(selected_cv_group_col) && selected_cv_group_col %in% colnames(md)) {
      values <- trimws(as.character(md[[selected_cv_group_col]]))
      sort(unique(values[!is.na(values) & nzchar(values)]))
    } else character(0)
    current <- isolate(input$cv_conditions)
    selected <- current[current %in% conditions]
    saved_cv <- unlist(saved_settings$cv_conditions, use.names = FALSE)
    saved_cv <- saved_cv[saved_cv %in% conditions]
    if (length(selected) == 0 && length(saved_cv) > 0) selected <- saved_cv
    if (length(selected) == 0) {
      spqc_conditions <- conditions[grepl("SPQC", conditions, ignore.case = TRUE)]
      if (length(spqc_conditions) > 0) selected <- spqc_conditions
    }
    updateSelectizeInput(session, "cv_conditions", choices = conditions, selected = selected, server = FALSE)
    candidate_columns <- table_s1_metadata_choices(stats_metadata_group_columns(md), downstream_metadata_columns())
    current_group_columns <- isolate(input$stats_group_columns)
    selected_group_columns <- current_group_columns[current_group_columns %in% candidate_columns]
    saved_group_columns <- unlist(saved_settings$stats_group_columns, use.names = FALSE)
    saved_group_columns <- saved_group_columns[saved_group_columns %in% candidate_columns]
    if (is.null(current_group_columns) && !length(selected_group_columns)) selected_group_columns <- saved_group_columns
    if (isTRUE(update_group_selector)) {
      updateSelectizeInput(session, "stats_group_columns", choices = candidate_columns, selected = selected_group_columns, server = FALSE)
    }

    comparison_ids <- unlist(lapply(selected_group_columns, function(group_col) {
      group_values <- unique(trimws(as.character(md[[group_col]])))
      group_values <- group_values[!is.na(group_values) & nzchar(group_values) & !grepl("SPQC", group_values, ignore.case = TRUE)]
      if (length(group_values) < 2L) return(character(0))
      unlist(lapply(group_values, function(numerator) {
        vapply(group_values[group_values != numerator], function(denominator) make_stats_comparison_id(group_col, numerator, denominator), character(1))
      }))
    }))
    current_comparisons <- as.character(isolate(input$stats_comparisons))
    saved_comparisons <- unlist(saved_settings$stats_comparisons, use.names = FALSE)
    legacy_comparisons <- unique(c(current_comparisons, saved_comparisons))
    legacy_comparisons <- legacy_comparisons[vapply(legacy_comparisons, function(value) {
      parsed <- tryCatch(parse_stats_comparison_id(value), error = function(e) NULL)
      if (is.null(parsed) || !isTRUE(parsed$legacy) || !"Condition" %in% selected_group_columns) return(FALSE)
      parsed$numerator %in% as.character(md$Condition) && parsed$denominator %in% as.character(md$Condition)
    }, logical(1))]
    if (length(legacy_comparisons)) {
      legacy_new_ids <- vapply(legacy_comparisons, function(value) {
        parsed <- parse_stats_comparison_id(value)
        make_stats_comparison_id("Condition", parsed$numerator, parsed$denominator)
      }, character(1))
      for (i in seq_along(legacy_comparisons)) comparison_ids[comparison_ids == legacy_new_ids[[i]]] <- legacy_comparisons[[i]]
    }
    comparison_ids <- unique(as.character(comparison_ids))
    comparison_labels <- vapply(comparison_ids, stats_comparison_label, character(1))
    comparison_choices <- stats::setNames(comparison_ids, comparison_labels)
    chosen <- current_comparisons[current_comparisons %in% comparison_ids]
    saved_comparisons <- saved_comparisons[saved_comparisons %in% comparison_ids]
    if (length(chosen) == 0 && length(saved_comparisons) > 0) chosen <- saved_comparisons
    if (length(chosen) == 0) {
      defaults <- c("Condition|||CKDu|||CKD", "Condition|||CKD|||Control", "Condition|||CKDu|||Control")
      chosen <- defaults[defaults %in% comparison_ids]
      if (length(chosen) == 0 && length(comparison_ids) > 0) chosen <- comparison_ids[1]
    }
    updateSelectizeInput(session, "stats_comparisons", choices = comparison_choices, selected = chosen, server = TRUE)
    paired_current <- isolate(input$stats_paired_comparisons)
    paired_selected <- paired_current[paired_current %in% chosen]
    saved_paired <- unlist(saved_settings$stats_paired_comparisons, use.names = FALSE)
    saved_paired <- saved_paired[saved_paired %in% chosen]
    if (length(paired_selected) == 0 && length(saved_paired) > 0) paired_selected <- saved_paired
    updateSelectizeInput(
      session,
      "stats_paired_comparisons",
      choices = stats::setNames(chosen, vapply(chosen, stats_comparison_label, character(1))),
      selected = paired_selected,
      server = FALSE
    )
  }

  observeEvent(list(built_metadata(), downstream_metadata_columns()), {
    update_protein_condition_controls(restored_project_settings())
  }, ignoreInit = FALSE)

  observeEvent(list(built_metadata(), downstream_metadata_columns()), {
    choices <- table_s1_metadata_choices(cv_metadata_group_columns(active_metadata()), downstream_metadata_columns())
    current <- isolate(input$cv_plot_group_col)
    saved <- as.character(restored_project_settings()$cv_plot_group_col)[1L]
    selected <- retain_metadata_choice(current, choices, c(saved, "Condition"))
    updateSelectInput(session, "cv_plot_group_col", choices = c("Select grouping field..." = "", stats::setNames(choices, choices)), selected = selected)
  }, ignoreInit = FALSE)

  observeEvent(input$stats_group_columns, {
    update_protein_condition_controls(restored_project_settings(), update_group_selector = FALSE)
  }, ignoreInit = TRUE)

  observeEvent(input$protein_cv_group_col, {
    update_protein_condition_controls(restored_project_settings())
  }, ignoreInit = TRUE)

  observeEvent(input$stats_comparisons, {
    comparisons <- input$stats_comparisons
    if (is.null(comparisons)) comparisons <- character(0)
    labels <- vapply(comparisons, stats_comparison_label, character(1))
    current <- isolate(input$stats_paired_comparisons)
    selected <- current[current %in% comparisons]
    updateSelectizeInput(
      session,
      "stats_paired_comparisons",
      choices = stats::setNames(comparisons, labels),
      selected = selected,
      server = FALSE
    )
  }, ignoreInit = FALSE)

  protein_non_data_columns_from_table <- function(report) {
    proteomics_annotation_columns(colnames(report))
  }

  protein_non_data_columns <- function(file_info) {
    protein_non_data_columns_from_table(read_uploaded_table(file_info))
  }

  observeEvent(list(imported_project_files(), input$protein_no_impute_file), {
    if (is.null(project_file("protein_no_impute_file"))) return()
    columns <- protein_non_data_columns(project_file("protein_no_impute_file"))
    current <- isolate(input$s2_non_data_columns)
    selected <- current[current %in% columns]
    if (length(selected) == 0) selected <- columns
    updateSelectizeInput(session, "s2_non_data_columns", choices = columns, selected = selected, server = TRUE)
  }, ignoreInit = TRUE)

  observeEvent(list(imported_project_files(), input$protein_imputed_file), {
    if (is.null(project_file("protein_imputed_file"))) return()
    columns <- protein_non_data_columns(project_file("protein_imputed_file"))
    current <- isolate(input$s3_non_data_columns)
    selected <- current[current %in% columns]
    if (length(selected) == 0) selected <- columns
    updateSelectizeInput(session, "s3_non_data_columns", choices = columns, selected = selected, server = TRUE)
  }, ignoreInit = TRUE)

  sync_protein_non_data_inputs <- function() {
    if (isolate(protein_source_available("no_impute"))) {
      columns <- protein_non_data_columns_from_table(isolate(protein_source_table("no_impute")))
      current <- isolate(input$s2_non_data_columns)
      selected <- current[current %in% columns]
      if (length(selected) == 0) selected <- columns
      updateSelectizeInput(session, "s2_non_data_columns", choices = columns, selected = selected, server = TRUE)
    }
    if (isolate(protein_source_available("imputed"))) {
      columns <- protein_non_data_columns_from_table(isolate(protein_source_table("imputed")))
      current <- isolate(input$s3_non_data_columns)
      selected <- current[current %in% columns]
      if (length(selected) == 0) selected <- columns
      updateSelectizeInput(session, "s3_non_data_columns", choices = columns, selected = selected, server = TRUE)
    }
  }

  sync_metadata_column_inputs <- function(settings = list(), prefer_saved = FALSE) {
    built <- built_metadata()
    columns <- colnames(built)
    current_metadata <- isolate(input$metadata_export_columns)
    selected_metadata <- current_metadata[current_metadata %in% columns]
    saved_metadata <- unlist(settings$metadata_export_columns, use.names = FALSE)
    saved_metadata <- saved_metadata[saved_metadata %in% columns]
    if (isTRUE(prefer_saved) && length(saved_metadata) > 0) selected_metadata <- saved_metadata
    if (is.null(current_metadata) && length(selected_metadata) == 0) selected_metadata <- metadata_default_columns(columns)
    updateSelectizeInput(session, "metadata_export_columns", choices = columns, selected = selected_metadata, server = FALSE)

    update_protein_header_label_columns(columns, settings$protein_header_label_columns)
    update_protein_quantity_order_columns(columns, settings$protein_quantity_order_columns, prefer_saved = prefer_saved)
  }

  sync_condition_dependent_inputs <- function(settings = list(), prefer_saved_pair = FALSE) {
    update_protein_condition_controls(settings, prefer_saved_pair = prefer_saved_pair)
  }

  sync_restored_project_inputs <- function(settings = list()) {
    try(restore_project_settings(settings), silent = TRUE)
    try(sync_metadata_column_inputs(settings, prefer_saved = TRUE), silent = TRUE)
    try(sync_condition_dependent_inputs(settings, prefer_saved_pair = TRUE), silent = TRUE)
    try(sync_protein_non_data_inputs(), silent = TRUE)
  }

  refresh_restored_project_controls <- function(settings = restored_project_settings()) {
    sync_restored_project_inputs(settings)
    session$onFlushed(function() {
      isolate(sync_restored_project_inputs(settings))
      session$onFlushed(function() {
        isolate(sync_restored_project_inputs(settings))
      }, once = TRUE)
    }, once = TRUE)
  }

  observeEvent(project_restore_token(), {
    refresh_restored_project_controls(restored_project_settings())
  }, ignoreInit = TRUE)

  observeEvent(input$refresh_project_selectors, {
    refresh_restored_project_controls(restored_project_settings())
    project_bundle_message(paste0(project_bundle_message(), "\nRefreshed restored project controls."))
  })

  observeEvent(input$workflow_tabs, {
    if (!identical(input$workflow_tabs, "Protein tables")) return()
    settings <- restored_project_settings()
    update_protein_header_label_columns(downstream_metadata_columns(), settings$protein_header_label_columns)
    update_protein_quantity_order_columns(downstream_metadata_columns(), settings$protein_quantity_order_columns)
    sync_condition_dependent_inputs(settings)
    session$onFlushed(function() {
      isolate(update_protein_header_label_columns(downstream_metadata_columns(), settings$protein_header_label_columns))
      isolate(update_protein_quantity_order_columns(downstream_metadata_columns(), settings$protein_quantity_order_columns))
      isolate(sync_condition_dependent_inputs(settings))
    }, once = TRUE)
  }, ignoreInit = TRUE)

  metadata_sorted_samples <- function(md, order_columns) {
    order_columns <- order_columns[order_columns %in% colnames(md)]
    if (length(order_columns) == 0) return(character(0))
    sortable <- md
    sortable$.input_order <- seq_len(nrow(sortable))
    order_args <- lapply(order_columns, function(column_name) {
      values <- sortable[[column_name]]
      numeric_values <- suppressWarnings(as.numeric(as.character(values)))
      finite_count <- sum(is.finite(numeric_values))
      if (finite_count >= max(2, floor(0.8 * length(values)))) {
        numeric_values[!is.finite(numeric_values)] <- Inf
        numeric_values
      } else {
        text_values <- as.character(values)
        text_values[is.na(text_values) | !nzchar(text_values)] <- "zzzz_missing"
        text_values
      }
    })
    order_index <- do.call(order, c(order_args, list(sortable$.input_order, na.last = TRUE)))
    as.character(sortable$Sample[order_index])
  }

  reorder_batch_corrected_protein_table <- function(table, sample_map, suffix = "_batch_corrected_Protein_group_abundance") {
    sample_order <- metadata_sorted_samples(active_metadata(), input$protein_quantity_order_columns)
    if (length(sample_order) == 0 || is.null(sample_map) || nrow(sample_map) == 0) return(table)
    sample_map <- sample_map[match(sample_order, sample_map$Sample), , drop = FALSE]
    sample_map <- sample_map[!is.na(sample_map$Sample), , drop = FALSE]
    ordered_abundance <- paste0(sample_map$HeaderLabel, suffix)
    ordered_abundance <- ordered_abundance[ordered_abundance %in% colnames(table)]
    if (length(ordered_abundance) == 0) return(table)
    all_abundance <- colnames(table)[endsWith(colnames(table), suffix)]
    derived_cols <- grep("(_percent_CV|_log2_fold_change|_(paired|unpaired)_t_test_p_value|_BH_FDR)$", colnames(table), value = TRUE)
    info_cols <- setdiff(colnames(table), c(all_abundance, derived_cols))
    table[, c(info_cols, ordered_abundance, derived_cols), drop = FALSE]
  }

  renamed_protein_table <- function(file_info, add_stats = FALSE, selected_non_data = NULL) {
    if (isTRUE(add_stats)) protein_stats_refresh_revision()
    report <- if (is.data.frame(file_info)) {
      as.data.frame(file_info, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      read_uploaded_table(file_info)
    }
    md <- active_metadata()
    validate(need("Sample" %in% colnames(md), "Build metadata before adding protein tables."))
    validate(need(!anyDuplicated(md$Sample), "Metadata Run Label values must be unique before renaming protein columns."))

    metadata_text <- function(column, default = "") {
      if (column %in% colnames(md)) {
        values <- as.character(md[[column]])
        values[is.na(values)] <- default
        if (identical(column, "SampleName")) {
          missing_sample_name <- !nzchar(values)
          if (any(missing_sample_name) && "Replicate" %in% colnames(md)) {
            replicate_values <- as.character(md$Replicate)
            replicate_values[is.na(replicate_values)] <- ""
            values[missing_sample_name & nzchar(replicate_values)] <- replicate_values[missing_sample_name & nzchar(replicate_values)]
          }
        }
        values
      } else {
        rep(default, nrow(md))
      }
    }

    header_columns <- input$protein_header_label_columns
    header_columns <- header_columns[header_columns %in% colnames(md)]
    if (length(header_columns) == 0) {
      header_mode <- input$protein_header_label_mode
      header_columns <- if (identical(header_mode, "condition")) {
        c("Condition", "Replicate")
      } else if ("SampleName" %in% colnames(md)) {
        "SampleName"
      } else if ("AnalysisLabel" %in% colnames(md)) {
        "AnalysisLabel"
      } else {
        "Sample"
      }
    }

    header_parts <- lapply(header_columns, metadata_text)
    header_labels <- vapply(seq_len(nrow(md)), function(row_index) {
      parts <- vapply(header_parts, `[`, character(1), row_index)
      parts <- trimws(parts)
      parts <- parts[!is.na(parts) & nzchar(parts)]
      if (length(parts) == 0) return(as.character(md$Sample[row_index]))
      paste(parts, collapse = "_")
    }, character(1))
    header_labels <- trimws(header_labels)
    missing_header_labels <- is.na(header_labels) | !nzchar(header_labels)
    header_labels[missing_header_labels] <- as.character(md$Sample[missing_header_labels])
    validate(need(!anyDuplicated(header_labels), paste0(
      "Protein header labels are not unique when using ",
      paste(header_columns, collapse = "_"),
      ". Add another metadata variable, such as Replicate, or update metadata labels."
    )))

    label_by_run <- setNames(header_labels, as.character(md$Sample))
    original_names <- normalize_proteomics_text(colnames(report))
    validate(need(!anyDuplicated(original_names), "Protein table headers become duplicated after normalizing encoded ampersands."))
    colnames(report) <- original_names
    names_without_index <- sub("^\\[[0-9]+\\][[:space:]]*", "", original_names)
    measurement_kind <- proteomics_measurement_kind(names_without_index)
    precursor <- !is.na(measurement_kind) & measurement_kind == "precursor"
    abundance <- !is.na(measurement_kind) & measurement_kind == "abundance"
    processed_measurement <- grepl("_quantified_precursors$|_Protein_group_abundance$", names_without_index)
    derived <- grepl("(_percent_CV|_log2_fold_change|_(paired|unpaired)_t_test_p_value|_BH_FDR)$", names_without_index)
    run_labels <- names_without_index
    run_labels[precursor] <- sub("\\.PG\\.NrOfPrecursorsUsedForQuantification$", "", run_labels[precursor])
    run_labels[precursor] <- sub("_quantified_precursors$", "", run_labels[precursor])
    run_labels[abundance] <- sub("\\.PG\\.Quantity$", "", run_labels[abundance])
    run_labels[abundance] <- sub("_Protein_group_abundance$", "", run_labels[abundance])
    sample_ids <- run_labels
    sample_ids[processed_measurement] <- resolve_proteomics_processed_sample_ids(
      run_labels[processed_measurement],
      md,
      header_labels,
      project_db_cache()$processed_sample_map
    )
    matched <- (precursor | abundance) & !is.na(sample_ids) & sample_ids %in% names(label_by_run)

    replacement_names <- original_names
    replacement_names[matched & precursor] <- paste0(
      unname(label_by_run[sample_ids[matched & precursor]]),
      "_quantified_precursors"
    )
    replacement_names[matched & abundance] <- paste0(
      unname(label_by_run[sample_ids[matched & abundance]]),
      "_Protein_group_abundance"
    )
    validate(need(!anyDuplicated(replacement_names), "Renamed protein table headers are not unique. Review AnalysisLabel values."))
    colnames(report) <- replacement_names

    available_non_data <- original_names[!(precursor | abundance | derived)]
    selected_non_data <- selected_non_data[selected_non_data %in% available_non_data]
    if (length(selected_non_data) == 0) selected_non_data <- available_non_data
    measurement_index <- which(matched & (precursor | abundance))
    measurement_info <- data.frame(
      renamed = replacement_names[measurement_index],
      sample_id = sample_ids[measurement_index],
      measurement_type = ifelse(precursor[measurement_index], "precursor", "abundance"),
      source_order = seq_along(measurement_index),
      stringsAsFactors = FALSE
    )
    sample_order <- metadata_sorted_samples(md, input$protein_quantity_order_columns)
    if (length(sample_order) == 0) sample_order <- unique(measurement_info$sample_id)
    measurement_info$sample_rank <- match(measurement_info$sample_id, sample_order)
    ordered_measurements <- order_proteomics_measurement_columns(measurement_info$renamed, measurement_info$sample_rank)
    measurement_info <- measurement_info[match(ordered_measurements, measurement_info$renamed), , drop = FALSE]
    measurement_cols <- measurement_info$renamed
    report <- report[, c(selected_non_data, measurement_cols), drop = FALSE]

    selected_conditions <- input$cv_conditions
    if (is.null(selected_conditions)) selected_conditions <- character(0)
    cv_group_col <- input$protein_cv_group_col
    if (length(selected_conditions)) {
      validate(need(!is.null(cv_group_col) && length(cv_group_col) && nzchar(cv_group_col) && cv_group_col %in% colnames(md), "Select a metadata column for %CV calculation."))
    }
    protein_annotation_cols <- seq_along(selected_non_data)
    cv_df <- NULL
    if (length(selected_conditions) > 0) {
      cv_data <- lapply(selected_conditions, function(condition) {
        value_cols <- paste0(
          as.character(header_labels[as.character(md[[cv_group_col]]) == condition]),
          "_Protein_group_abundance"
        )
        value_cols <- value_cols[value_cols %in% colnames(report)]
        if (length(value_cols) < 2) return(rep(NaN, nrow(report)))

        abundances <- as.data.frame(report[, value_cols, drop = FALSE], stringsAsFactors = FALSE)
        abundances[] <- lapply(abundances, function(column) suppressWarnings(as.numeric(as.character(column))))
        apply(abundances, 1, function(values) {
          values <- values[is.finite(values)]
          mean_value <- mean(values)
          if (length(values) < 2 || !is.finite(mean_value) || mean_value <= 0) return(NaN)
          100 * stats::sd(values) / mean_value
        })
      })
      cv_df <- as.data.frame(cv_data, check.names = FALSE)
      colnames(cv_df) <- paste0(selected_conditions, "_percent_CV")
    }

    comparisons <- input$stats_comparisons
    if (is.null(comparisons)) comparisons <- character(0)
    stats_df <- NULL
    stats_comparison_labels <- character(0)
    stats_method_notes <- character(0)
    if (isTRUE(add_stats) && length(comparisons) > 0) {
      stats_blocks <- lapply(comparisons, function(comparison) {
        parsed_comparison <- parse_stats_comparison_id(comparison)
        group_col <- parsed_comparison$group_col
        numerator <- parsed_comparison$numerator
        denominator <- parsed_comparison$denominator
        validate(need(group_col %in% colnames(md), paste0("Statistics metadata column not found: ", group_col)))
        validate(need(numerator != denominator, "Each statistics comparison must have different numerator and denominator groups."))
        requested_paired <- comparison %in% input$stats_paired_comparisons
        pairing <- replicate_pair_plan(md, header_labels, numerator, denominator, group_col = group_col, pair_col = input$stats_pair_col)
        effective_paired <- requested_paired && pairing$balanced
        fallback_reason <- ""

        if (effective_paired) {
          numerator_cols <- paste0(pairing$numerator_labels, "_Protein_group_abundance")
          denominator_cols <- paste0(pairing$denominator_labels, "_Protein_group_abundance")
          if (!all(numerator_cols %in% colnames(report)) || !all(denominator_cols %in% colnames(report))) {
            effective_paired <- FALSE
            fallback_reason <- "matched abundance columns are missing"
          }
        }
        if (!effective_paired) {
          numerator_cols <- paste0(as.character(header_labels[as.character(md[[group_col]]) == numerator]), "_Protein_group_abundance")
          denominator_cols <- paste0(as.character(header_labels[as.character(md[[group_col]]) == denominator]), "_Protein_group_abundance")
          if (requested_paired && !pairing$balanced) fallback_reason <- pairing$reason
        }
        numerator_cols <- numerator_cols[numerator_cols %in% colnames(report)]
        denominator_cols <- denominator_cols[denominator_cols %in% colnames(report)]
        validate(need(length(numerator_cols) >= 2 && length(denominator_cols) >= 2, paste0("Comparison ", numerator, " vs ", denominator, " requires at least two matched abundance columns per group.")))

        result <- calculate_protein_comparison(report, numerator_cols, denominator_cols, paired = effective_paired)
        log2_fc <- result$data$log2_fc
        p_value <- result$data$p_value
        comparison_prefix <- stats_comparison_prefix(comparison)
        comparison_df <- data.frame(log2_fc, p_value, check.names = FALSE)
        colnames(comparison_df) <- c(
          paste0(comparison_prefix, "_log2_fold_change"),
          paste0(comparison_prefix, "_", result$method, "_t_test_p_value")
        )
        if (isTRUE(input$stats_bh_fdr)) {
          comparison_df[[paste0(comparison_prefix, "_BH_FDR")]] <- stats::p.adjust(p_value, method = "BH")
        }
        method_note <- paste0(stats_comparison_label(comparison), ": ", result$method)
        if (requested_paired && !effective_paired) {
          method_note <- paste0(method_note, " (automatic fallback: ", fallback_reason, ")")
        }
        attr(comparison_df, "method_note") <- method_note
        comparison_df
      })
      stats_method_notes <- vapply(stats_blocks, function(block) attr(block, "method_note"), character(1))
      stats_df <- do.call(cbind, stats_blocks)
      stats_comparison_labels <- vapply(comparisons, stats_comparison_label, character(1))
    }

    derived_order <- input$protein_derived_column_order
    if (is.null(derived_order) || !derived_order %in% c("cv_before_stats", "stats_before_cv")) {
      derived_order <- "cv_before_stats"
    }
    derived_blocks <- if (derived_order == "stats_before_cv") {
      list(stats_df, cv_df)
    } else {
      list(cv_df, stats_df)
    }
    derived_blocks <- derived_blocks[!vapply(derived_blocks, is.null, logical(1))]
    report_blocks <- c(
      list(report[, protein_annotation_cols, drop = FALSE]),
      derived_blocks,
      list(report[, setdiff(seq_len(ncol(report)), protein_annotation_cols), drop = FALSE])
    )
    report <- do.call(data.frame, c(report_blocks, list(check.names = FALSE)))

    attr(report, "matched_precursor") <- sum(matched & precursor)
    attr(report, "matched_abundance") <- sum(matched & abundance)
    attr(report, "sample_measurement_columns") <- sum(precursor | abundance)
    attr(report, "cv_conditions") <- selected_conditions
    attr(report, "stats_comparison") <- stats_comparison_labels
    attr(report, "stats_methods") <- stats_method_notes
    mapped_samples <- unique(as.character(measurement_info$sample_id))
    mapped_samples <- mapped_samples[!is.na(mapped_samples) & mapped_samples %in% as.character(md$Sample)]
    attr(report, "sample_map") <- data.frame(
      Sample = mapped_samples,
      HeaderLabel = unname(label_by_run[mapped_samples]),
      stringsAsFactors = FALSE
    )
    report
  }

  protein_no_impute_table <- reactive({
    cache <- project_db_cache()
    source <- protein_source_table("S2")
    revision <- protein_stats_refresh_revision()
    restored_cache <- is.null(project_file("protein_no_impute_file")) && !is.null(cache$processed_s2)
    if (restored_cache && use_cached_proteomics_table(source, revision)) return(source)
    renamed_protein_table(source, add_stats = "S2" %in% input$stats_tables, selected_non_data = input$s2_non_data_columns)
  })

  observe({
    md <- tryCatch(active_metadata(), error = function(e) NULL)
    if (is.null(md) || !nrow(md)) return()
    choices <- scoped_metadata_columns(md)
    current <- isolate(input$protein_knn_group_col)
    selected <- retain_metadata_choice(current, choices, c("Batch", "Condition"))
    updateSelectInput(
      session,
      "protein_knn_group_col",
      choices = c("Select metadata field..." = "", stats::setNames(choices, choices)),
      selected = selected
    )
  })

  observeEvent(list(input$protein_no_impute_file, imported_project_files()), {
    generated_s3_result(NULL)
  }, ignoreInit = TRUE)

  observeEvent(list(input$protein_knn_k, input$protein_knn_max_missing_percent, input$protein_knn_scope, input$protein_knn_group_col), {
    generated_s3_result(NULL)
  }, ignoreInit = TRUE)

  observeEvent(input$run_protein_knn, {
    tryCatch({
      source <- protein_source_table("S2")
      abundance_columns <- proteomics_abundance_columns(source)
      validate(need(length(abundance_columns) >= 2L, "Table S2 requires at least two protein abundance columns for kNN imputation."))
      column_groups <- NULL
      scope <- input$protein_knn_scope
      if (identical(scope, "metadata")) {
        md <- active_metadata()
        group_col <- input$protein_knn_group_col
        validate(need(!is.null(group_col) && nzchar(group_col) && group_col %in% colnames(md), "Select a metadata field for within-group kNN imputation."))
        sample_labels <- proteomics_abundance_sample_names(abundance_columns)
        processed <- all(endsWith(abundance_columns, "_Protein_group_abundance"))
        metadata_keys <- if (processed) protein_header_labels_from_metadata(md, input$protein_header_label_columns) else as.character(md$Sample)
        matched_rows <- match(sample_labels, metadata_keys)
        group_values <- trimws(as.character(md[[group_col]][matched_rows]))
        validate(need(!anyNA(matched_rows) && all(nzchar(group_values)), paste0("Every Table S2 abundance column must match metadata with a non-empty ", group_col, " value.")))
        column_groups <- split(abundance_columns, group_values)
      }
      result <- withProgress(message = "Generating Table S3 with kNN", value = 0.2, {
        result <- knn_impute_protein_report(
          source,
          abundance_columns = abundance_columns,
          k = input$protein_knn_k,
          max_missing_percent = input$protein_knn_max_missing_percent,
          column_groups = column_groups
        )
        incProgress(0.8, detail = "Table S3 ready")
        result
      })
      result$scope <- if (identical(scope, "metadata")) paste0("within ", input$protein_knn_group_col) else "all samples"
      generated_s3_result(result)
      columns <- protein_non_data_columns_from_table(result$data)
      current <- isolate(input$s3_non_data_columns)
      selected <- current[current %in% columns]
      if (!length(selected)) selected <- columns
      updateSelectizeInput(session, "s3_non_data_columns", choices = columns, selected = selected, server = TRUE)
      showNotification("kNN-imputed Table S3 generated.", type = "message", duration = 4)
    }, error = function(e) {
      generated_s3_result(NULL)
      showNotification(conditionMessage(e), type = "error", duration = NULL)
    })
  })

  output$protein_imputation_status <- renderText({
    method <- input$s3_imputation_source
    if (identical(method, "s2")) return("Table S3 source: Table S2 without additional imputation.")
    if (identical(method, "spectronaut")) {
      if (protein_source_available("S3")) return("Table S3 source: Spectronaut-imputed upload or saved project table.")
      return("Table S3 source: Spectronaut imputation. Select an imputed report above.")
    }
    result <- generated_s3_result()
    if (is.null(result)) {
      if (is.null(project_file("protein_no_impute_file")) && !is.null(project_db_cache()$processed_s3)) {
        return("Table S3 source: saved app-generated kNN table restored from the project.")
      }
      return("Table S3 source: app-generated kNN. Click Generate Table S3 with kNN.")
    }
    paste0(
      "Table S3 source: app-generated kNN (k = ", result$k, ", ", result$scope, "). ",
      nrow(result$data), " proteins retained; ", result$dropped_features, " excluded by the missingness limit; ",
      result$missing_before, " missing values imputed; ", result$missing_after, " remain."
    )
  })

  protein_imputed_table <- reactive({
    source <- protein_source_table("S3")
    revision <- protein_stats_refresh_revision()
    method <- input$s3_imputation_source
    restored_cache <- !is.null(project_db_cache()$processed_s3) &&
      ((identical(method, "spectronaut") && is.null(project_file("protein_imputed_file"))) ||
         (identical(method, "knn") && is.null(generated_s3_result()) && is.null(project_file("protein_no_impute_file"))))
    if (restored_cache && use_cached_proteomics_table(source, revision)) return(source)
    renamed_protein_table(source, add_stats = "S3" %in% input$stats_tables, selected_non_data = input$s3_non_data_columns)
  })

  output$protein_stats_status <- renderText({
    cache <- project_db_cache()
    has_cached_tables <- !is.null(cache$processed_s2) || !is.null(cache$processed_s3)
    if (has_cached_tables && protein_stats_refresh_revision() <= 0L) {
      "Showing saved protein tables and statistics from the opened project. Click Recalculate statistics after changing metadata or comparison settings."
    } else if (isTRUE(protein_stats_paused())) {
      "Statistics recalculation is stopped."
    } else {
      "Protein tables use the current statistics settings."
    }
  })

  observe({
    md <- tryCatch(built_metadata(), error = function(e) NULL)
    if (is.null(md)) return()
    choices <- setdiff(scoped_metadata_columns(md), c("SampleDetailsID"))
    current_batch <- isolate(input$batch_correction_batch_col)
    selected_batch <- if (!is.null(current_batch) && current_batch %in% choices) {
      current_batch
    } else {
      preferred_metadata_column(choices, c("Batch", "batch", "RunBatch", "run_batch", "RunOrder"))
    }
    current_group <- isolate(input$batch_correction_group_col)
    selected_group <- if (!is.null(current_group) && current_group %in% choices) {
      current_group
    } else {
      preferred_metadata_column(choices, c("Group", "group", "Condition", "condition"))
    }
    if (is.na(selected_batch) || !nzchar(selected_batch)) selected_batch <- choices[1]
    if (is.na(selected_group) || !nzchar(selected_group)) selected_group <- choices[1]
    selected_batch <- retain_metadata_choice(current_batch, choices, selected_batch)
    selected_group <- retain_metadata_choice(current_group, choices, selected_group)
    field_choices <- c("Select metadata field..." = "", stats::setNames(choices, choices))
    updateSelectInput(session, "batch_correction_batch_col", choices = field_choices, selected = selected_batch)
    updateSelectInput(session, "batch_correction_group_col", choices = field_choices, selected = selected_group)
  })

  batch_corrected_s3_computed <- eventReactive(input$run_batch_correction, {
    validate(
      need(requireNamespace("HarmonizR", quietly = TRUE) || requireNamespace("sva", quietly = TRUE), "Install HarmonizR or sva to run ComBat batch correction: BiocManager::install(c('HarmonizR', 'sva'))."),
      need(protein_source_available("S3"), "Upload Table S3 or open a DuckDB project containing it before running batch correction.")
    )
    source <- protein_source_table("S3")
    md <- active_metadata()
    prepared <- batch_correct_prepare_input(
      report = source,
      md = md,
      batch_col = input$batch_correction_batch_col,
      group_col = input$batch_correction_group_col,
      header_label_columns = input$protein_header_label_columns,
      feature_col_requested = input$report_feature_col,
      source_label = "S3",
      pseudocount = input$batch_correction_pseudocount,
      min_batches_for_feature = input$batch_correction_min_batches_for_feature
    )
    description <- prepared$description_as_input
    description$batch <- as.integer(factor(description$batch))
    combat_mode <- suppressWarnings(as.integer(input$batch_correction_combat_mode))
    if (!is.finite(combat_mode) || is.na(combat_mode)) combat_mode <- 1L
    par_prior <- !identical(combat_mode, 2L)
    correction_method <- if (par_prior) "sva ComBat parametric" else "sva ComBat non-parametric"
    harmonizr_error <- NULL
    corrected <- NULL
    sva_error <- NULL
    if (requireNamespace("sva", quietly = TRUE)) {
      design <- stats::model.matrix(~ group, data = description)
      corrected <- tryCatch(
        sva::ComBat(
          dat = as.matrix(prepared$data_as_input),
          batch = description$batch,
          mod = design,
          par.prior = par_prior,
          prior.plots = FALSE
        ),
        error = function(e) {
          sva_error <<- conditionMessage(e)
          NULL
        }
      )
    }
    if (is.null(corrected) && requireNamespace("HarmonizR", quietly = TRUE)) {
      corrected <- tryCatch(
        HarmonizR::harmonizR(
          data_as_input = prepared$data_as_input,
          description_as_input = description,
          algorithm = "ComBat",
          ComBat_mode = combat_mode,
          cores = 1
        ),
        error = function(e) {
          harmonizr_error <<- conditionMessage(e)
          NULL
        }
      )
      if (!is.null(corrected)) {
        correction_method <- "HarmonizR ComBat fallback"
      }
    }
    if (is.null(corrected)) {
      failure_details <- paste(
        c(
          if (!is.null(sva_error)) paste0("sva ComBat failed: ", sva_error) else if (!requireNamespace("sva", quietly = TRUE)) "sva is not installed",
          if (!is.null(harmonizr_error)) paste0("HarmonizR failed: ", harmonizr_error) else if (!requireNamespace("HarmonizR", quietly = TRUE)) "HarmonizR is not installed"
        ),
        collapse = "; "
      )
      validate(need(FALSE, paste0("Batch correction failed. ", failure_details)))
    }
    corrected_log2 <- as.data.frame(corrected, check.names = FALSE)
    corrected_log2[] <- lapply(corrected_log2, function(column) suppressWarnings(as.numeric(as.character(column))))
    corrected_log2 <- corrected_log2[, prepared$sample_map$Sample, drop = FALSE]
    corrected_abundance <- as.data.frame(lapply(corrected_log2, function(column) 2^column), check.names = FALSE)
    colnames(corrected_abundance) <- colnames(corrected_log2)
    expression <- data.frame(
      Feature = prepared$feature_info$Feature,
      corrected_log2,
      check.names = FALSE
    )
    table <- batch_correct_rebuild_table(source, corrected_abundance, prepared, suffix = "_batch_corrected_Protein_group_abundance")
    attr(table, "matched_precursor") <- 0
    attr(table, "matched_abundance") <- ncol(corrected_abundance)
    attr(table, "cv_conditions") <- character(0)
    result <- list(
      table = table,
      expression = expression,
      corrected_log2 = corrected_log2,
      corrected_abundance = corrected_abundance,
      prepared = prepared,
      combat_mode = combat_mode,
      correction_method = correction_method,
      cache_signature = batch_correction_cache_signature()
    )
    cached_batch_corrected_s3_result(result)
    batch_cache_message("Batch correction was run and is cached in memory. Save to the project database to reuse it after reopening the project.")
    result
  })

  batch_corrected_s3_result <- reactive({
    computed <- tryCatch(batch_corrected_s3_computed(), error = function(e) NULL)
    if (!is.null(computed)) return(computed)
    cached <- cached_batch_corrected_s3_result()
    validate(need(!is.null(cached), "Run batch correction first, or load a database project that contains a saved batch correction."))
    cached
  })

  batch_corrected_s3_table <- reactive({
    result <- batch_corrected_s3_result()
    table <- result$table
    if ("S3_batch_corrected" %in% input$stats_tables) {
      protein_stats_refresh_revision()
      log2_stats_base <- batch_correct_rebuild_table(
        protein_source_table("S3"),
        result$corrected_log2,
        result$prepared,
        suffix = "_batch_corrected_log2"
      )
      log2_stats_table <- append_log2_stats_to_protein_table(
        log2_stats_base,
        md = active_metadata(),
        sample_map = result$prepared$sample_map,
        comparisons = input$stats_comparisons,
        paired_comparisons = input$stats_paired_comparisons,
        pair_col = input$stats_pair_col,
        include_fdr = isTRUE(input$stats_bh_fdr),
        abundance_suffix = "_batch_corrected_log2"
      )
      stats_cols <- setdiff(colnames(log2_stats_table), colnames(log2_stats_base))
      if (length(stats_cols) > 0) {
        table <- data.frame(table, log2_stats_table[, stats_cols, drop = FALSE], check.names = FALSE)
      }
      attr(table, "stats_comparison") <- attr(log2_stats_table, "stats_comparison")
      attr(table, "stats_methods") <- attr(log2_stats_table, "stats_methods")
    }
    table <- reorder_batch_corrected_protein_table(table, result$prepared$sample_map)
    attr(table, "matched_precursor") <- 0
    attr(table, "matched_abundance") <- ncol(result$expression) - 1
    attr(table, "cv_conditions") <- character(0)
    table
  })

  batch_corrected_s3_expression_data <- reactive({
    batch_corrected_s3_result()$expression
  })

  observe({
    comparisons <- input$stats_comparisons
    if (length(comparisons) == 0) return()
    labels <- vapply(comparisons, stats_comparison_label, character(1))
    selected <- isolate(input$volcano_comparison)
    if (is.null(selected) || !selected %in% comparisons) selected <- comparisons[1]
    updateSelectInput(session, "volcano_comparison", choices = stats::setNames(comparisons, labels), selected = selected)

    gsea_selected <- isolate(input$gsea_comparison)
    if (is.null(gsea_selected) || !gsea_selected %in% comparisons) gsea_selected <- comparisons[1]
    updateSelectInput(session, "gsea_comparison", choices = stats::setNames(comparisons, labels), selected = gsea_selected)
  })

  volcano_source_file <- reactive({
    source <- input$volcano_source
    if (is.null(source) || !nzchar(source)) source <- "S3"
    validate(need(!identical(source, "S3_batch_corrected"), "Batch-corrected S3 is held in memory and does not have a source file."))
    if (identical(source, "S2")) {
      req(project_file("protein_no_impute_file"))
      project_file("protein_no_impute_file")
    } else {
      req(project_file("protein_imputed_file"))
      project_file("protein_imputed_file")
    }
  })

  volcano_source_table <- reactive({
    source <- input$volcano_source
    if (is.null(source) || !nzchar(source)) source <- "S3"
    if (identical(source, "S3_batch_corrected")) {
      batch_corrected_s3_table()
    } else {
      protein_source_table(source)
    }
  })

  observeEvent(list(input$volcano_source, project_restore_token(), imported_project_files(), input$protein_no_impute_file, input$protein_imputed_file, input$run_batch_correction), {
    label_options <- c(
      "Gene name" = "PG.Genes",
      "Protein name" = "PG.ProteinNames",
      "Protein accession / group" = "PG.ProteinGroups",
      "Protein accession" = "PG.ProteinAccessions",
      "Protein description" = "PG.ProteinDescriptions"
    )
    if (identical(input$volcano_source, "S3_batch_corrected")) {
      corrected <- tryCatch(batch_corrected_s3_table(), error = function(e) NULL)
      if (is.null(corrected)) {
        columns <- unname(label_options)
      } else {
        columns <- intersect(unname(label_options), colnames(corrected))
        if (length(columns) == 0) columns <- colnames(corrected)[1]
      }
    } else {
      columns <- intersect(unname(label_options), colnames(protein_source_table(input$volcano_source)))
    }
    selected <- isolate(input$volcano_label_col)
    if (is.null(selected) || !selected %in% columns) {
      preferred <- c("PG.Genes", "PG.ProteinGroups", "PG.ProteinNames")
      selected <- preferred[preferred %in% columns][1]
      if (is.na(selected)) selected <- columns[1]
    }
    choices <- label_options[label_options %in% columns]
    if (!length(choices)) choices <- stats::setNames(columns, columns)
    updateSelectInput(session, "volcano_label_col", choices = choices, selected = selected)
  }, ignoreInit = FALSE)

  volcano_stats_table <- reactive({
    if (input$volcano_source == "S3_batch_corrected") {
      batch_corrected_s3_table()
    } else if (input$volcano_source == "S2") {
      protein_no_impute_table()
    } else {
      protein_imputed_table()
    }
  })

  gsea_source_file <- reactive({
    source <- input$gsea_source
    if (is.null(source) || !nzchar(source)) source <- "S3"
    if (identical(source, "S2")) {
      req(project_file("protein_no_impute_file"))
      project_file("protein_no_impute_file")
    } else {
      req(project_file("protein_imputed_file"))
      project_file("protein_imputed_file")
    }
  })

  observeEvent(list(input$gsea_source, imported_project_files(), input$protein_no_impute_file, input$protein_imputed_file), {
    columns <- protein_non_data_columns_from_table(protein_source_table(input$gsea_source))
    selected <- isolate(input$gsea_gene_col)
    if (is.null(selected) || !selected %in% columns) {
      preferred <- c("PG.Genes", "Genes", "Gene", "PG.ProteinNames", "PG.ProteinGroups")
      selected <- preferred[preferred %in% columns][1]
      if (is.na(selected)) selected <- columns[1]
    }
    updateSelectInput(session, "gsea_gene_col", choices = columns, selected = selected)
  }, ignoreInit = FALSE)

  gsea_stats_table <- reactive({
    if (input$gsea_source == "S2") {
      protein_no_impute_table()
    } else {
      protein_imputed_table()
    }
  })

  msigdb_gene_sets <- function(species, collection_id) {
    validate(need(requireNamespace("msigdbr", quietly = TRUE), "Install the R package 'msigdbr' to load human or mouse gene sets."))
    parts <- strsplit(collection_id, ":", fixed = TRUE)[[1]]
    collection <- parts[1]
    subcollection <- if (length(parts) > 1) paste(parts[-1], collapse = ":") else NULL
    msig <- tryCatch(
      {
        if (is.null(subcollection)) {
          msigdbr::msigdbr(species = species, collection = collection)
        } else {
          msigdbr::msigdbr(species = species, collection = collection, subcollection = subcollection)
        }
      },
      error = function(e) NULL
    )
    if (is.null(msig)) {
      msig <- tryCatch(
        {
          if (is.null(subcollection)) {
            msigdbr::msigdbr(species = species, category = collection)
          } else {
            msigdbr::msigdbr(species = species, category = collection, subcategory = subcollection)
          }
        },
        error = function(e) NULL
      )
    }
    validate(need(!is.null(msig) && all(c("gs_name", "gene_symbol") %in% colnames(msig)), "Could not load the selected MSigDB gene set collection. Try a different collection or update 'msigdbr'."))
    lapply(split(msig$gene_symbol, msig$gs_name), unique)
  }

  gene_set_enrichment_data <- eventReactive(input$run_gsea, {
    validate(
      need(requireNamespace("fgsea", quietly = TRUE), "Install the R package 'fgsea' to run preranked gene set enrichment."),
      need(!is.null(input$gsea_comparison) && nzchar(input$gsea_comparison), "Select a protein-group comparison before running enrichment.")
    )
    report <- gsea_stats_table()
    source <- protein_source_table(input$gsea_source)
    gene_col <- input$gsea_gene_col
    if (is.null(gene_col) || !gene_col %in% colnames(source)) gene_col <- "PG.Genes"
    validate(need(gene_col %in% colnames(source), "Select a gene symbol column for enrichment."))

    prefix <- stats_comparison_prefix(input$gsea_comparison)
    fold_change_col <- paste0(prefix, "_log2_fold_change")
    p_value_col <- comparison_p_value_column(report, prefix)
    fdr_col <- paste0(prefix, "_BH_FDR")
    validate(need(!is.na(p_value_col) && fold_change_col %in% colnames(report), "Statistics are unavailable for the selected enrichment comparison."))

    p_value <- suppressWarnings(as.numeric(report[[p_value_col]]))
    fdr <- if (fdr_col %in% colnames(report)) suppressWarnings(as.numeric(report[[fdr_col]])) else stats::p.adjust(p_value, method = "BH")
    proteins <- if ("PG.ProteinGroups" %in% colnames(source)) as.character(source[["PG.ProteinGroups"]]) else as.character(seq_len(nrow(source)))
    ranks <- ranked_gene_statistics(
      proteins = proteins,
      genes = source[[gene_col]],
      log2_fc = suppressWarnings(as.numeric(report[[fold_change_col]])),
      p_value = p_value,
      fdr = fdr,
      metric = input$gsea_rank_metric
    )
    validate(need(length(ranks) >= 10, "Fewer than 10 ranked genes are available. Check the gene-symbol column and statistics table."))
    ranks <- sort(ranks, decreasing = TRUE)
    pathways <- msigdb_gene_sets(input$gsea_species, input$gsea_collection)
    result <- suppressWarnings(fgsea::fgsea(
      pathways = pathways,
      stats = ranks,
      minSize = max(1, as.integer(input$gsea_min_size)),
      maxSize = max(1, as.integer(input$gsea_max_size))
    ))
    result <- as.data.frame(result, stringsAsFactors = FALSE)
    if ("leadingEdge" %in% colnames(result)) {
      result$leadingEdge <- vapply(result$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
    }
    if ("pathway" %in% colnames(result)) names(result)[names(result) == "pathway"] <- "GeneSet"
    if ("padj" %in% colnames(result)) names(result)[names(result) == "padj"] <- "BH_FDR"
    result$Direction <- ifelse(result$NES >= 0, "Enriched in numerator", "Enriched in denominator")
    result <- result[order(result$BH_FDR, -abs(result$NES), na.last = TRUE), , drop = FALSE]
    rownames(result) <- NULL
    attr(result, "ranked_genes") <- length(ranks)
    attr(result, "gene_sets") <- length(pathways)
    result
  })

  gene_set_enrichment_plot_obj <- reactive({
    data <- gene_set_enrichment_data()
    validate(need(nrow(data) > 0, "No enriched gene sets were returned. Try a different collection or relax gene-set size filters."))
    top_n <- max(1, as.integer(input$gsea_top_n))
    plot_data <- data %>%
      dplyr::filter(is.finite(NES), is.finite(BH_FDR)) %>%
      dplyr::arrange(BH_FDR, dplyr::desc(abs(NES))) %>%
      dplyr::slice_head(n = top_n)
    validate(need(nrow(plot_data) > 0, "No finite enrichment results are available to plot."))
    plot_data$GeneSet <- factor(plot_data$GeneSet, levels = rev(plot_data$GeneSet))
    ggplot(plot_data, aes(x = GeneSet, y = NES, fill = Direction)) +
      geom_col(width = 0.75) +
      coord_flip() +
      scale_fill_manual(values = c("Enriched in numerator" = "#B2182B", "Enriched in denominator" = "#2166AC")) +
      labs(
        title = "Gene Set Enrichment",
        subtitle = paste0(stats_comparison_label(input$gsea_comparison), "; ", input$gsea_species),
        x = NULL,
        y = "Normalized enrichment score",
        fill = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top",
        panel.grid.major.y = element_blank()
      )
  })

  volcano_plot_data <- reactive({
    req(input$volcano_comparison)
    report <- volcano_stats_table()
    prefix <- stats_comparison_prefix(input$volcano_comparison)
    fold_change_col <- paste0(prefix, "_log2_fold_change")
    p_value_col <- comparison_p_value_column(report, prefix)
    fdr_col <- paste0(prefix, "_BH_FDR")
    validate(need(!is.na(p_value_col) && fold_change_col %in% colnames(report), "Select a statistics comparison for the volcano plot."))

    p_value <- suppressWarnings(as.numeric(report[[p_value_col]]))
    fdr <- if (!isTRUE(input$stats_bh_fdr)) {
      rep(NA_real_, length(p_value))
    } else if (fdr_col %in% colnames(report)) {
      suppressWarnings(as.numeric(report[[fdr_col]]))
    } else {
      stats::p.adjust(p_value, method = "BH")
    }
    metric_name <- if (isTRUE(input$stats_bh_fdr)) input$volcano_significance_metric else "p_value"
    metric <- if (identical(metric_name, "BH_FDR")) fdr else p_value
    source <- volcano_source_table()
    label_col <- input$volcano_label_col
    if (is.null(label_col) || !label_col %in% colnames(source)) label_col <- colnames(source)[1]
    labels <- as.character(source[[label_col]])
    labels[is.na(labels) | labels == ""] <- as.character(source[[1]])[is.na(labels) | labels == ""]
    protein_names <- if ("PG.ProteinNames" %in% colnames(source)) as.character(source[["PG.ProteinNames"]]) else rep(NA_character_, nrow(source))
    protein_descriptions <- if ("PG.ProteinDescriptions" %in% colnames(source)) as.character(source[["PG.ProteinDescriptions"]]) else rep(NA_character_, nrow(source))

    out <- data.frame(
      Protein = labels,
      ProteinGroupID = if ("PG.ProteinGroups" %in% colnames(source)) as.character(source[["PG.ProteinGroups"]]) else as.character(seq_len(nrow(source))),
      ProteinName = protein_names,
      ProteinDescription = protein_descriptions,
      Log2FoldChange = suppressWarnings(as.numeric(report[[fold_change_col]])),
      PValue = p_value,
      BH_FDR = fdr,
      Significance = metric,
      stringsAsFactors = FALSE
    )
    out$MinusLog10Significance <- -log10(pmax(out$Significance, .Machine$double.xmin))
    out$Status <- "Not significant"
    significant <- is.finite(out$Significance) &
      out$Significance <= input$volcano_sig_cutoff &
      is.finite(out$Log2FoldChange) &
      abs(out$Log2FoldChange) >= input$volcano_fc_cutoff
    out$Status[significant & out$Log2FoldChange > 0] <- "Increased"
    out$Status[significant & out$Log2FoldChange < 0] <- "Decreased"
    out
  })

  volcano_hits_data <- reactive({
    volcano_plot_data() %>%
      dplyr::filter(Status != "Not significant") %>%
      dplyr::arrange(Significance, dplyr::desc(abs(Log2FoldChange)))
  })

  highlighted_volcano_hits <- reactive({
    rows <- input$volcano_hits_table_rows_selected
    if (is.null(rows) || !length(rows)) return(data.frame())
    hits <- volcano_hits_data()
    rows <- rows[rows >= 1L & rows <= nrow(hits)]
    hits[rows, , drop = FALSE]
  })

  highlighted_volcano_features <- reactive({
    hits <- highlighted_volcano_hits()
    if (!nrow(hits)) return(character(0))
    source <- volcano_source_to_feature_source(input$volcano_source)
    info <- protein_info_for_source(source, input$report_feature_col)
    map_volcano_hits_to_features(hits, info)
  })

  pending_feature_selection <- reactiveVal(NULL)
  pending_box_selection <- reactiveVal(NULL)

  observeEvent(input$volcano_hits_table_rows_selected, {
    features <- highlighted_volcano_features()
    if (!length(features)) return()
    source <- volcano_source_to_feature_source(input$volcano_source)
    pending_feature_selection(list(source = source, features = features))
    pending_box_selection(list(source = source, features = features))
    updateRadioButtons(session, "feature_data_source", selected = source)
    updateRadioButtons(session, "script_box_source", selected = source)
  }, ignoreInit = TRUE)

  volcano_counts_data <- reactive({
    comparisons <- input$stats_comparisons
    validate(need(length(comparisons) > 0, "Select at least one protein-group comparison to calculate significant-protein counts."))
    report <- volcano_stats_table()
    rows <- lapply(comparisons, function(comparison) {
      parsed_comparison <- parse_stats_comparison_id(comparison)
      groups <- c(parsed_comparison$numerator, parsed_comparison$denominator)
      prefix <- stats_comparison_prefix(comparison)
      fold_change_col <- paste0(prefix, "_log2_fold_change")
      p_value_col <- comparison_p_value_column(report, prefix)
      fdr_col <- paste0(prefix, "_BH_FDR")
      validate(need(!is.na(p_value_col) && fold_change_col %in% colnames(report), paste0("Statistics unavailable for comparison: ", paste(groups, collapse = " vs "))))
      fold_change <- suppressWarnings(as.numeric(report[[fold_change_col]]))
      p_value <- suppressWarnings(as.numeric(report[[p_value_col]]))
      fdr <- if (fdr_col %in% colnames(report)) suppressWarnings(as.numeric(report[[fdr_col]])) else stats::p.adjust(p_value, method = "BH")
      metric_values <- volcano_significance_metric_values(p_value, fdr, input$stats_bh_fdr)
      dplyr::bind_rows(lapply(names(metric_values), function(metric_name) {
        significance <- metric_values[[metric_name]]
        is_significant <- is.finite(significance) &
          significance <= input$volcano_sig_cutoff &
          is.finite(fold_change) &
          abs(fold_change) >= input$volcano_fc_cutoff
        data.frame(
          ComparisonID = comparison,
          MetadataColumn = parsed_comparison$group_col,
          Comparison = stats_comparison_label(comparison),
          Numerator = groups[1],
          Denominator = groups[2],
          SignificanceMetric = metric_name,
          FoldChangeCutoff = input$volcano_fc_cutoff,
          SignificanceCutoff = input$volcano_sig_cutoff,
          Increased = sum(is_significant & fold_change > 0, na.rm = TRUE),
          Decreased = sum(is_significant & fold_change < 0, na.rm = TRUE),
          TotalSignificant = sum(is_significant, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }))
    })
    dplyr::bind_rows(rows)
  })

  volcano_count_subtitle <- function(data) {
    paste0(
      "Significant proteins: ",
      sum(data$Status == "Increased", na.rm = TRUE),
      " increased; ",
      sum(data$Status == "Decreased", na.rm = TRUE),
      " decreased"
    )
  }

  observeEvent(input$volcano_counts_table_rows_selected, {
    selected_row <- input$volcano_counts_table_rows_selected[1]
    req(selected_row)
    counts <- volcano_counts_data()
    req(selected_row <= nrow(counts))
    selected <- counts[selected_row, , drop = FALSE]
    updateSelectInput(session, "volcano_comparison", selected = selected$ComparisonID)
    updateSelectInput(
      session,
      "volcano_significance_metric",
      selected = if (selected$SignificanceMetric == "BH FDR") "BH_FDR" else "p_value"
    )
    updateNumericInput(session, "volcano_fc_cutoff", value = selected$FoldChangeCutoff)
    updateNumericInput(session, "volcano_sig_cutoff", value = selected$SignificanceCutoff)
  })

  volcano_plot_obj <- reactive({
    data <- volcano_plot_data() %>%
      dplyr::filter(is.finite(Log2FoldChange), is.finite(MinusLog10Significance))
    validate(need(nrow(data) > 0, "No finite statistics are available for this volcano plot."))
    metric_label <- if (isTRUE(input$stats_bh_fdr) && input$volcano_significance_metric == "BH_FDR") "BH FDR" else "p-value"
    comparison_label <- stats_comparison_label(input$volcano_comparison)
    count_subtitle <- volcano_count_subtitle(data)
    cutoff_subtitle <- paste0(comparison_label, "; cutoffs: |log2 FC| >= ", input$volcano_fc_cutoff, " and ", metric_label, " <= ", input$volcano_sig_cutoff)

    plot <- ggplot(data, aes(x = Log2FoldChange, y = MinusLog10Significance, color = Status)) +
      geom_point(alpha = 0.75, size = input$volcano_point_size) +
      geom_vline(xintercept = c(-input$volcano_fc_cutoff, input$volcano_fc_cutoff), linetype = "dashed", color = "grey45") +
      geom_hline(yintercept = -log10(input$volcano_sig_cutoff), linetype = "dashed", color = "grey45") +
      scale_color_manual(
        values = c("Increased" = "#B2182B", "Decreased" = "#2166AC", "Not significant" = "#BDBDBD"),
        breaks = c("Increased", "Decreased", "Not significant")
      ) +
      labs(
        title = input$volcano_title,
        subtitle = paste(count_subtitle, cutoff_subtitle, sep = "\n"),
        x = paste0("log2 fold-change (", comparison_label, ")"),
        y = paste0("-log10(", metric_label, ")"),
        color = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top",
        panel.grid.minor = element_blank()
      )

    labels <- NULL
    if (identical(input$volcano_label_mode, "selected")) {
      selected_rows <- input$volcano_hits_table_rows_selected
      hits <- volcano_hits_data()
      if (length(selected_rows) > 0) {
        selected_rows <- selected_rows[selected_rows <= nrow(hits)]
        labels <- hits[selected_rows, , drop = FALSE]
      }
    } else if (identical(input$volcano_label_mode, "automatic") && input$volcano_max_labels > 0) {
      labels <- volcano_hits_data() %>%
        dplyr::filter(Status %in% c("Increased", "Decreased")) %>%
        dplyr::group_by(Status) %>%
        dplyr::slice_head(n = input$volcano_max_labels) %>%
        dplyr::ungroup()
    }
    if (!is.null(labels) && nrow(labels) > 0) {
      labels <- labels %>%
        dplyr::filter(!is.na(Protein), Protein != "")
      plot <- plot + geom_text(
        data = labels,
        aes(label = Protein),
        size = input$volcano_label_size,
        vjust = -0.5,
        check_overlap = TRUE,
        show.legend = FALSE
      )
    }
    plot
  })

  volcano_interactive_obj <- reactive({
    validate(need(requireNamespace("plotly", quietly = TRUE), "Install the R package 'plotly' to show or export the interactive volcano plot."))
    data <- volcano_plot_data() %>%
      dplyr::filter(is.finite(Log2FoldChange), is.finite(MinusLog10Significance))
    metric_label <- if (isTRUE(input$stats_bh_fdr) && input$volcano_significance_metric == "BH_FDR") "BH FDR" else "p-value"
    comparison_label <- stats_comparison_label(input$volcano_comparison)
    count_subtitle <- volcano_count_subtitle(data)
    cutoff_subtitle <- paste0(comparison_label, "; cutoffs: |log2 FC| >= ", input$volcano_fc_cutoff, " and ", metric_label, " <= ", input$volcano_sig_cutoff)
    colors <- c("Increased" = "#B2182B", "Decreased" = "#2166AC", "Not significant" = "#BDBDBD")
    data$Hover <- paste0(
      "Protein: ", data$Protein,
      "<br>Protein name: ", ifelse(is.na(data$ProteinName), "", data$ProteinName),
      "<br>Status: ", data$Status,
      "<br>log2 fold-change: ", signif(data$Log2FoldChange, 4),
      "<br>p-value: ", signif(data$PValue, 4),
      if (isTRUE(input$stats_bh_fdr)) paste0("<br>BH FDR: ", signif(data$BH_FDR, 4)) else ""
    )
    widget <- plotly::plot_ly(
      data,
      x = ~Log2FoldChange,
      y = ~MinusLog10Significance,
      type = "scatter",
      mode = "markers",
      color = ~Status,
      colors = colors,
      text = ~Hover,
      hoverinfo = "text",
      marker = list(size = max(4, input$volcano_point_size * 3), opacity = 0.75)
    )
    widget <- plotly::layout(
      widget,
      title = list(text = paste0(input$volcano_title, "<br><sup>", count_subtitle, "</sup><br><sup>", cutoff_subtitle, "</sup>")),
      xaxis = list(title = paste0("log2 fold-change (", comparison_label, ")")),
      yaxis = list(title = paste0("-log10(", metric_label, ")")),
      legend = list(orientation = "h", x = 0, y = 1.1),
      shapes = list(
        list(type = "line", x0 = -input$volcano_fc_cutoff, x1 = -input$volcano_fc_cutoff, y0 = 0, y1 = 1, yref = "paper", line = list(dash = "dash", color = "#777777")),
        list(type = "line", x0 = input$volcano_fc_cutoff, x1 = input$volcano_fc_cutoff, y0 = 0, y1 = 1, yref = "paper", line = list(dash = "dash", color = "#777777")),
        list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = -log10(input$volcano_sig_cutoff), y1 = -log10(input$volcano_sig_cutoff), line = list(dash = "dash", color = "#777777"))
      )
    )
    widget
  })

  feature_expression_data <- reactive({
    feature_source <- input$feature_data_source
    if (is.null(feature_source) || !nzchar(feature_source)) feature_source <- "imputed"
    expression_data_for_protein_source(feature_source)
  })

  feature_matrix_data <- reactive({
    feature_source <- input$feature_data_source
    if (is.null(feature_source) || !nzchar(feature_source)) feature_source <- "imputed"
    if (identical(feature_source, "S3_batch_corrected")) {
      expr <- expression_data_for_protein_source("S3_batch_corrected")
      values <- as.matrix(expr[, -1, drop = FALSE])
      values[!is.finite(values)] <- NA_real_
      keep_rows <- rowSums(!is.na(values)) > 0
      z_values <- t(apply(values[keep_rows, , drop = FALSE], 1, function(row_values) {
        observed <- is.finite(row_values)
        row_z <- rep(NA_real_, length(row_values))
        row_sd <- stats::sd(row_values[observed], na.rm = TRUE)
        if (sum(observed) == 1 || !is.finite(row_sd) || row_sd == 0) {
          row_z[observed] <- 0
        } else {
          row_z[observed] <- as.numeric(scale(row_values[observed]))
        }
        row_z
      }))
      sample_names <- colnames(expr)[-1]
      if (!is.null(project_file("meta_file"))) {
        md <- active_metadata()
        label_by_sample <- stats::setNames(as.character(md$AnalysisLabel), as.character(md$Sample))
        matched <- sample_names %in% names(label_by_sample)
        sample_names[matched] <- unname(label_by_sample[sample_names[matched]])
      }
      result <- batch_corrected_s3_result()
      info <- result$prepared$feature_info[keep_rows, , drop = FALSE]
      feature_names <- make.unique(as.character(expr$Feature[keep_rows]))
      protein_names <- if ("PG.ProteinNames" %in% colnames(info)) as.character(info[["PG.ProteinNames"]]) else rep(NA_character_, sum(keep_rows))
      protein_descriptions <- if ("PG.ProteinDescriptions" %in% colnames(info)) as.character(info[["PG.ProteinDescriptions"]]) else rep(NA_character_, sum(keep_rows))
      protein_group_id <- if ("PG.ProteinGroups" %in% colnames(info)) as.character(info[["PG.ProteinGroups"]]) else as.character(seq_len(sum(keep_rows)))
      significant <- rep("No", sum(keep_rows))
      volcano_ready <- !is.null(input$volcano_comparison) && identical(input$volcano_source, "S3_batch_corrected")
      if (volcano_ready) {
        volcano_data <- tryCatch(volcano_plot_data(), error = function(e) NULL)
        if (!is.null(volcano_data) && "ProteinGroupID" %in% colnames(volcano_data)) {
          matched <- match(protein_group_id, volcano_data$ProteinGroupID)
          status <- volcano_data$Status[matched]
          significant[status == "Increased" & !is.na(status)] <- "Yes: Increased"
          significant[status == "Decreased" & !is.na(status)] <- "Yes: Decreased"
        }
      }
      log2_table <- as.data.frame(round(values[keep_rows, , drop = FALSE], 1), check.names = FALSE)
      colnames(log2_table) <- make.unique(sample_names)
      display_data <- data.frame(
        Significant = significant,
        Protein = feature_names,
        ProteinName = protein_names,
        ProteinDescription = protein_descriptions,
        log2_table,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      colnames(z_values) <- colnames(log2_table)
      attr(display_data, "z_scores") <- as.data.frame(z_values, check.names = FALSE)
      return(display_data)
    }
    source <- protein_source_table(feature_source)
    quantity_cols <- proteomics_abundance_columns(source)
    validate(need(length(quantity_cols) >= 1, "The selected feature source does not contain abundance columns."))

    feature_col <- resolve_report_feature_col(source, input$report_feature_col)
    validate(need(!is.na(feature_col), "No feature identifier column found in the report."))
    raw_values <- source[, quantity_cols, drop = FALSE]
    raw_values[] <- lapply(raw_values, function(x) suppressWarnings(as.numeric(as.character(x))))
    values <- as.matrix(raw_values)
    values[!is.finite(values) | values <= 0] <- NA_real_
    log2_values <- log2(values)
    keep_rows <- rowSums(!is.na(log2_values)) > 0
    z_values <- t(apply(values[keep_rows, , drop = FALSE], 1, function(row_values) {
      observed <- is.finite(row_values)
      row_z <- rep(NA_real_, length(row_values))
      row_sd <- stats::sd(row_values[observed], na.rm = TRUE)
      if (sum(observed) == 1 || !is.finite(row_sd) || row_sd == 0) {
        row_z[observed] <- 0
      } else {
        row_z[observed] <- as.numeric(scale(row_values[observed]))
      }
      row_z
    }))

    sample_names <- protein_sample_names_for_table(source, quantity_cols)
    if (!is.null(project_db_cache()$metadata) || !is.null(project_file("meta_file"))) {
      md <- active_metadata()
      label_by_sample <- stats::setNames(as.character(md$AnalysisLabel), as.character(md$Sample))
      matched <- sample_names %in% names(label_by_sample)
      sample_names[matched] <- unname(label_by_sample[sample_names[matched]])
    }

    feature_names <- as.character(source[[feature_col]])
    feature_names[is.na(feature_names) | feature_names == ""] <- as.character(source[[1]])[is.na(feature_names) | feature_names == ""]
    feature_names <- make.unique(feature_names[keep_rows])
    protein_names <- if ("PG.ProteinNames" %in% colnames(source)) as.character(source[["PG.ProteinNames"]][keep_rows]) else rep(NA_character_, sum(keep_rows))
    protein_descriptions <- if ("PG.ProteinDescriptions" %in% colnames(source)) as.character(source[["PG.ProteinDescriptions"]][keep_rows]) else rep(NA_character_, sum(keep_rows))
    protein_group_id <- if ("PG.ProteinGroups" %in% colnames(source)) as.character(source[["PG.ProteinGroups"]][keep_rows]) else as.character(which(keep_rows))

    significant <- rep("No", sum(keep_rows))
    volcano_ready <- !is.null(input$volcano_comparison) &&
      ((identical(input$volcano_source, "S2") && protein_source_available("S2")) ||
       (identical(input$volcano_source, "S3") && protein_source_available("S3")))
    if (volcano_ready) {
      volcano_data <- tryCatch(volcano_plot_data(), error = function(e) NULL)
      if (!is.null(volcano_data) && "ProteinGroupID" %in% colnames(volcano_data)) {
        matched <- match(protein_group_id, volcano_data$ProteinGroupID)
        status <- volcano_data$Status[matched]
        significant[status == "Increased" & !is.na(status)] <- "Yes: Increased"
        significant[status == "Decreased" & !is.na(status)] <- "Yes: Decreased"
      }
    }

    log2_table <- as.data.frame(round(log2_values[keep_rows, , drop = FALSE], 1), check.names = FALSE)
    colnames(log2_table) <- make.unique(sample_names)
    display_data <- data.frame(
      Significant = significant,
      Protein = feature_names,
      ProteinName = protein_names,
      ProteinDescription = protein_descriptions,
      log2_table,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    colnames(z_values) <- colnames(log2_table)
    attr(display_data, "z_scores") <- as.data.frame(z_values, check.names = FALSE)
    display_data
  })

  feature_name_for_row <- function(df, row_index) {
    quantity_cols <- proteomics_abundance_columns(df)
    feature_col <- resolve_report_feature_col(df, input$report_feature_col)
    validate(need(!is.na(feature_col), "No feature identifier column found in the report."))
    values <- df[, quantity_cols, drop = FALSE]
    values[] <- lapply(values, function(x) suppressWarnings(as.numeric(as.character(x))))
    keep_rows <- rowSums(!is.na(values)) > 0
    feature_names <- as.character(df[[feature_col]])
    feature_names[is.na(feature_names) | feature_names == ""] <- as.character(df[[1]])[is.na(feature_names) | feature_names == ""]
    unique_names <- rep(NA_character_, nrow(df))
    unique_names[keep_rows] <- make.unique(feature_names[keep_rows])
    unique_names[row_index]
  }

  observeEvent(feature_expression_data(), {
    expr <- feature_expression_data()
    choices <- protein_feature_choices(input$feature_data_source, expr)
    choice_values <- unname(choices)
    current <- isolate(input$feature_select)
    current <- current[current %in% choice_values]
    selected <- if (length(current)) current else choice_values[1]
    freezeReactiveValue(input, "feature_select")
    updateSelectizeInput(session, "feature_select", choices = choices, selected = selected, server = TRUE)
  }, ignoreInit = FALSE)

  observe({
    pending <- pending_feature_selection()
    req(!is.null(pending), identical(input$feature_data_source, pending$source))
    expr <- feature_expression_data()
    choices <- protein_feature_choices(input$feature_data_source, expr)
    selected <- intersect(pending$features, unname(choices))
    pending_feature_selection(NULL)
    req(length(selected))
    session$onFlushed(function() {
      freezeReactiveValue(input, "feature_select")
      updateSelectizeInput(session, "feature_select", choices = choices, selected = selected, server = TRUE)
    }, once = TRUE)
  })

  observe({
    md <- tryCatch(active_metadata(), error = function(e) NULL)
    if (is.null(md)) {
      updateSelectizeInput(session, "feature_order_columns", choices = character(0), selected = character(0), server = FALSE)
      return()
    }
    choices <- setdiff(scoped_metadata_columns(md), c("SampleDetailsID", "Excluded", "ExclusionReason"))
    current <- isolate(input$feature_order_columns)
    selected <- current[current %in% choices]
    saved <- unlist(restored_project_settings()$feature_order_columns, use.names = FALSE)
    saved <- saved[saved %in% choices]
    if (is.null(current) && length(saved) > 0) selected <- saved
    if (is.null(current) && length(selected) == 0) {
      defaults <- c("Condition", "Group", "Batch", "Subject", "Replicate", "RunOrder")
      selected <- defaults[defaults %in% choices]
    }
    updateSelectizeInput(
      session,
      "feature_order_columns",
      choices = choices,
      selected = selected,
      server = FALSE
    )
    current_group <- isolate(input$feature_group_by)
    selected_group <- retain_metadata_choice(current_group, choices, "Condition")
    updateSelectInput(session, "feature_group_by", choices = c("Select grouping field..." = "", stats::setNames(choices, choices)), selected = selected_group)
    label_choices <- c("Select label field..." = "", stats::setNames(choices, choices))
    current_label <- isolate(input$feature_label_by)
    selected_label <- retain_metadata_choice(current_label, choices, c("AnalysisLabel", "Sample"))
    updateSelectInput(session, "feature_label_by", choices = label_choices, selected = selected_label)
  })

  observeEvent(input$protein_no_impute_preview_rows_selected, {
    row_index <- input$protein_no_impute_preview_rows_selected[1]
    req(row_index, protein_source_available("S2"))
    feature <- feature_name_for_row(protein_source_table("S2"), row_index)
    req(!is.na(feature))
    updateRadioButtons(session, "feature_data_source", selected = "no_impute")
    session$onFlushed(function() {
      updateSelectizeInput(session, "feature_select", selected = feature)
      updateTabsetPanel(session, "workflow_tabs", selected = "Feature plot")
    }, once = TRUE)
  })

  observeEvent(input$protein_imputed_preview_rows_selected, {
    row_index <- input$protein_imputed_preview_rows_selected[1]
    req(row_index, protein_source_available("S3"))
    feature <- feature_name_for_row(protein_source_table("S3"), row_index)
    req(!is.na(feature))
    updateRadioButtons(session, "feature_data_source", selected = "imputed")
    session$onFlushed(function() {
      updateSelectizeInput(session, "feature_select", selected = feature)
      updateTabsetPanel(session, "workflow_tabs", selected = "Feature plot")
    }, once = TRUE)
  })

  observeEvent(input$feature_matrix_table_rows_selected, {
    row_index <- input$feature_matrix_table_rows_selected[1]
    req(row_index)
    table <- feature_matrix_data()
    req(row_index <= nrow(table))
    updateRadioButtons(session, "script_box_source", selected = input$feature_data_source)
    updateSelectizeInput(session, "script_box_features", selected = table$Protein[row_index])
    updateTabsetPanel(session, "workflow_tabs", selected = "Boxplot")
  })

  cv_from_report <- function(report) {
    md <- active_metadata()
    group_col <- input$cv_plot_group_col
    validate(need(length(group_col) == 1L && group_col %in% colnames(md), "Select a metadata field for CV grouping."))

    abundance_columns <- proteomics_abundance_columns(report)
    run_labels <- proteomics_abundance_sample_names(abundance_columns)
    value_cols <- match(abundance_columns, colnames(report))
    restored_processed <- any(endsWith(abundance_columns, "_Protein_group_abundance"))
    header_labels <- if (restored_processed) protein_header_labels_from_metadata(md, input$protein_header_label_columns) else as.character(md$Sample)
    group_values <- as.character(md[[group_col]])
    condition_data <- unique(group_values[!is.na(group_values) & group_values != ""])
    density_rows <- list()
    median_rows <- list()

    for (condition in condition_data) {
      condition_runs <- header_labels[group_values == condition]
      cols <- value_cols[run_labels %in% condition_runs]
      if (length(cols) < 2) next
      values <- as.data.frame(report[, cols, drop = FALSE], stringsAsFactors = FALSE)
      values[] <- lapply(values, function(column) suppressWarnings(as.numeric(as.character(column))))
      cv <- apply(values, 1, function(row_values) {
        row_values <- row_values[is.finite(row_values)]
        mean_value <- mean(row_values)
        if (length(row_values) < 2 || !is.finite(mean_value) || mean_value <= 0) return(NA_real_)
        100 * stats::sd(row_values) / mean_value
      })
      cv <- cv[is.finite(cv)]
      if (length(cv) < 2) next
      curve <- stats::density(cv, from = 0, n = 1024, adjust = input$cv_density_adjust)
      density_rows[[condition]] <- data.frame(Condition = condition, CV = curve$x, Density = curve$y, stringsAsFactors = FALSE)
      median_rows[[condition]] <- data.frame(Condition = condition, MedianCV = stats::median(cv), stringsAsFactors = FALSE)
    }

    validate(need(length(density_rows) > 0, "Not enough abundance values to calculate CV distributions."))
      list(
        density = dplyr::bind_rows(density_rows),
        medians = dplyr::bind_rows(median_rows),
        mode = paste0("Calculated from protein abundance values, grouped by ", group_col, " (density smoothness ", input$cv_density_adjust, ")")
      )
  }

  cv_from_batch_corrected_s3 <- function() {
    result <- batch_corrected_s3_result()
    expression <- result$expression
    sample_map <- result$prepared$sample_map
    active_md <- active_metadata()
    group_col <- input$cv_plot_group_col
    validate(need(length(group_col) == 1L && group_col %in% colnames(active_md), "Select a metadata field for CV grouping."))
    if ("Sample" %in% colnames(active_md)) {
      sample_map <- sample_map[sample_map$Sample %in% as.character(active_md$Sample), , drop = FALSE]
    }
    validate(
      need(ncol(expression) > 2, "Run batch correction before calculating CV distributions from batch-corrected S3."),
      need("Sample" %in% colnames(sample_map), "Batch-corrected sample metadata is missing sample labels.")
    )

    value_cols <- intersect(sample_map$Sample, colnames(expression))
    validate(need(length(value_cols) >= 2, "Batch-corrected S3 has fewer than two corrected sample columns."))
    values <- as.data.frame(expression[, value_cols, drop = FALSE], stringsAsFactors = FALSE)
    values[] <- lapply(values, function(column) 2^suppressWarnings(as.numeric(as.character(column))))

    condition_by_sample <- cv_sample_groups(active_md, group_col, as.character(sample_map$Sample))
    condition_data <- unique(condition_by_sample[value_cols])
    condition_data <- condition_data[!is.na(condition_data) & nzchar(condition_data)]
    density_rows <- list()
    median_rows <- list()

    for (condition in condition_data) {
      cols <- value_cols[condition_by_sample[value_cols] == condition]
      if (length(cols) < 2) next
      cv <- apply(values[, cols, drop = FALSE], 1, function(row_values) {
        row_values <- row_values[is.finite(row_values)]
        mean_value <- mean(row_values)
        if (length(row_values) < 2 || !is.finite(mean_value) || mean_value <= 0) return(NA_real_)
        100 * stats::sd(row_values) / mean_value
      })
      cv <- cv[is.finite(cv)]
      if (length(cv) < 2) next
      curve <- stats::density(cv, from = 0, n = 1024, adjust = input$cv_density_adjust)
      density_rows[[condition]] <- data.frame(Condition = condition, CV = curve$x, Density = curve$y, stringsAsFactors = FALSE)
      median_rows[[condition]] <- data.frame(Condition = condition, MedianCV = stats::median(cv), stringsAsFactors = FALSE)
    }

    validate(need(length(density_rows) > 0, "Not enough batch-corrected abundance values to calculate CV distributions. Each plotted condition needs at least two corrected samples."))
    list(
      density = dplyr::bind_rows(density_rows),
      medians = dplyr::bind_rows(median_rows),
      mode = paste0("Calculated from batch-corrected S3 log2 values converted back to abundance scale, grouped by ", group_col, " (density smoothness ", input$cv_density_adjust, ")")
    )
  }

  cv_from_distribution_file <- reactive({
    req(project_file("cv_distribution_file"))
    source <- read_uploaded_table(project_file("cv_distribution_file"))
    x_cols <- grep("^\\(x\\)% CV_\\[", colnames(source), value = TRUE)
    validate(need(length(x_cols) > 0, "CV distribution table must contain columns such as '(x)% CV_[CKD]'."))
    density_rows <- lapply(x_cols, function(x_col) {
      condition <- sub("^\\(x\\)% CV_\\[([^]]+)\\]$", "\\1", x_col)
      y_col <- paste0("(y)Density_[", condition, "]")
      validate(need(y_col %in% colnames(source), paste0("Missing density column: ", y_col)))
      data.frame(
        Condition = condition,
        CV = suppressWarnings(as.numeric(source[[x_col]])),
        Density = suppressWarnings(as.numeric(source[[y_col]])),
        stringsAsFactors = FALSE
      )
    })
    density_df <- dplyr::bind_rows(density_rows) %>% dplyr::filter(is.finite(CV), is.finite(Density))
    median_rows <- density_df %>%
      dplyr::group_by(Condition) %>%
      dplyr::summarise(
        MedianCV = {
          order_index <- order(CV)
          x <- CV[order_index]
          y <- Density[order_index]
          cumulative <- cumsum(c(0, diff(x) * (head(y, -1) + tail(y, -1)) / 2))
          cumulative <- cumulative / max(cumulative, na.rm = TRUE)
          x[which.min(abs(cumulative - 0.5))]
        },
        .groups = "drop"
      )
    list(density = density_df, medians = median_rows, mode = "Uploaded density distribution (median approximated from curve)")
  })

  cv_plot_data <- reactive({
    cv_source <- input$cv_plot_source
    if (is.null(cv_source) || !nzchar(cv_source)) cv_source <- "uploaded"
    switch(
      cv_source,
      "uploaded" = cv_from_distribution_file(),
      "no_impute" = {
        cv_from_report(protein_source_table("no_impute"))
      },
      "imputed" = {
        cv_from_report(protein_source_table("imputed"))
      },
      "S3_batch_corrected" = {
        cv_from_batch_corrected_s3()
      }
    )
  })

  observeEvent(cv_plot_data(), {
    conditions <- unique(as.character(cv_plot_data()$density$Condition))
    current <- isolate(input$cv_plot_conditions)
    selected <- current[current %in% conditions]
    if (length(selected) == 0) selected <- conditions
    updateSelectizeInput(session, "cv_plot_conditions", choices = conditions, selected = selected, server = TRUE)
  }, ignoreInit = FALSE)

  cv_plot_obj <- reactive({
    plot_data <- cv_plot_data()
    selected <- input$cv_plot_conditions
    if (length(selected) == 0) selected <- unique(plot_data$density$Condition)
    density_df <- plot_data$density %>% dplyr::filter(Condition %in% selected)
    medians <- plot_data$medians %>% dplyr::filter(Condition %in% selected)
    validate(need(nrow(density_df) > 0, "Select at least one condition for the CV plot."))

    colors <- c("CKD" = "#00B050", "CKDu" = "#6D8CFF", "Control" = "#FF3333", "SPQC" = "#FF00FF")
    extra <- setdiff(unique(density_df$Condition), names(colors))
    if (length(extra) > 0) colors <- c(colors, stats::setNames(scales::hue_pal()(length(extra)), extra))
    medians$Label <- paste0("Median: ", format(round(medians$MedianCV, 1), nsmall = 1), "%")
    max_y <- max(density_df$Density, na.rm = TRUE)
    medians$LabelY <- max_y * seq(0.86, 0.62, length.out = nrow(medians))

    p <- ggplot(density_df, aes(x = CV, y = Density, color = Condition, fill = Condition, group = Condition))
    if (isTRUE(input$cv_fill_density)) {
      p <- p + geom_area(alpha = 0.08, position = "identity", color = NA)
    }
    p +
      geom_line(linewidth = input$cv_line_width) +
      geom_vline(data = medians, aes(xintercept = MedianCV, color = Condition), linetype = "dashed", linewidth = 0.5) +
      geom_text(
        data = medians,
        aes(x = MedianCV, y = LabelY, label = Label, color = Condition),
        hjust = -0.02,
        size = input$cv_median_text_size,
        show.legend = FALSE
      ) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      coord_cartesian(xlim = c(0, input$cv_x_cutoff), expand = FALSE) +
      labs(title = input$cv_plot_title, x = "% CV", y = "Density", color = NULL, fill = NULL) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = input$cv_title_size),
        axis.title = element_text(size = input$cv_axis_title_size),
        axis.text = element_text(size = input$cv_axis_text_size),
        legend.text = element_text(size = input$cv_legend_text_size),
        legend.position = "top",
        legend.justification = "left",
        legend.box.just = "left",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  })

  identifications_overview_data <- reactive({
    data <- project_input_table(
      project_file("identifications_overview_file"),
      project_db_cache()$identifications_overview
    )
    validate(need(!is.null(data), "Select an IdentificationsOverview table or open a project containing one."))
    required <- c("Condition", "Replicate", "Precursors", "ProteinGroups")
    validate(need(all(required %in% colnames(data)), "IdentificationsOverview table must contain Condition, Replicate, Precursors, and ProteinGroups columns."))
    data$SourceLabel <- paste(data$Condition, data$Replicate, sep = ".")
    data$RunLabel <- data$SourceLabel
    if (!is.null(project_file("meta_file"))) {
      md <- active_metadata()
      validate(need(all(c("Condition", "Replicate") %in% colnames(md)), "Metadata must contain Condition and Replicate to label identification plots."))
      metadata_labels <- md[, c("Condition", "Replicate"), drop = FALSE]
      metadata_labels$Condition <- as.character(metadata_labels$Condition)
      metadata_labels$Replicate <- as.character(metadata_labels$Replicate)
      metadata_labels <- metadata_labels[!duplicated(metadata_labels[, c("Condition", "Replicate")]), , drop = FALSE]
      metadata_labels$.metadata_match <- TRUE
      data <- data %>%
        dplyr::mutate(Condition = as.character(Condition), Replicate = as.character(Replicate)) %>%
        dplyr::left_join(metadata_labels, by = c("Condition", "Replicate"))
      matched <- !is.na(data$.metadata_match)
      preferred_labels <- condition_replicate_label(data$Condition, data$Replicate, data$RunLabel)
      data$RunLabel[matched] <- preferred_labels[matched]
      data$.metadata_match <- NULL
    }
    data$RunLabel <- factor(data$RunLabel, levels = data$RunLabel)
    data
  })

  identification_overview_plot_obj <- reactive({
    data <- identifications_overview_data()
    metric <- input$identification_metric
    validate(need(metric %in% colnames(data), paste0("Overview metric not found: ", metric)))
    data$Value <- suppressWarnings(as.numeric(data[[metric]])) / 1000
    metric_label <- if (metric == "ProteinGroups") "Protein Groups" else "Precursors"
    ggplot(data, aes(x = RunLabel, y = Value, fill = Condition)) +
      geom_col(width = 0.72, color = "black", linewidth = 0.25) +
      labs(
        title = input$identification_overview_title,
        x = "Run",
        y = paste0("Nr. of ", metric_label, " (10^3)"),
        fill = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = input$identification_title_size),
        axis.text.x = element_text(angle = 50, hjust = 1, size = input$identification_axis_text_size),
        axis.text.y = element_text(size = input$identification_axis_text_size),
        legend.text = element_text(size = input$identification_legend_size),
        legend.position = "top",
        panel.grid.minor = element_blank()
      )
  })

  run_identifications_data <- reactive({
    data <- if (input$identification_metric == "ProteinGroups") {
      project_input_table(
        project_file("run_identifications_protein_file"),
        project_db_cache()$run_identifications_protein
      )
    } else {
      project_input_table(
        project_file("run_identifications_precursor_file"),
        project_db_cache()$run_identifications_precursor
      )
    }
    validate(need(!is.null(data), "Select the corresponding run-identifications table or open a project containing one."))
    validate(need(ncol(data) >= 5, "Run identifications table must contain a run-label column and four identification category columns."))
    names(data)[1:5] <- c("SourceLabel", "Complete Identifications", "Shared in >=50% of the Runs", "Sparse Identifications", "Unique Identifications")
    overview <- identifications_overview_data()
    label_by_source <- stats::setNames(as.character(overview$RunLabel), as.character(overview$SourceLabel))
    run_labels <- unname(label_by_source[as.character(data$SourceLabel)])
    missing_labels <- is.na(run_labels) | run_labels == ""
    run_labels[missing_labels] <- as.character(data$SourceLabel)[missing_labels]
    category_names <- names(data)[2:5]
    long <- data.frame(
      RunLabel = rep(run_labels, times = length(category_names)),
      Category = rep(category_names, each = nrow(data)),
      Value = as.numeric(unlist(data[, category_names, drop = FALSE])),
      stringsAsFactors = FALSE
    )
    long$RunLabel <- factor(long$RunLabel, levels = as.character(overview$RunLabel))
    long$Category <- factor(
      long$Category,
      levels = c("Unique Identifications", "Sparse Identifications", "Shared in >=50% of the Runs", "Complete Identifications")
    )
    long
  })

  run_identifications_plot_obj <- reactive({
    data <- run_identifications_data()
    metric_label <- if (input$identification_metric == "ProteinGroups") "Protein Groups" else "Precursors"
    fill_colors <- c(
      "Complete Identifications" = "#5E5E5E",
      "Shared in >=50% of the Runs" = "#A6A6A6",
      "Sparse Identifications" = "#FFFFFF",
      "Unique Identifications" = "#C8102E"
    )
    ggplot(data, aes(x = RunLabel, y = Value / 1000, fill = Category)) +
      geom_col(width = 0.5, color = "black", linewidth = 0.35) +
      scale_fill_manual(
        values = fill_colors,
        breaks = c("Complete Identifications", "Shared in >=50% of the Runs", "Sparse Identifications", "Unique Identifications"),
        drop = FALSE
      ) +
      labs(
        title = input$run_identifications_title,
        x = "Run",
        y = paste0("Nr. of ", metric_label, " (10^3)"),
        fill = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = input$identification_title_size),
        axis.text.x = element_text(angle = 50, hjust = 1, size = input$identification_axis_text_size),
        axis.text.y = element_text(size = input$identification_axis_text_size),
        legend.text = element_text(size = input$identification_legend_size),
        legend.position = "top",
        legend.justification = "left",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      )
  })

  make_display_label <- function(df) {
    out <- rep(NA_character_, nrow(df))
    has_group <- "Group" %in% colnames(df)
    has_subject <- "Subject" %in% colnames(df)
    has_timepoint <- "Timepoint" %in% colnames(df)

    if (has_group && has_subject && has_timepoint) {
      is_spqc <- !is.na(df$Group) & grepl("spqc", df$Group, ignore.case = TRUE)
      out[!is_spqc] <- paste(df$Group[!is_spqc], df$Subject[!is_spqc], df$Timepoint[!is_spqc], sep = "_")

      if (any(is_spqc)) {
        samp <- df$Sample[is_spqc]

        # extract replicate numbers from sample names while preserving input order
        rep_num <- stringr::str_match(samp, "SPQC[_-]?(\\d+)")[, 2]

        # fallback: use current order if replicate number is not present
        missing_rep <- is.na(rep_num) | rep_num == ""
        if (any(missing_rep)) {
          rep_num[missing_rep] <- as.character(seq_along(rep_num)[missing_rep])
        }

        out[is_spqc] <- paste0("SPQC_", rep_num)
      }
    }

    missing_out <- is.na(out) | out == "" | grepl("^NA", out)
    if (any(missing_out) && has_subject && has_timepoint) {
      out[missing_out] <- paste(df$Subject[missing_out], df$Timepoint[missing_out], sep = "_")
    }

    out[is.na(out) | out == ""] <- df$Sample[is.na(out) | out == ""]
    out
  }

  observe({
    md <- tryCatch(built_metadata(), error = function(e) NULL)
    if (!is.null(md)) {
      meta_choices <- setdiff(scoped_metadata_columns(md), "Sample")
      label_choices <- scoped_metadata_columns(md)
      current_label <- isolate(input$label_mode)
      selected_label <- if (!is.null(current_label) && current_label %in% label_choices) {
        current_label
      } else if ("AnalysisLabel" %in% label_choices) {
        "AnalysisLabel"
      } else {
        label_choices[1]
      }
      selected_label <- retain_metadata_choice(current_label, label_choices, "AnalysisLabel")
      updateSelectInput(session, "label_mode", choices = c("Select label field..." = "", stats::setNames(label_choices, label_choices)), selected = selected_label)
      updateSelectInput(session, "color_by",
        choices = c("Select metadata field..." = "", "None" = "None", stats::setNames(meta_choices, meta_choices)),
        selected = retain_metadata_choice(isolate(input$color_by), c("None", meta_choices), c("Condition", "None")))
      updateSelectInput(session, "shape_by",
        choices = c("Select metadata field..." = "", "None" = "None", stats::setNames(meta_choices, meta_choices)),
        selected = retain_metadata_choice(isolate(input$shape_by), c("None", meta_choices), "None"))
    } else {
      updateSelectInput(session, "label_mode", choices = c("Sample"), selected = "Sample")
      updateSelectInput(session, "color_by", choices = "None", selected = "None")
      updateSelectInput(session, "shape_by", choices = "None", selected = "None")
    }
  })

  observe({
    md <- tryCatch(built_metadata(), error = function(e) NULL)
    choices <- if (!is.null(md)) scoped_metadata_columns(md) else character(0)
    metadata_choices <- c("Select metadata field..." = "", "None" = "None", stats::setNames(choices, choices))
    selected_color <- retain_metadata_choice(isolate(input$clustvis_pca_color_by), c("None", choices), c("Condition", "None"))
    selected_shape <- retain_metadata_choice(isolate(input$clustvis_pca_shape_by), c("None", choices), c("Replicate", "None"))
    selected_label <- retain_metadata_choice(isolate(input$clustvis_pca_label_by), c("None", choices), "None")
    selected_subset <- retain_metadata_choice(isolate(input$clustvis_pca_subset_column), c("", choices), "")
    updateSelectInput(session, "clustvis_pca_color_by", choices = metadata_choices, selected = selected_color)
    updateSelectInput(session, "clustvis_pca_shape_by", choices = metadata_choices, selected = selected_shape)
    updateSelectInput(session, "clustvis_pca_label_by", choices = metadata_choices, selected = selected_label)
    updateSelectInput(
      session,
      "clustvis_pca_subset_column",
      choices = c("All samples" = "", stats::setNames(choices, choices)),
      selected = selected_subset
    )
  })

  observe({
    md <- tryCatch(active_metadata(), error = function(e) NULL)
    subset_column <- input$clustvis_pca_subset_column
    if (is.null(md) || is.null(subset_column) || !nzchar(subset_column) || !subset_column %in% colnames(md)) {
      updateSelectizeInput(session, "clustvis_pca_subset_values", choices = character(0), selected = character(0), server = TRUE)
      return()
    }
    values <- sort(unique(trimws(as.character(md[[subset_column]]))))
    values <- values[!is.na(values) & nzchar(values)]
    current <- isolate(input$clustvis_pca_subset_values)
    selected <- current[current %in% values]
    updateSelectizeInput(
      session,
      "clustvis_pca_subset_values",
      choices = stats::setNames(values, values),
      selected = selected,
      server = TRUE
    )
  })

  observe({
    scores <- tryCatch(clustvis_pca_results()$scores, error = function(e) NULL)
    color_by <- input$clustvis_pca_color_by
    if (is.null(scores) || is.null(color_by) || identical(color_by, "None") || !color_by %in% colnames(scores)) {
      updateSelectizeInput(session, "clustvis_pca_opacity_override_groups", choices = character(0), selected = character(0), server = TRUE)
      return()
    }
    groups <- sort(unique(as.character(scores[[color_by]])))
    groups <- groups[!is.na(groups) & nzchar(groups)]
    current <- isolate(input$clustvis_pca_opacity_override_groups)
    selected <- current[current %in% groups]
    updateSelectizeInput(
      session,
      "clustvis_pca_opacity_override_groups",
      choices = stats::setNames(groups, groups),
      selected = selected,
      server = TRUE
    )
  })

  clustvis_pca_results <- eventReactive(list(input$run_clustvis_pca, metadata_apply_revision()), {
    req(input$run_clustvis_pca > 0)
    expr <- expression_data_for_protein_source(input$clustvis_pca_source)
    X <- t(as.matrix(expr[, -1, drop = FALSE]))
    X <- normalize_pca_missing_values(X)
    rownames(X) <- colnames(expr)[-1]
    subset_result <- subset_pca_samples(
      X,
      sample_metadata_for_plotting(rownames(X)),
      input$clustvis_pca_subset_column,
      input$clustvis_pca_subset_values
    )
    X <- subset_result$expression
    min_observed <- ceiling(nrow(X) * input$clustvis_pca_min_observed_percent / 100)
    keep_features <- colSums(!is.na(X)) >= min_observed
    retained_features <- sum(keep_features)
    total_features <- ncol(X)
    X <- X[, keep_features, drop = FALSE]

    validate(
      need(nrow(X) >= 3, "Need at least 3 samples for ClustVis PCA."),
      need(ncol(X) >= 2, "Fewer than 2 proteins meet the minimum observed-samples filter. Lower the threshold.")
    )

    rank_use <- min(max(2, as.integer(input$clustvis_pca_npcs)), min(dim(X)) - 1)
    imputation_method <- "Median fallback imputation"
    X_imp <- NULL
    if (requireNamespace("pcaMethods", quietly = TRUE)) {
      pca_model <- tryCatch(
        suppressWarnings(pcaMethods::pca(X, method = "svdImpute", nPcs = rank_use, scale = "none")),
        error = function(e) NULL
      )
      if (!is.null(pca_model)) {
        X_imp <- pca_model@completeObs
        imputation_method <- paste0("pcaMethods::svdImpute, nPcs = ", rank_use)
      }
    }
    if (is.null(X_imp)) {
      X_imp <- X
      for (column_index in seq_len(ncol(X_imp))) {
        values <- X_imp[, column_index]
        replacement <- stats::median(values, na.rm = TRUE)
        if (!is.finite(replacement)) replacement <- 0
        values[is.na(values)] <- replacement
        X_imp[, column_index] <- values
      }
    }
    X_imp[!is.finite(X_imp)] <- NA_real_
    dimnames(X_imp) <- dimnames(X)
    X_imp[X_imp < 0] <- 0.00001
    validate(need(!anyNA(X_imp), "ClustVis PCA matrix still contains missing values after imputation."))

    variable_features <- apply(X_imp, 2, function(values) {
      values <- values[is.finite(values)]
      length(values) >= 2 && stats::sd(values) > 0
    })
    X_imp <- X_imp[, variable_features, drop = FALSE]
    validate(need(ncol(X_imp) >= 2, "Fewer than 2 variable proteins remain after imputation."))

    pca <- stats::prcomp(X_imp, scale. = isTRUE(input$clustvis_pca_scale), rank. = min(rank_use, ncol(X_imp), nrow(X_imp) - 1))
    variance <- (pca$sdev^2) / sum(pca$sdev^2) * 100
    scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    colnames(scores) <- c("PC1", "PC2")
    scores$Sample <- rownames(scores)
    scores <- scores %>% left_join(sample_metadata_for_plotting(scores$Sample), by = "Sample")

    list(
      scores = scores,
      pca = pca,
      protein_info = protein_info_for_source(input$clustvis_pca_source, input$report_feature_col),
      variance = variance,
      retained_features = retained_features,
      variable_features = ncol(X_imp),
      total_features = total_features,
      included_samples = nrow(X),
      excluded_by_subset = subset_result$excluded_samples,
      missing_values_before_imputation = sum(is.na(X)),
      imputation_method = imputation_method
    )
  })

  clustvis_pca_plot_obj <- reactive({
    res <- clustvis_pca_results()
    df <- res$scores
    color_var <- if (!is.null(input$clustvis_pca_color_by) && input$clustvis_pca_color_by != "None" && input$clustvis_pca_color_by %in% colnames(df)) input$clustvis_pca_color_by else NULL
    shape_var <- if (!is.null(input$clustvis_pca_shape_by) && input$clustvis_pca_shape_by != "None" && input$clustvis_pca_shape_by %in% colnames(df)) input$clustvis_pca_shape_by else NULL
    label_var <- if (!is.null(input$clustvis_pca_label_by) && input$clustvis_pca_label_by != "None" && input$clustvis_pca_label_by %in% colnames(df)) input$clustvis_pca_label_by else NULL
    shape_plot_var <- NULL
    if (!is.null(shape_var)) {
      shape_plot_var <- ".ShapeValue"
      df[[shape_plot_var]] <- factor(as.character(df[[shape_var]]))
    }
    default_alpha <- suppressWarnings(as.numeric(input$clustvis_pca_default_opacity))
    override_alpha <- suppressWarnings(as.numeric(input$clustvis_pca_override_opacity))
    if (!is.finite(default_alpha)) default_alpha <- 0.65
    if (!is.finite(override_alpha)) override_alpha <- 0.95
    default_alpha <- min(max(default_alpha, 0), 1)
    override_alpha <- min(max(override_alpha, 0), 1)
    override_groups <- input$clustvis_pca_opacity_override_groups
    if (is.null(override_groups)) override_groups <- character(0)
    override_groups <- as.character(override_groups)
    if (!is.null(color_var)) {
      color_levels <- sort(unique(as.character(df[[color_var]])))
      color_levels <- color_levels[!is.na(color_levels) & nzchar(color_levels)]
      override_groups <- intersect(override_groups, color_levels)
      df$.PCAColorGroup <- factor(as.character(df[[color_var]]), levels = color_levels)
      df$.PCAAlpha <- ifelse(as.character(df[[color_var]]) %in% override_groups, override_alpha, default_alpha)
    } else {
      df$.PCAAlpha <- default_alpha
    }
    hover_parts <- list(
      paste0("Sample: ", df$Sample),
      paste0("PC1: ", signif(df$PC1, 4)),
      paste0("PC2: ", signif(df$PC2, 4))
    )
    if (!is.null(color_var)) hover_parts <- c(hover_parts, list(paste0(color_var, ": ", df[[color_var]])))
    if (!is.null(shape_var) && !identical(shape_var, color_var)) hover_parts <- c(hover_parts, list(paste0(shape_var, ": ", df[[shape_var]])))
    if (!is.null(label_var) && !label_var %in% c(color_var, shape_var, "Sample")) hover_parts <- c(hover_parts, list(paste0(label_var, ": ", df[[label_var]])))
    df$.PCAHover <- do.call(paste, c(hover_parts, sep = "<br>"))

    p <- ggplot(df, aes(x = PC1, y = PC2, text = .PCAHover))
    if (isTRUE(input$clustvis_pca_ellipses) && !is.null(color_var)) {
      p <- p + stat_ellipse(aes(color = .PCAColorGroup, group = .PCAColorGroup), type = "norm", linewidth = 0.6, show.legend = FALSE)
    }
    no_fill_shapes <- isTRUE(input$clustvis_pca_no_fill_shapes)
    point_shape <- if (no_fill_shapes) 1 else 19
    if (!is.null(color_var) && !is.null(shape_plot_var)) {
      p <- p + geom_point(aes(color = .PCAColorGroup, shape = .data[[shape_plot_var]], alpha = .PCAAlpha), size = input$clustvis_pca_point_size, show.legend = TRUE)
    } else if (!is.null(color_var)) {
      p <- p + geom_point(aes(color = .PCAColorGroup, alpha = .PCAAlpha), shape = point_shape, size = input$clustvis_pca_point_size, show.legend = TRUE)
    } else if (!is.null(shape_plot_var)) {
      p <- p + geom_point(aes(shape = .data[[shape_plot_var]], alpha = .PCAAlpha), size = input$clustvis_pca_point_size, show.legend = TRUE)
    } else {
      p <- p + geom_point(aes(alpha = .PCAAlpha), shape = point_shape, size = input$clustvis_pca_point_size)
    }
    p <- p + scale_alpha_identity()
    if (!is.null(shape_plot_var)) {
      shape_values <- if (no_fill_shapes) {
        rep(c(1, 2, 0, 5, 6, 3, 7, 8, 9, 10, 11, 12, 13, 14), length.out = nlevels(df[[shape_plot_var]]))
      } else {
        rep(c(16, 17, 15, 3, 7, 8, 0, 1, 2, 5, 6, 9, 10, 11, 12, 13, 14), length.out = nlevels(df[[shape_plot_var]]))
      }
      p <- p + scale_shape_manual(values = shape_values)
    }
    if (!is.null(label_var)) {
      label_df <- df
      label_df$Label <- as.character(label_df[[label_var]])
      label_df <- label_df[!is.na(label_df$Label) & label_df$Label != "", , drop = FALSE]
      p <- p + geom_text(data = label_df, aes(label = Label), vjust = -0.75, size = input$clustvis_pca_label_size, show.legend = FALSE)
    }
    p +
      theme_bw(base_size = 12) +
      guides(
        color = guide_legend(order = 1, override.aes = list(alpha = 1)),
        shape = guide_legend(order = 2, override.aes = list(alpha = 1))
      ) +
      labs(
        title = input$clustvis_pca_title,
        x = paste0("PC1 (", round(res$variance[1], 2), "%)"),
        y = paste0("PC2 (", round(res$variance[2], 2), "%)"),
        color = color_var,
        shape = shape_var
      ) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold")
      )
  })

  split_plotly_pca_legend <- function(plot_obj, color_var, shape_var) {
    if (is.null(color_var) || is.null(shape_var) || !requireNamespace("plotly", quietly = TRUE)) return(plot_obj)
    traces <- plot_obj$x$data
    color_map <- list()
    shape_map <- list()
    for (i in seq_along(traces)) {
      trace_name <- traces[[i]]$name
      if (is.null(trace_name) || !nzchar(trace_name) || !grepl(",", trace_name, fixed = TRUE)) next
      parts <- trimws(strsplit(trace_name, ",", fixed = TRUE)[[1]])
      if (length(parts) < 2) next
      color_level <- parts[[1]]
      shape_level <- parts[[2]]
      marker <- traces[[i]]$marker
      if (!is.null(marker$color) && is.null(color_map[[color_level]])) color_map[[color_level]] <- marker$color
      if (!is.null(marker$symbol) && is.null(shape_map[[shape_level]])) shape_map[[shape_level]] <- marker$symbol
    }
    if (length(color_map) == 0 && length(shape_map) == 0) return(plot_obj)
    for (i in seq_along(traces)) {
      trace_name <- traces[[i]]$name
      if (is.null(trace_name) || !nzchar(trace_name) || !grepl(",", trace_name, fixed = TRUE)) next
      plot_obj$x$data[[i]]$showlegend <- FALSE
    }
    if (length(color_map) > 0) {
      for (level in names(color_map)) {
        plot_obj <- plotly::add_trace(
          plot_obj,
          x = NA_real_, y = NA_real_,
          type = "scatter", mode = "markers",
          name = paste0(color_var, ": ", level),
          marker = list(color = color_map[[level]], symbol = "circle", size = 10),
          hoverinfo = "skip",
          showlegend = TRUE,
          inherit = FALSE
        )
      }
    }
    if (length(shape_map) > 0) {
      for (level in names(shape_map)) {
        plot_obj <- plotly::add_trace(
          plot_obj,
          x = NA_real_, y = NA_real_,
          type = "scatter", mode = "markers",
          name = paste0(shape_var, ": ", level),
          marker = list(color = "black", symbol = shape_map[[level]], size = 10),
          hoverinfo = "skip",
          showlegend = TRUE,
          inherit = FALSE
        )
      }
    }
    plotly::layout(plot_obj, showlegend = TRUE)
  }

  pca_loadings_data <- reactive({
    res <- clustvis_pca_results()
    rank_pca_loadings(
      pca_model = res$pca,
      protein_info = res$protein_info,
      pc_mode = input$pca_loading_rank_by,
      top_n = input$pca_loading_top_n
    )
  })

  pca_loadings_plot_obj <- reactive({
    data <- pca_loadings_data()
    validate(need(nrow(data) > 0, "Run PCA to calculate PC1/PC2 protein loadings."))
    plot_data <- data
    plot_data$Label <- best_protein_display_label(plot_data, fallback_col = "Protein")
    plot_data$Contribution <- switch(
      input$pca_loading_rank_by,
      "PC1 positive" = plot_data$PC1_Loading,
      "PC1 negative" = plot_data$PC1_Loading,
      "PC2 positive" = plot_data$PC2_Loading,
      "PC2 negative" = plot_data$PC2_Loading,
      plot_data$Combined_PC1_PC2
    )
    plot_data <- head(plot_data, min(25, nrow(plot_data)))
    plot_data$Label <- factor(plot_data$Label, levels = rev(plot_data$Label))
    ggplot(plot_data, aes(x = Label, y = Contribution, fill = Contribution >= 0)) +
      geom_col(width = 0.75) +
      coord_flip() +
      scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"), guide = "none") +
      labs(
        title = "Top PCA Protein Contributors",
        subtitle = input$pca_loading_rank_by,
        x = NULL,
        y = if (identical(input$pca_loading_rank_by, "combined")) "Combined PC1/PC2 loading" else "Loading"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major.y = element_blank()
      )
  })

  pca_results <- eventReactive(input$run_pca, {
    df <- parsed_data()
    X <- as.matrix(df[, -1, drop = FALSE])
    X <- normalize_pca_missing_values(X)
    rownames(X) <- df$Sample
    min_observed <- ceiling(nrow(X) * input$pca_min_observed_percent / 100)
    keep_features <- colSums(!is.na(X)) >= min_observed
    retained_features <- sum(keep_features)
    total_features <- ncol(X)
    X <- X[, keep_features, drop = FALSE]

    validate(
      need(retained_features >= 2, "Fewer than 2 proteins meet the minimum observed-samples filter. Lower the PCA threshold."),
      need(!any(rowSums(!is.na(X)) == 0), "At least one sample has no retained protein values after filtering.")
    )

    ncp_limit <- max(1, min(dim(X)) - 1)
    max_ncp <- min(input$max_ncp, ncp_limit)
    manual_ncp <- suppressWarnings(as.integer(input$manual_ncp))
    if (length(manual_ncp) == 0 || is.na(manual_ncp) || !is.finite(manual_ncp)) manual_ncp <- min(2, ncp_limit)
    fallback_ncp <- min(max(1, manual_ncp), ncp_limit)
    if (isTRUE(input$use_estimated_ncp)) {
      est <- tryCatch(
        suppressWarnings(missMDA::estim_ncpPCA(X, method = "Regularized", ncp.max = max_ncp)),
        error = function(e) NULL
      )
      ncp_use <- if (is.list(est) && !is.null(est$ncp)) est$ncp else suppressWarnings(as.numeric(est)[1])
      ncp_use <- suppressWarnings(as.numeric(ncp_use)[1])
      if (is.null(ncp_use) || length(ncp_use) == 0 || !is.finite(ncp_use) || is.na(ncp_use) || ncp_use < 1) {
        ncp_use <- fallback_ncp
      }
    } else {
      ncp_use <- fallback_ncp
    }
    ncp_use <- min(max(1, as.integer(round(ncp_use))), ncp_limit)

    imp <- tryCatch(
      suppressWarnings(missMDA::imputePCA(X, ncp = ncp_use)),
      error = function(e) NULL
    )
    if (is.null(imp) || is.null(imp$completeObs) || any(!is.finite(imp$completeObs))) {
      ncp_use <- 1
      imp <- missMDA::imputePCA(X, ncp = ncp_use)
    }
    X_imp <- imp$completeObs
    variable_features <- apply(X_imp, 2, function(values) {
      values <- values[is.finite(values)]
      length(values) >= 2 && stats::sd(values) > 0
    })
    X_imp <- X_imp[, variable_features, drop = FALSE]
    validate(
      need(ncol(X_imp) >= 2, "Fewer than 2 variable proteins remain after PCA imputation. Lower the PCA threshold or use a less filtered abundance table.")
    )
    X_proc <- if (isTRUE(input$zscore)) scale(X_imp) else X_imp
    X_proc <- as.matrix(X_proc)
    validate(
      need(all(is.finite(X_proc)), "PCA matrix contains non-finite values after preprocessing. Try turning off z-scoring or lowering the minimum observed-samples threshold.")
    )
    pca <- FactoMineR::PCA(X_proc, graph = FALSE)

    scores <- as.data.frame(pca$ind$coord[, 1:2, drop = FALSE])
    colnames(scores) <- c("PC1", "PC2")
    scores$Sample <- rownames(scores)
    scores$RowOrder <- seq_len(nrow(scores))

    parsed_label <- stringr::str_extract(scores$Sample, "ID[0-9]+_[^_]+")
    scores$ParsedLabel <- ifelse(is.na(parsed_label), scores$Sample, parsed_label)
    scores$Label <- scores$Sample

    ann <- detected_annotations()
    if (!is.null(ann)) {
      scores <- scores %>% left_join(ann, by = "Sample")
      scores$DisplayLabel <- make_display_label(scores)
    }

    if (input$color_by == "Auto row groups") {
      scores$AutoColorGroup <- paste0(input$color_group_prefix, "_", ceiling(scores$RowOrder / input$rows_per_color_group))
    } else {
      scores$AutoColorGroup <- NA_character_
    }

    if (input$shape_by == "Auto row groups") {
      scores$AutoShapeGroup <- paste0(input$shape_group_prefix, "_", ceiling(scores$RowOrder / input$rows_per_shape_group))
    } else {
      scores$AutoShapeGroup <- NA_character_
    }

    if (!is.null(project_file("meta_file")) || !is.null(project_file("run_order_file")) || !is.null(project_file("sample_details_file"))) {
      scores <- scores %>% left_join(active_metadata(), by = "Sample")
    }
    if (!is.null(input$label_mode) && input$label_mode %in% colnames(scores)) {
      labels <- as.character(scores[[input$label_mode]])
      missing_labels <- is.na(labels) | labels == ""
      if (isTRUE(input$label_missing_points)) {
        labels[missing_labels] <- scores$Sample[missing_labels]
      } else {
        labels[missing_labels] <- NA_character_
      }
      scores$Label <- labels
    }

    list(
      scores = scores,
      pca = pca,
      retained_features = retained_features,
      variable_features = ncol(X_imp),
      total_features = total_features,
      ncp_used = ncp_use,
      missing_values_before_imputation = sum(is.na(X))
    )
  })

  pca_plot_obj <- reactive({
    res <- pca_results()
    df <- res$scores
    color_var <- NULL
    shape_var <- NULL

    if (input$color_by == "Auto row groups") {
      color_var <- "AutoColorGroup"
    } else if (startsWith(input$color_by, "Auto-detected: ")) {
      candidate <- sub("^Auto-detected: ", "", input$color_by)
      if (candidate %in% colnames(df)) color_var <- candidate
    } else if (input$color_by != "None" && input$color_by %in% colnames(df)) {
      color_var <- input$color_by
    }

    if (input$shape_by == "Auto row groups") {
      shape_var <- "AutoShapeGroup"
    } else if (startsWith(input$shape_by, "Auto-detected: ")) {
      candidate <- sub("^Auto-detected: ", "", input$shape_by)
      if (candidate %in% colnames(df)) shape_var <- candidate
    } else if (input$shape_by != "None" && input$shape_by %in% colnames(df)) {
      shape_var <- input$shape_by
    }

    p <- ggplot(df, aes(x = PC1, y = PC2))
    if (!is.null(color_var) && !is.null(shape_var)) {
      p <- p + geom_point(aes_string(color = color_var, shape = shape_var), size = input$point_size, alpha = 0.85)
    } else if (!is.null(color_var)) {
      p <- p + geom_point(aes_string(color = color_var), size = input$point_size, alpha = 0.85)
    } else if (!is.null(shape_var)) {
      p <- p + geom_point(aes_string(shape = shape_var), size = input$point_size, alpha = 0.85)
    } else {
      p <- p + geom_point(size = input$point_size, alpha = 0.85)
    }

    if (isTRUE(input$label_points)) {
      label_df <- df[!is.na(df$Label) & df$Label != "", , drop = FALSE]
      p <- p + geom_text(data = label_df, aes(label = Label), vjust = -0.7, size = input$label_size)
    }

    p + theme_bw(base_size = 14) + coord_fixed() +
      labs(
        title = input$plot_title,
        subtitle = input$plot_subtitle,
        x = paste0("PC1 (", round(res$pca$eig[1, 2], 2), "%)"),
        y = paste0("PC2 (", round(res$pca$eig[2, 2], 2), "%)")
      ) +
      theme(
        axis.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = input$feature_title_size, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(face = "bold")
      )
  })

  feature_plot_df <- reactive({
    features <- unique(as.character(input$feature_select))
    features <- features[!is.na(features) & nzchar(features)]
    req(length(features) > 0)
    expr <- feature_expression_data()
    features <- intersect(features, as.character(expr$Feature))
    validate(need(length(features) > 0, "Selected feature(s) not found."))

    source_is_log2 <- identical(input$feature_data_source, "S3_batch_corrected")
    value_scale <- input$feature_value_scale
    if (is.null(value_scale) || !value_scale %in% c("raw", "log2")) value_scale <- "raw"

    base <- data.frame(Sample = colnames(expr)[-1], stringsAsFactors = FALSE)
    base <- base %>% left_join(sample_metadata_for_plotting(base$Sample), by = "Sample")
    label_col <- input$feature_label_by
    base$DisplayLabel <- if (!is.null(label_col) && label_col %in% colnames(base)) as.character(base[[label_col]]) else as.character(base$Sample)
    missing_labels <- is.na(base$DisplayLabel) | !nzchar(base$DisplayLabel)
    base$DisplayLabel[missing_labels] <- as.character(base$Sample[missing_labels])
    group_col <- input$feature_group_by
    base$GroupValue <- if (!is.null(group_col) && group_col %in% colnames(base)) as.character(base[[group_col]]) else "All samples"
    base$GroupValue[is.na(base$GroupValue) | !nzchar(base$GroupValue)] <- "Missing"

    base$.input_order <- seq_len(nrow(base))
    order_columns <- input$feature_order_columns
    order_columns <- order_columns[order_columns %in% colnames(base)]
    if (length(order_columns) > 0) {
      order_args <- lapply(order_columns, function(column) {
        values <- base[[column]]
        if (inherits(values, c("numeric", "integer", "Date", "POSIXct", "POSIXlt"))) return(values)
        numeric_values <- suppressWarnings(as.numeric(as.character(values)))
        if (sum(is.finite(numeric_values)) >= max(2, floor(0.8 * length(values)))) {
          numeric_values[!is.finite(numeric_values)] <- Inf
          return(numeric_values)
        }
        values <- as.character(values)
        values[is.na(values) | !nzchar(values)] <- "zzzz_missing"
        values
      })
      order_index <- do.call(order, c(order_args, list(base$.input_order, na.last = TRUE)))
      base <- base[order_index, , drop = FALSE]
    }
    base$.plot_order <- seq_len(nrow(base))
    base$.BarID <- factor(base$Sample, levels = unique(base$Sample))
    display_levels <- unique(as.character(base$DisplayLabel))

    rows <- lapply(features, function(feature) {
      row <- expr[expr$Feature == feature, , drop = FALSE]
      vals <- suppressWarnings(as.numeric(row[1, -1, drop = TRUE]))
      names(vals) <- colnames(expr)[-1]
      out <- base
      out$Feature <- feature
      out$Value <- unname(vals[out$Sample])
      out <- out[!is.na(out$Value), , drop = FALSE]
      if (!nrow(out)) return(out)
      out$RawValue <- if (source_is_log2) {
        2^out$Value
      } else {
        out$Value
      }
      out$Log2Value <- if (source_is_log2) {
        out$Value
      } else {
        values <- suppressWarnings(log2(out$Value))
        values[!is.finite(values)] <- NA_real_
        values
      }
      out$PlotValue <- if (identical(value_scale, "log2")) out$Log2Value else out$RawValue
      out <- out[is.finite(out$PlotValue), , drop = FALSE]
      if (!nrow(out)) return(out)
      out$PlotCenteredValue <- out$PlotValue - mean(out$PlotValue, na.rm = TRUE)
      val_sd <- stats::sd(out$PlotValue, na.rm = TRUE)
      out$PlotZValue <- if (is.na(val_sd) || val_sd == 0) 0 else as.numeric(scale(out$PlotValue))
      out
    })
    out <- do.call(rbind, rows)
    validate(need(!is.null(out) && nrow(out) > 0, "No finite values available for the selected feature(s) at this plot scale."))
    out$Feature <- factor(as.character(out$Feature), levels = features)
    out$DisplayLabel <- factor(as.character(out$DisplayLabel), levels = display_levels)
    out
  })

  observe({
    md <- tryCatch(active_metadata(), error = function(e) NULL)
    if (is.null(md)) {
      updateSelectInput(session, "box_group_by", choices = c("Sample"), selected = "Sample")
      updateSelectInput(session, "script_box_group_by", choices = c("Sample"), selected = "Sample")
      updateSelectInput(session, "script_box_label_by", choices = c("None", "Sample"), selected = "None")
      return()
    }

    choices <- setdiff(scoped_metadata_columns(md), c("SampleDetailsID", "File Name"))
    current <- isolate(input$box_group_by)
    selected <- if (!is.null(current) && current %in% choices) {
      current
    } else if ("Condition" %in% choices) {
      "Condition"
    } else {
      choices[1]
    }
    selected <- retain_metadata_choice(current, choices, "Condition")
    field_choices <- c("Select grouping field..." = "", stats::setNames(choices, choices))
    updateSelectInput(session, "box_group_by", choices = field_choices, selected = selected)
    current_script <- isolate(input$script_box_group_by)
    selected_script <- if (!is.null(current_script) && current_script %in% choices) {
      current_script
    } else if ("Condition" %in% choices) {
      "Condition"
    } else {
      choices[1]
    }
    selected_script <- retain_metadata_choice(current_script, choices, "Condition")
    updateSelectInput(session, "script_box_group_by", choices = field_choices, selected = selected_script)
    label_choices <- c("None" = "None", stats::setNames(choices, choices))
    current_label <- isolate(input$script_box_label_by)
    selected_label <- if (!is.null(current_label) && current_label %in% unname(label_choices)) current_label else "None"
    updateSelectInput(session, "script_box_label_by", choices = label_choices, selected = selected_label)
  })

  script_box_expression_data <- reactive({
    expression_data_for_protein_source(input$script_box_source)
  })

  observeEvent(script_box_expression_data(), {
    expr <- script_box_expression_data()
    choices <- protein_feature_choices(input$script_box_source, expr)
    current <- isolate(input$script_box_features)
    choice_values <- unname(choices)
    current <- current[current %in% choice_values]
    selected <- if (length(current)) current else choice_values[1]
    freezeReactiveValue(input, "script_box_features")
    updateSelectizeInput(session, "script_box_features", choices = choices, selected = selected, server = TRUE)
  }, ignoreInit = FALSE)

  observe({
    pending <- pending_box_selection()
    req(!is.null(pending), identical(input$script_box_source, pending$source))
    expr <- script_box_expression_data()
    choices <- protein_feature_choices(input$script_box_source, expr)
    selected <- intersect(pending$features, unname(choices))
    pending_box_selection(NULL)
    req(length(selected))
    session$onFlushed(function() {
      freezeReactiveValue(input, "script_box_features")
      updateSelectizeInput(session, "script_box_features", choices = choices, selected = selected, server = TRUE)
    }, once = TRUE)
  })

  script_box_all_groups <- reactive({
    expr <- script_box_expression_data()
    df <- data.frame(Sample = colnames(expr)[-1], stringsAsFactors = FALSE)
    df <- df %>% left_join(sample_metadata_for_plotting(df$Sample), by = "Sample")
    group_col <- input$script_box_group_by
    validate(need(!is.null(group_col) && group_col %in% colnames(df), "Choose an available metadata column for grouping."))
    groups <- as.character(df[[group_col]])
    groups[is.na(groups) | !nzchar(groups)] <- "Missing"
    unique(groups)
  })

  observeEvent(script_box_all_groups(), {
    groups <- script_box_all_groups()
    current <- isolate(input$script_box_conditions)
    selected <- current[current %in% groups]
    if (length(selected) == 0) selected <- groups
    updateSelectizeInput(session, "script_box_conditions", choices = groups, selected = selected, server = TRUE)
  }, ignoreInit = FALSE)

  script_box_plot_df <- reactive({
    req(input$script_box_features)
    expr <- script_box_expression_data()
    selected_features <- intersect(as.character(input$script_box_features), as.character(expr$Feature))
    validate(need(length(selected_features) > 0, "Selected proteins were not found."))
    df <- build_multifeature_boxplot_data(expr, selected_features)
    df <- df %>% left_join(sample_metadata_for_plotting(df$Sample), by = "Sample")
    df <- deduplicate_feature_samples(df)
    group_col <- input$script_box_group_by
    validate(need(!is.null(group_col) && group_col %in% colnames(df), "Choose an available metadata column for grouping."))
    df$GroupValue <- as.character(df[[group_col]])
    df$GroupValue[is.na(df$GroupValue) | !nzchar(df$GroupValue)] <- "Missing"
    selected_groups <- input$script_box_conditions
    if (!is.null(selected_groups) && length(selected_groups) > 0) {
      df <- df[df$GroupValue %in% selected_groups, , drop = FALSE]
      df$GroupValue <- factor(df$GroupValue, levels = selected_groups)
    } else {
      df$GroupValue <- factor(df$GroupValue, levels = unique(df$GroupValue))
    }
    label_col <- input$script_box_label_by
    df$PointLabel <- if (!is.null(label_col) && !identical(label_col, "None") && label_col %in% colnames(df)) as.character(df[[label_col]]) else ""
    df$PointLabel[is.na(df$PointLabel)] <- ""
    source_is_log2 <- identical(input$script_box_source, "S3_batch_corrected")
    df$PlotValue <- if (identical(input$script_box_value_scale, "log2")) {
      if (source_is_log2) {
        suppressWarnings(as.numeric(df$Value))
      } else {
        out <- suppressWarnings(log2(as.numeric(df$Value)))
        out[!is.finite(out)] <- NA_real_
        out
      }
    } else {
      if (source_is_log2) {
        2^suppressWarnings(as.numeric(df$Value))
      } else {
        suppressWarnings(as.numeric(df$Value))
      }
    }
    df <- df[is.finite(df$PlotValue), , drop = FALSE]
    validate(need(nrow(df) > 0, "No finite abundance values are available for this script-style boxplot."))
    df
  })

  box_plot_df <- reactive({
    df <- feature_plot_df()
    group_col <- input$box_group_by
    validate(need(!is.null(group_col) && group_col %in% colnames(df), "Choose an available metadata column for grouping."))

    df$GroupValue <- as.character(df[[group_col]])
    df$GroupValue[is.na(df$GroupValue) | !nzchar(df$GroupValue)] <- "Missing"
    df$GroupValue <- factor(df$GroupValue, levels = unique(df$GroupValue))
    source_is_log2 <- identical(input$feature_data_source, "S3_batch_corrected")
    df$PlotValue <- if (identical(input$box_value_scale, "log2")) {
      if (source_is_log2) {
        suppressWarnings(as.numeric(df$Value))
      } else {
        values <- suppressWarnings(log2(as.numeric(df$Value)))
        values[!is.finite(values)] <- NA_real_
        values
      }
    } else {
      if (source_is_log2) {
        2^suppressWarnings(as.numeric(df$Value))
      } else {
        suppressWarnings(as.numeric(df$Value))
      }
    }

    df <- df[is.finite(df$PlotValue), , drop = FALSE]
    validate(need(nrow(df) > 0, "No finite abundance values are available for this protein group plot."))
    df
  })

  box_plot_obj <- reactive({
    df <- box_plot_df()
    title_text <- if (nzchar(input$box_plot_title)) input$box_plot_title else input$feature_select
    y_title <- if (nzchar(input$box_y_axis_title)) {
      input$box_y_axis_title
    } else if (identical(input$box_value_scale, "log2")) {
      "Log2 protein abundance"
    } else {
      "Protein abundance"
    }

    group_levels <- levels(df$GroupValue)
    condition_colors <- c(CKD = "#00A651", CKDu = "#597DFF", Control = "#ED3438", SPQC = "#FF00CB")
    other_groups <- setdiff(group_levels, names(condition_colors))
    if (length(other_groups)) {
      condition_colors <- c(condition_colors, stats::setNames(scales::hue_pal()(length(other_groups)), other_groups))
    }

    p <- ggplot(df, aes(x = GroupValue, y = PlotValue, color = GroupValue))
    if (identical(input$box_plot_style, "boxplot")) {
      p <- p + geom_boxplot(width = 0.5, fill = NA, linewidth = 0.7, outlier.shape = NA)
    } else {
      mean_sd <- function(values) {
        average <- mean(values, na.rm = TRUE)
        spread <- stats::sd(values, na.rm = TRUE)
        if (!is.finite(spread)) spread <- 0
        data.frame(y = average, ymin = average - spread, ymax = average + spread)
      }
      p <- p +
        stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.17, linewidth = 0.8) +
        stat_summary(fun = mean, geom = "crossbar", width = 0.30, linewidth = 0.8)
    }

    p +
      geom_jitter(width = 0.1, height = 0, shape = 21, fill = "white", stroke = 0.9, size = input$box_point_size) +
      scale_color_manual(values = condition_colors, drop = FALSE) +
      labs(title = title_text, x = NULL, y = y_title, color = NULL) +
      theme_classic(base_size = input$box_text_size) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text = element_text(face = "bold", color = "black"),
        axis.text.x = element_text(angle = 50, hjust = 1),
        axis.title.y = element_text(face = "bold"),
        legend.position = "none"
      )
  })

  script_box_plot_obj <- reactive({
    df <- script_box_plot_df()
    script_box_ncol_value <- suppressWarnings(as.integer(input$script_box_ncol))
    if (is.na(script_box_ncol_value) || script_box_ncol_value < 1) script_box_ncol_value <- 1L
    text_scale <- facet_text_scale(script_box_ncol_value)
    scaled_text_size <- max(8, input$script_box_text_size * text_scale)
    scaled_label_size <- max(1.8, input$script_box_label_size * text_scale)
    selected_features <- unique(as.character(df$Feature))
    title_text <- if (nzchar(input$script_box_title)) input$script_box_title else if (length(selected_features) == 1L) selected_features else "Selected protein groups"
    y_title <- if (nzchar(input$script_box_y_axis_title)) {
      input$script_box_y_axis_title
    } else if (identical(input$script_box_value_scale, "log2")) {
      "Log2 protein abundance"
    } else {
      "Protein abundance"
    }
    group_levels <- levels(df$GroupValue)
    condition_colors <- c(
      CKD = "#00A651",
      CKDu = "#597DFF",
      Control = "#ED3438",
      SPQC = "#FF00CB",
      CTL = "blue",
      IFN = "orange",
      IFN_Pred = "grey60",
      Pred = "darkred"
    )
    other_groups <- setdiff(group_levels, names(condition_colors))
    if (length(other_groups)) {
      condition_colors <- c(condition_colors, stats::setNames(scales::hue_pal()(length(other_groups)), other_groups))
    }

    p <- ggplot(df, aes(x = GroupValue, y = PlotValue, fill = GroupValue))
    if (identical(input$script_box_plot_style, "mean_sd")) {
      mean_sd <- function(values) {
        average <- mean(values, na.rm = TRUE)
        spread <- stats::sd(values, na.rm = TRUE)
        if (!is.finite(spread)) spread <- 0
        data.frame(y = average, ymin = average - spread, ymax = average + spread)
      }
      p <- p +
        stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.18, linewidth = 0.9, color = "black") +
        stat_summary(fun = mean, geom = "crossbar", width = 0.35, linewidth = 0.9, color = "black", fill = "white")
    } else {
      p <- p +
        geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black", width = 0.55)
    }

    point_position <- position_jitter(width = 0.12, height = 0, seed = 1)
    p <- p +
      geom_point(
        aes(fill = GroupValue),
        shape = 21,
        color = "black",
        size = input$script_box_point_size,
        position = point_position,
        stroke = 0.8,
        alpha = input$script_box_point_opacity
      ) +
      scale_fill_manual(values = condition_colors, drop = FALSE) +
      facet_wrap(~Feature, scales = "free_y", ncol = script_box_ncol_value) +
      labs(title = title_text, y = y_title, x = input$script_box_group_by, fill = NULL) +
      theme_classic(base_size = scaled_text_size) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, face = "bold", color = "black"),
        axis.text.y = element_text(face = "bold", color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none"
      )
    if (any(nzchar(df$PointLabel))) {
      p <- p + geom_text(aes(label = PointLabel), position = point_position, vjust = -0.8, size = scaled_label_size, check_overlap = TRUE, show.legend = FALSE)
    }
    if (requireNamespace("ggprism", quietly = TRUE)) {
      p <- p + ggprism::theme_prism(base_size = scaled_text_size) +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)
        )
    }
    p
  })

  correlation_expression_data <- reactive({
    expr <- expression_data_for_protein_source(input$correlation_source)
    source_is_log2 <- identical(input$correlation_source, "S3_batch_corrected")
    out <- expr
    numeric_cols <- setdiff(colnames(out), "Feature")
    out[numeric_cols] <- lapply(out[numeric_cols], function(values) {
      values <- suppressWarnings(as.numeric(as.character(values)))
      if (identical(input$correlation_value_scale, "log2")) {
        if (source_is_log2) {
          values
        } else {
          values <- log2(values)
          values[!is.finite(values)] <- NA_real_
          values
        }
      } else {
        if (source_is_log2) {
          2^values
        } else {
          values
        }
      }
    })
    out
  })

  observeEvent(correlation_expression_data(), {
    expr <- correlation_expression_data()
    choices <- protein_feature_choices(input$correlation_source, expr)
    current <- isolate(input$correlation_feature)
    choice_values <- unname(choices)
    selected <- if (!is.null(current) && current %in% choice_values) current else choice_values[1]
    updateSelectizeInput(session, "correlation_feature", choices = choices, selected = selected, server = TRUE)
  }, ignoreInit = FALSE)

  correlation_sample_metadata <- reactive({
    expr <- correlation_expression_data()
    samples <- setdiff(colnames(expr), "Feature")
    sample_metadata_for_plotting(samples)
  })

  observeEvent(list(correlation_sample_metadata(), downstream_metadata_columns()), {
    md <- correlation_sample_metadata()
    choices <- setdiff(scoped_metadata_columns(md), c("SampleDetailsID", "File Name"))
    choices <- choices[vapply(md[choices], function(values) length(unique(na.omit(as.character(values)))) > 0, logical(1))]
    group_choices <- c("Select grouping field..." = "", "All samples" = "None", stats::setNames(choices, choices))
    current_group_by <- isolate(input$correlation_group_by)
    selected_group_by <- if (!is.null(current_group_by) && current_group_by %in% group_choices) {
      current_group_by
    } else if ("Group" %in% choices) {
      "Group"
    } else if ("Condition" %in% choices) {
      "Condition"
    } else {
      "None"
    }
    selected_group_by <- retain_metadata_choice(current_group_by, unname(group_choices), c("Group", "Condition", "None"))
    updateSelectInput(session, "correlation_group_by", choices = group_choices, selected = selected_group_by)

    covariate_choices <- choices[!choices %in% c("Sample", "AnalysisLabel")]
    current_covariates <- isolate(input$correlation_covariates)
    selected_covariates <- current_covariates[current_covariates %in% covariate_choices]
    updateSelectizeInput(session, "correlation_covariates", choices = covariate_choices, selected = selected_covariates, server = TRUE)
  }, ignoreInit = FALSE)

  observeEvent(list(correlation_sample_metadata(), input$correlation_group_by), {
    md <- correlation_sample_metadata()
    group_col <- input$correlation_group_by
    if (is.null(group_col) || identical(group_col, "None") || !group_col %in% colnames(md)) {
      updateSelectizeInput(session, "correlation_groups", choices = character(0), selected = character(0), server = TRUE)
      return()
    }
    groups <- unique(as.character(md[[group_col]]))
    groups <- groups[!is.na(groups) & nzchar(groups)]
    current <- isolate(input$correlation_groups)
    selected <- current[current %in% groups]
    updateSelectizeInput(session, "correlation_groups", choices = groups, selected = selected, server = TRUE)
  }, ignoreInit = FALSE)

  correlation_filtered_inputs <- reactive({
    expr <- correlation_expression_data()
    md <- correlation_sample_metadata()
    samples <- setdiff(colnames(expr), "Feature")
    md <- md[match(samples, md$Sample), , drop = FALSE]

    keep <- rep(TRUE, length(samples))
    group_col <- input$correlation_group_by
    selected_groups <- input$correlation_groups
    if (!is.null(group_col) && !identical(group_col, "None") && group_col %in% colnames(md) &&
        !is.null(selected_groups) && length(selected_groups) > 0) {
      keep <- as.character(md[[group_col]]) %in% selected_groups
    }

    samples <- samples[keep]
    validate(need(length(samples) >= 3, "Need at least 3 samples after within-group filtering."))
    list(
      expr = expr[, c("Feature", samples), drop = FALSE],
      metadata = md[keep, , drop = FALSE],
      samples = samples,
      group_filter = if (!is.null(group_col) && !identical(group_col, "None") && length(selected_groups) > 0) {
        paste0(group_col, " in ", paste(selected_groups, collapse = ", "))
      } else {
        "all samples"
      }
    )
  })

  correlation_results_data <- reactive({
    req(input$correlation_feature)
    inputs <- correlation_filtered_inputs()
    expr <- inputs$expr
    validate(
      need(nrow(expr) >= 2, "Need at least 2 protein features for correlation analysis."),
      need(input$correlation_feature %in% expr$Feature, "Selected reference protein was not found.")
    )

    samples <- setdiff(colnames(expr), "Feature")
    validate(need(length(samples) >= 3, "Need at least 3 matched sample abundance columns for regression."))

    mat <- as.matrix(expr[, samples, drop = FALSE])
    storage.mode(mat) <- "numeric"
    rownames(mat) <- expr$Feature
    ref <- as.numeric(mat[input$correlation_feature, , drop = TRUE])
    validate(need(sum(is.finite(ref)) >= 3, "Reference protein needs at least 3 finite sample values."))

    covariates <- input$correlation_covariates
    covariates <- covariates[covariates %in% colnames(inputs$metadata)]
    method <- input$correlation_method

    build_model_frame <- function(y, ok) {
      model_data <- data.frame(
        Target = y[ok],
        Reference = ref[ok],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      used_covariates <- character(0)
      if (identical(method, "linear") && length(covariates) > 0) {
        for (covariate in covariates) {
          values <- inputs$metadata[[covariate]][ok]
          if (all(is.na(values))) next
          unique_values <- unique(na.omit(as.character(values)))
          if (length(unique_values) < 2) next
          safe_name <- make.names(covariate, unique = TRUE)
          if (is.numeric(values) || is.integer(values)) {
            model_data[[safe_name]] <- as.numeric(values)
          } else {
            model_data[[safe_name]] <- factor(as.character(values))
          }
          used_covariates <- c(used_covariates, safe_name)
        }
      }
      model_data <- model_data[stats::complete.cases(model_data), , drop = FALSE]
      list(data = model_data, covariates = used_covariates)
    }

    fit_one <- function(y) {
      ok <- is.finite(ref) & is.finite(y)
      n <- sum(ok)
      if (n < 3 || stats::sd(ref[ok]) == 0 || stats::sd(y[ok]) == 0) {
        return(c(Beta = NA_real_, R2 = NA_real_, PValue = NA_real_, Correlation = NA_real_, N = n))
      }

      if (identical(method, "pearson") || identical(method, "spearman")) {
        test <- tryCatch(stats::cor.test(ref[ok], y[ok], method = method, exact = FALSE), error = function(e) NULL)
        if (is.null(test)) {
          return(c(Beta = NA_real_, R2 = NA_real_, PValue = NA_real_, Correlation = NA_real_, N = n))
        }
        correlation <- unname(test$estimate)
        beta <- if (stats::var(ref[ok]) == 0) NA_real_ else stats::cov(ref[ok], y[ok]) / stats::var(ref[ok])
        return(c(Beta = beta, R2 = correlation^2, PValue = unname(test$p.value), Correlation = correlation, N = n))
      }

      model <- build_model_frame(y, ok)
      if (nrow(model$data) < 3 || stats::sd(model$data$Reference) == 0 || stats::sd(model$data$Target) == 0) {
        return(c(Beta = NA_real_, R2 = NA_real_, PValue = NA_real_, Correlation = NA_real_, N = nrow(model$data)))
      }
      formula_text <- if (length(model$covariates) > 0) {
        paste("Target ~ Reference +", paste(model$covariates, collapse = " + "))
      } else {
        "Target ~ Reference"
      }
      fit <- tryCatch(stats::lm(stats::as.formula(formula_text), data = model$data), error = function(e) NULL)
      if (is.null(fit)) {
        return(c(Beta = NA_real_, R2 = NA_real_, PValue = NA_real_, Correlation = NA_real_, N = nrow(model$data)))
      }
      fit_summary <- summary(fit)
      coefficient_table <- fit_summary$coefficients
      if (!"Reference" %in% rownames(coefficient_table)) {
        return(c(Beta = NA_real_, R2 = NA_real_, PValue = NA_real_, Correlation = NA_real_, N = nrow(model$data)))
      }
      beta <- unname(coefficient_table["Reference", "Estimate"])
      p_value <- unname(coefficient_table["Reference", "Pr(>|t|)"])
      r_squared <- unname(fit_summary$r.squared)
      correlation <- suppressWarnings(stats::cor(model$data$Reference, model$data$Target))
      c(Beta = beta, R2 = r_squared, PValue = p_value, Correlation = correlation, N = nrow(model$data))
    }

    stats_mat <- t(apply(mat, 1, fit_one))
    out <- data.frame(
      Protein = rownames(stats_mat),
      Beta = as.numeric(stats_mat[, "Beta"]),
      R2 = as.numeric(stats_mat[, "R2"]),
      PValue = as.numeric(stats_mat[, "PValue"]),
      QValue = stats::p.adjust(as.numeric(stats_mat[, "PValue"]), method = "BH"),
      Correlation = as.numeric(stats_mat[, "Correlation"]),
      N = as.integer(stats_mat[, "N"]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    out$AbsBeta <- abs(out$Beta)
    out$Method <- switch(method, pearson = "Pearson correlation", spearman = "Spearman correlation", "Linear regression")
    out$GroupFilter <- inputs$group_filter
    out$Covariates <- if (identical(method, "linear") && length(covariates) > 0) paste(covariates, collapse = "; ") else ""

    if (isTRUE(input$correlation_exclude_reference)) {
      out <- out[out$Protein != input$correlation_feature, , drop = FALSE]
    }

    info <- protein_info_for_source(input$correlation_source, input$report_feature_col)
    info_cols <- intersect(c("Protein", "PG.Genes", "PG.ProteinGroups", "PG.ProteinNames", "PG.ProteinDescriptions", "RawFeatureID"), colnames(info))
    info <- info[, info_cols, drop = FALSE]
    info <- info[!duplicated(info$Protein), , drop = FALSE]
    out <- out %>% left_join(info, by = "Protein")
    out$DisplayProtein <- best_protein_display_label(out, fallback_col = "Protein")
    out <- out[order(out$QValue, -out$R2, -out$AbsBeta, na.last = TRUE), , drop = FALSE]
    rownames(out) <- NULL
    out
  })

  correlation_lollipop_data <- reactive({
    data <- correlation_results_data()
    data <- data[is.finite(data$Beta) & is.finite(data$R2), , drop = FALSE]
    validate(need(nrow(data) > 0, "No valid regression results are available for plotting."))

    rank_by <- input$correlation_rank_by
    if (identical(rank_by, "r_squared")) {
      data <- data[order(-data$R2, data$QValue, na.last = TRUE), , drop = FALSE]
    } else if (identical(rank_by, "abs_beta")) {
      data <- data[order(-data$AbsBeta, data$QValue, na.last = TRUE), , drop = FALSE]
    } else if (identical(rank_by, "p_value")) {
      data <- data[order(data$PValue, -data$R2, na.last = TRUE), , drop = FALSE]
    } else {
      data <- data[order(data$QValue, -data$R2, na.last = TRUE), , drop = FALSE]
    }
    top_n <- min(max(1, input$correlation_top_n), nrow(data))
    data <- data[seq_len(top_n), , drop = FALSE]
    data$PlotLabel <- make.unique(as.character(data$DisplayProtein))
    data$NegLog10Q <- -log10(pmax(data$QValue, .Machine$double.xmin))
    data$PlotLabel <- factor(data$PlotLabel, levels = rev(data$PlotLabel))
    data
  })

  correlation_lollipop_plot_obj <- reactive({
    data <- correlation_lollipop_data()
    title_text <- if (nzchar(input$correlation_plot_title)) input$correlation_plot_title else "Top Protein Correlations"
    subtitle_text <- paste0("Reference: ", input$correlation_feature)
    x_title <- if (identical(input$correlation_value_scale, "log2")) {
      "Beta coefficient (log2 abundance)"
    } else {
      "Beta coefficient (raw abundance)"
    }

    ggplot(data, aes(y = PlotLabel, x = Beta)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
      geom_segment(aes(x = 0, xend = Beta, yend = PlotLabel, color = R2), linewidth = 0.9) +
      geom_point(aes(color = R2, size = NegLog10Q), alpha = 0.9) +
      scale_color_gradient(low = "#4575B4", high = "#D73027", name = "R2") +
      scale_size_continuous(name = "-log10(q)", range = c(2.5, 7)) +
      labs(title = title_text, subtitle = subtitle_text, x = x_title, y = NULL) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(color = "black"),
        legend.position = "right"
      )
  })

  feature_plot_obj <- reactive({
    df <- feature_plot_df()
    validate(need(nrow(df) > 0, "No numeric abundance values available for the selected feature(s)."))

    has_group <- "GroupValue" %in% colnames(df)
    selected_features <- levels(df$Feature)
    if (is.null(selected_features) || !length(selected_features)) selected_features <- unique(as.character(df$Feature))
    title_text <- if (nzchar(input$feature_plot_title)) {
      input$feature_plot_title
    } else if (length(selected_features) == 1L) {
      selected_features
    } else {
      "Selected proteins"
    }

    fill_var <- switch(
      input$feature_color_mode,
      "Z-score" = "PlotZValue",
      "Centered value" = "PlotCenteredValue",
      "Raw value" = "PlotValue",
      "PlotZValue"
    )
    y_axis_title <- if (identical(input$feature_value_scale, "log2")) "Log2 protein abundance" else "Protein abundance"
    group_text <- if (has_group) as.character(df$GroupValue) else "Not available"
    batch_text <- if ("Batch" %in% colnames(df)) as.character(df$Batch) else "Not available"
    condition_text <- if ("Condition" %in% colnames(df)) as.character(df$Condition) else "Not available"
    df$.FeatureHover <- paste0(
      "Protein: ", as.character(df$Feature),
      "<br>Sample: ", as.character(df$Sample),
      "<br>Label: ", as.character(df$DisplayLabel),
      "<br>Value: ", signif(df$PlotValue, 5),
      "<br>Group: ", group_text,
      "<br>Condition: ", condition_text,
      "<br>Batch: ", batch_text
    )

    p <- ggplot(df, aes(x = .BarID, y = PlotValue, group = .plot_order, text = .FeatureHover))
    bar_labels <- stats::setNames(as.character(df$DisplayLabel), as.character(df$.BarID))
    p <- p + scale_x_discrete(labels = bar_labels)

    if (has_group && input$feature_group_style == "Outline color") {
      p <- p + geom_col(aes_string(fill = fill_var, color = "GroupValue"),
                        width = input$feature_bar_width, linewidth = 0.7, na.rm = TRUE)
    } else {
      p <- p + geom_col(aes_string(fill = fill_var),
                        width = input$feature_bar_width, color = "black", linewidth = 0.4, na.rm = TRUE)
    }

    if (isTRUE(input$show_feature_mean) && any(!is.na(df$PlotValue))) {
      mean_df <- df %>%
        dplyr::group_by(Feature) %>%
        dplyr::summarise(MeanValue = mean(PlotValue, na.rm = TRUE), .groups = "drop")
      p <- p +
        geom_hline(data = mean_df, aes(yintercept = MeanValue), linetype = "dashed", inherit.aes = FALSE) +
        geom_text(data = mean_df, aes(x = -Inf, y = MeanValue, label = "Mean"), hjust = -0.1, vjust = -0.4, size = 3.5, inherit.aes = FALSE)
    }

    if (input$feature_color_mode %in% c("Z-score", "Centered value")) {
      if (isTRUE(input$feature_symmetric_scale)) {
        zmax <- max(abs(df[[fill_var]]), na.rm = TRUE)
        if (!is.finite(zmax) || zmax == 0) zmax <- 1
        p <- p + scale_fill_gradient2(
          low = "#2166AC",
          mid = "white",
          high = "#B2182B",
          midpoint = 0,
          limits = c(-zmax, zmax),
          name = input$feature_color_mode
        )
      } else {
        p <- p + scale_fill_gradient2(
          low = "#2166AC",
          mid = "white",
          high = "#B2182B",
          midpoint = 0,
          name = input$feature_color_mode
        )
      }
    } else {
      p <- p + scale_fill_gradient(low = "white", high = "#B2182B", name = input$feature_color_mode)
    }

    if (has_group) {
      group_vals <- unique(as.character(df$GroupValue))
      if (all(c("non_smoker", "cigarette", "ecig", "SPQC") %in% group_vals)) {
        p <- p + scale_color_manual(values = c(
          "non_smoker" = "firebrick3",
          "cigarette" = "black",
          "ecig" = "dodgerblue3",
          "SPQC" = "gray50"
        ), drop = FALSE)
      }
    }

    ncol_value <- suppressWarnings(as.integer(input$feature_plot_ncol))
    if (is.na(ncol_value) || ncol_value < 1) ncol_value <- 1L
    text_scale <- facet_text_scale(ncol_value)
    scaled_feature_text_size <- max(1.8, input$feature_text_size * text_scale)
    scaled_base_size <- max(8, 13 * text_scale)

    if (has_group && input$feature_group_style == "Colored x labels") {
      label_df <- df %>%
        dplyr::group_by(Feature) %>%
        dplyr::mutate(LabelY = min(PlotValue, na.rm = TRUE) - 0.08 * max(diff(range(PlotValue, na.rm = TRUE)), 1)) %>%
        dplyr::ungroup()

      p <- p +
        geom_text(
          data = label_df,
          aes(x = .BarID, y = LabelY, label = DisplayLabel, color = GroupValue),
          inherit.aes = FALSE,
          angle = if (isTRUE(input$rotate_feature_labels)) 90 else 0,
          hjust = 1,
          vjust = 0.5,
          size = scaled_feature_text_size,
          show.legend = TRUE
        ) +
        scale_x_discrete(labels = rep("", nlevels(df$.BarID))) +
        coord_cartesian(clip = "off")
    }

    p + facet_wrap(~Feature, scales = "free_y", ncol = ncol_value) +
      theme_bw(base_size = scaled_base_size) +
      labs(title = title_text, x = "", y = y_axis_title) +
      theme(
        axis.text.x = element_text(
          angle = if (isTRUE(input$rotate_feature_labels)) 90 else 0,
          hjust = 1,
          vjust = 0.5,
          size = scaled_feature_text_size
        ),
        axis.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(size = input$feature_title_size, face = "bold", hjust = 0.5),
        legend.title = element_text(face = "bold"),
        plot.margin = margin(10, 10, 35, 10),
        strip.text = element_text(face = "bold")
      )
  })

  output$feature_plot_note <- renderText({
    order_columns <- input$feature_order_columns
    order_columns <- order_columns[!is.na(order_columns) & nzchar(order_columns)]
    order_note <- if (length(order_columns) > 0) {
      paste0(" Samples ordered by: ", paste(order_columns, collapse = " -> "), ".")
    } else {
      " Samples ordered by current data/metadata order."
    }
    scale_note <- if (identical(input$feature_value_scale, "log2")) {
      " Plot scale: log2 abundance."
    } else {
      " Plot scale: raw abundance."
    }
    interactive_note <- if (isTRUE(input$feature_interactive)) {
      " Interactive zoom/pan is enabled; use the plot toolbar to reset axes."
    } else {
      ""
    }
    paste0(
      "Feature source: ",
      switch(
        input$feature_data_source,
        "no_impute" = "Table S2 non-imputed protein report.",
        "S3_batch_corrected" = "Table S3 batch-corrected protein report (values are corrected log2 abundance).",
        "Table S3 imputed protein report."
      ),
      scale_note,
      order_note,
      interactive_note,
      " Select a row in a protein table or the log2 abundance matrix to open that protein here."
    )
  })

  output$feature_matrix_table <- renderDT({
    data <- feature_matrix_data()
    abundance_cols <- colnames(data)[seq.int(5, ncol(data))]
    z_scores <- attr(data, "z_scores")
    validate(need(!is.null(z_scores), "Feature matrix z-score styling values are unavailable."))
    z_cols <- paste0(".z_", seq_along(abundance_cols))
    colnames(z_scores) <- z_cols
    styled_data <- data.frame(data, z_scores, check.names = FALSE)
    hidden_targets <- seq.int(ncol(data), ncol(styled_data) - 1)
    z_breaks <- c(-2, -1, -0.01, 0.01, 1, 2)
    z_colors <- c("#2166AC", "#67A9CF", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#EF8A62", "#B2182B")

    table <- datatable(
      styled_data,
      rownames = FALSE,
      selection = "single",
      filter = "top",
      options = list(
        scrollX = TRUE,
        pageLength = 25,
        deferRender = TRUE,
        order = list(list(0, "desc")),
        columnDefs = list(list(targets = hidden_targets, visible = FALSE, searchable = FALSE, orderable = FALSE))
      )
    )
    table <- formatStyle(
      table,
      "Significant",
      backgroundColor = styleEqual(c("Yes: Increased", "Yes: Decreased", "No"), c("#F4B4B4", "#B4D3EC", "#FFFFFF")),
      fontWeight = styleEqual(c("Yes: Increased", "Yes: Decreased", "No"), c("bold", "bold", "normal"))
    )
    table <- formatRound(table, abundance_cols, digits = 1)
    for (index in seq_along(abundance_cols)) {
      table <- formatStyle(
        table,
        abundance_cols[index],
        valueColumns = z_cols[index],
        backgroundColor = styleInterval(z_breaks, z_colors)
      )
    }
    table
  }, server = TRUE)

  output$annotation_table <- renderDT({
    ann <- detected_annotations()
    req(ann)
    datatable(ann, options = list(scrollX = TRUE, pageLength = 10))
  })

  output$metadata_build_note <- renderText({
    built <- draft_metadata()
    built <- apply_proteomics_metadata_cell_edits(built, spqc_metadata_draft_edits())
    if (is.null(project_file("meta_file")) && is.null(project_file("run_order_file")) && is.null(project_file("sample_details_file"))) {
      return(paste0(
        nrow(built),
        " samples detected. Upload the condition setup, run-order table, and sample-details workbook to populate analysis metadata."
      ))
    }

    notes <- paste0(nrow(built), " samples detected.")
    if (!is.null(project_file("meta_file"))) {
      md <- metadata_df()
      notes <- paste0(notes, " ", nrow(md), " condition rows loaded as the metadata basis.")
    }
    if (!is.null(project_file("run_order_file"))) {
      notes <- paste0(notes, " ", sum(!is.na(built$RunOrder)), " matched run-order rows.")
    }
    if (!is.null(project_file("sample_details_file"))) {
      notes <- paste0(notes, " ", sum(!is.na(built$SampleName)), " received submitted sample names.")
    }
    spqc_rows <- grepl("SPQC", built$Sample, ignore.case = TRUE) |
      ("SampleName" %in% colnames(built) & grepl("SPQC", built$SampleName, ignore.case = TRUE))
    spqc_rows[is.na(spqc_rows)] <- FALSE
    if (any(spqc_rows)) {
      spqc_groups <- unique(as.character(built$Condition[spqc_rows]))
      spqc_groups <- spqc_groups[!is.na(spqc_groups) & nzchar(spqc_groups)]
      override_rules <- parse_spqc_batch_overrides(input$spqc_batch_overrides)
      notes <- paste0(
        notes,
        " ",
        sum(spqc_rows),
        " SPQC rows assigned to ",
        if (length(spqc_groups)) paste(spqc_groups, collapse = ", ") else "existing metadata groups",
        ".",
        if (nrow(override_rules) > 0) paste0(" SPQC batch override rules: ", nrow(override_rules), ".") else ""
      )
    }
    if ("Excluded" %in% colnames(built)) {
      notes <- paste0(notes, " ", sum(as.logical(built$Excluded), na.rm = TRUE), " samples excluded from analysis.")
    }
    notes
  })

  output$sample_exclusion_note <- renderText({
    terms <- parse_sample_exclusion_text(input$sample_exclusions_text)
    if (length(terms) == 0) return("No global sample exclusions entered.")
    built <- draft_metadata()
    built <- apply_proteomics_metadata_cell_edits(built, spqc_metadata_draft_edits())
    excluded_count <- if ("Excluded" %in% colnames(built)) sum(as.logical(built$Excluded), na.rm = TRUE) else 0
    paste0(
      length(terms),
      " exclusion term(s) entered; ",
      excluded_count,
      " metadata sample(s) currently matched. Exclusions apply to PCA, stats, CV plots, batch correction, feature plots, boxplots, and protein measurement exports."
    )
  })

  output$sample_exclusion_preview <- renderDT({
    built <- draft_metadata()
    built <- apply_proteomics_metadata_cell_edits(built, spqc_metadata_draft_edits())
    if (!"Excluded" %in% colnames(built)) return(datatable(data.frame(Message = "No metadata loaded."), options = list(dom = "t")))
    excluded <- as.logical(built$Excluded)
    excluded[is.na(excluded)] <- FALSE
    preview_cols <- intersect(c("Sample", "SampleName", "Condition", "Batch", "AnalysisLabel", "Run Label", "File Name", "ExclusionReason"), colnames(built))
    out <- built[excluded, preview_cols, drop = FALSE]
    if (nrow(out) == 0) out <- data.frame(Message = "No samples currently matched.", stringsAsFactors = FALSE)
    datatable(out, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5))
  })

  output$project_bundle_note <- renderText({
    project_bundle_message()
  })

  output$project_files_table <- renderDT({
    datatable(
      loaded_proteomics_project_status(project_file_status()),
      rownames = FALSE,
      options = list(dom = "t", paging = FALSE, scrollX = TRUE)
    )
  })

  output$metadata_replacement_note <- renderText(metadata_replacement_message())

  active_file_note <- function(ids) {
    rows <- project_file_status()
    rows <- rows[rows$InputID %in% ids, , drop = FALSE]
    loaded <- rows[rows$Status != "Not loaded", , drop = FALSE]
    if (nrow(loaded) == 0) return("Active restored/uploaded files: none yet.")
    paste(
      paste0(loaded$FileType, ": ", loaded$FileName, " (", loaded$Status, ")"),
      collapse = "\n"
    )
  }

  output$metadata_loaded_files_note <- renderText({
    active_file_note(c("meta_file", "run_order_file", "sample_details_file"))
  })

  output$protein_loaded_files_note <- renderText({
    active_file_note(c("protein_no_impute_file", "protein_imputed_file"))
  })

  output$condition_setup_note <- renderText({
    headers <- condition_setup_headers()
    suffix <- input$condition_setup_run_label_suffix
    if (is.null(suffix) || !nzchar(suffix)) suffix <- condition_setup_inferred_suffix()
    detail_text <- if (is.null(project_file("condition_setup_sample_details_file"))) {
      "Select an order sample details workbook to build the condition setup table."
    } else {
      details <- condition_setup_sample_details()
      condition_col <- input$condition_setup_condition_col
      replicate_col <- input$condition_setup_replicate_order_col
      inferred_col <- infer_condition_setup_label_col(details)
      resolved_condition_col <- if (is.null(condition_col) || identical(condition_col, "__infer__") || !condition_col %in% colnames(details)) {
        if (!is.na(inferred_col)) inferred_col else "(no suitable column found)"
      } else {
        condition_col
      }
      parsed_note <- ""
      if (!is.null(resolved_condition_col) && resolved_condition_col %in% colnames(details)) {
        parsed <- parse_condition_setup_replicate(details[[resolved_condition_col]])
        parsed_n <- sum(!is.na(parsed$Replicate) & nzchar(parsed$Label))
        if (parsed_n > 0) {
          parsed_note <- paste0(" Embedded replicate numbers detected in ", resolved_condition_col, " for ", parsed_n, " rows; those numbers will be used in the Replicate column.")
        }
      }
      paste0(
        nrow(details), " samples detected from sample details. Condition/Label column: ",
        resolved_condition_col,
        if (!is.na(inferred_col) && identical(condition_col, "__infer__")) paste0(" (inferred from ", inferred_col, ")") else "",
        "; replicate order column: ",
        if (!is.null(replicate_col) && nzchar(replicate_col)) replicate_col else "(select a column)",
        ".",
        parsed_note
      )
    }
    paste0(
      "Condition setup headers: ", paste(headers, collapse = " | "), "\n",
      "Run label suffix: ", if (nzchar(suffix)) suffix else "(blank; upload a template or type a suffix)", "\n",
      detail_text
    )
  })

  output$condition_setup_sample_details_preview <- renderDT({
    datatable(condition_setup_sample_details(), options = list(scrollX = TRUE, pageLength = 15))
  })

  output$condition_setup_preview <- renderDT({
    datatable(condition_setup_table(), options = list(scrollX = TRUE, pageLength = 15))
  })

  output$evosep_manifest_preview <- renderDT({
    datatable(evosep_manifest_details(), options = list(scrollX = TRUE, pageLength = 15))
  })

  output$evosep_queue_note <- renderText({
    if (is.null(project_file("evosep_manifest_file")) || is.null(project_file("evosep_template_file"))) {
      return("Select both an order sample details workbook and an Evosep .csl template.")
    }
    manifest <- evosep_manifest_details()
    plate <- evosep_assigned_plate()
    sample_positions <- parse_position_list(input$evosep_sample_positions)
    available_sample_positions <- setdiff(sample_positions, c(parse_position_list(input$evosep_adh_positions), parse_position_list(input$evosep_spqc_positions)))
    assigned_samples <- sum(plate$Assignment == "Sample")
    paste0(
      nrow(manifest), " manifest samples detected. ",
      assigned_samples, " samples assigned to source vial positions. ",
      sum(plate$Assignment == "ADH"), " ADH rows and ",
      sum(plate$Assignment == "SPQC"), " SPQC rows included. ",
      if (length(available_sample_positions) < nrow(manifest)) {
        paste0("Only ", length(available_sample_positions), " sample positions are available; remaining manifest samples are not queued.")
      } else {
        "All manifest samples fit into the selected sample positions."
      }
    )
  })

  output$evosep_plate_grid <- renderUI({
    plate <- evosep_assigned_plate()
    cell_style <- function(assignment) {
      color <- switch(
        assignment,
        "ADH" = "#FFF2CC",
        "SPQC" = "#FCE4D6",
        "Sample" = "#D9EAD3",
        "#F7F7F7"
      )
      paste0("border:1px solid #CCCCCC; padding:6px; min-width:88px; height:58px; background:", color, "; vertical-align:top; font-size:12px;")
    }
    rows <- lapply(LETTERS[1:8], function(row_label) {
      row_data <- plate[plate$Row == row_label, , drop = FALSE]
      tags$tr(c(
        list(tags$th(style = "padding:6px; text-align:center;", row_label)),
        lapply(seq_len(nrow(row_data)), function(i) {
          label <- if (nzchar(row_data$Assignment[i])) {
            paste0(row_data$Assignment[i], if (nzchar(row_data$ManifestID[i])) paste0(": ", row_data$ManifestID[i]) else "")
          } else {
            "Empty"
          }
          tags$td(
            style = cell_style(row_data$Assignment[i]),
            tags$strong(row_data$Position[i]),
            tags$br(),
            tags$span(label)
          )
        })
      ))
    })
    header <- tags$tr(c(list(tags$th("")), lapply(seq_len(12), function(i) tags$th(style = "padding:6px; text-align:center;", i))))
    tags$table(style = "border-collapse:collapse; margin-bottom:12px;", c(list(header), rows))
  })

  output$evosep_assignment_table <- renderDT({
    plate <- evosep_assigned_plate()
    datatable(
      plate[, c("Position", "Well", "Assignment", "ManifestID", "SampleName", "Filename"), drop = FALSE],
      rownames = FALSE,
      editable = list(target = "cell", disable = list(columns = c(0, 1))),
      filter = "top",
      options = list(
        scrollX = TRUE,
        pageLength = 24,
        order = list(list(0, "asc")),
        columnDefs = list(
          list(targets = c(0, 1), className = "dt-body-center"),
          list(targets = 2, className = "dt-body-center")
        )
      )
    )
  })

  output$evosep_queue_preview <- renderDT({
    preview <- evosep_queue_table()
    datatable(
      preview[, c("RunNumber", "SourceTray", "SourceVial", "SampleName", "XcaliburFilename", "OutputDir", "Comment", "Assignment", "ManifestID", "ManifestSampleName"), drop = FALSE],
      rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 25)
    )
  })

  locked_metadata_columns <- c("Sample", "RunOrder", "Run Order", "Filename", "FileName", "File Name")

  output$metadata_preview <- renderDT({
    metadata_editor_revision()
    table <- draft_metadata()
    table <- apply_proteomics_metadata_cell_edits(table, isolate(spqc_metadata_draft_edits()))
    selected <- input$metadata_export_columns
    view <- prepare_metadata_editor_view(table, selected, locked_metadata_columns)
    datatable(
      view$display,
      rownames = FALSE,
      editable = list(target = "cell", disable = list(columns = view$locked_indices)),
      options = list(scrollX = TRUE, pageLength = 15, stateSave = TRUE)
    )
  }, server = TRUE)

  output$metadata_apply_status <- renderText({
    count <- nrow(spqc_metadata_draft_edits())
    paste(
      metadata_apply_message(),
      if (count > 0L) paste0(count, " draft metadata cell change", if (count == 1L) "." else "s.") else NULL
    )
  })

  observeEvent(list(
    input$meta_file,
    input$run_order_file,
    input$sample_details_file,
    input$sample_exclusions_text,
    input$sample_exclusion_reason,
    input$spqc_assignment_mode,
    input$spqc_group_label,
    input$spqc_group_prefix,
    input$spqc_batch_overrides,
    input$metadata_export_columns
  ), {
    if (isTRUE(metadata_draft_tracking_enabled())) metadata_apply_message("Unapplied metadata changes.")
  }, ignoreInit = TRUE)

  observeEvent(input$metadata_preview_cell_edit, {
    info <- input$metadata_preview_cell_edit
    draft_edits <- isolate(spqc_metadata_draft_edits())
    table <- isolate(draft_metadata())
    table <- apply_proteomics_metadata_cell_edits(table, draft_edits)
    view <- prepare_metadata_editor_view(table, isolate(input$metadata_export_columns), locked_metadata_columns)
    table <- view$display
    column_index <- as.integer(info$col) + 1L
    if (!nrow(table) || info$row < 1 || info$row > nrow(table) || column_index < 1 || column_index > ncol(table)) return()
    column <- colnames(table)[column_index]
    if (column %in% locked_metadata_columns) return()
    sample <- view$sample_keys[info$row]
    draft_edits <- draft_edits[!(as.character(draft_edits$Sample) == sample & as.character(draft_edits$Column) == column), , drop = FALSE]
    draft_edits <- rbind(draft_edits, data.frame(Sample = sample, Column = column, Value = normalize_proteomics_text(info$value), stringsAsFactors = FALSE))
    spqc_metadata_draft_edits(draft_edits)
    metadata_apply_message("Unapplied metadata changes.")
  })

  observeEvent(input$apply_metadata_changes, {
    tryCatch({
      draft <- draft_metadata()
      columns <- input$metadata_export_columns
      if (is.null(columns) || !length(columns)) columns <- metadata_default_columns(colnames(draft))
      candidate <- build_proteomics_applied_state(
        metadata = draft,
        columns = columns,
        committed_edits = if (isTRUE(spqc_clear_pending())) empty_spqc_metadata_edits() else spqc_metadata_edits(),
        draft_edits = spqc_metadata_draft_edits()
      )
      cache <- invalidate_proteomics_project_cache(project_db_cache(), "meta_file")
      cache$metadata <- candidate$metadata
      project_db_cache(cache)
      cached_batch_corrected_s3_result(NULL)
      spqc_metadata_edits(candidate$spqc_metadata_edits)
      spqc_metadata_draft_edits(empty_spqc_metadata_edits())
      spqc_clear_pending(FALSE)
      applied_metadata_state(candidate)
      metadata_apply_revision(isolate(metadata_apply_revision()) + 1L)
      metadata_apply_message("All metadata changes applied.")
      autosave_active_project("applied metadata changes", include_derived = FALSE)
      showNotification("Metadata changes applied to all tabs.", type = "message")
    }, error = function(e) {
      metadata_apply_message(paste0("Metadata was not applied: ", conditionMessage(e)))
      showNotification(conditionMessage(e), type = "error", duration = NULL)
    })
  })

  observeEvent(input$discard_metadata_edits, {
    spqc_metadata_draft_edits(empty_spqc_metadata_edits())
    spqc_clear_pending(FALSE)
    metadata_editor_revision(isolate(metadata_editor_revision()) + 1L)
    metadata_apply_message("Draft metadata edits discarded.")
  })

  output$protein_table_note <- renderText({
    notes <- character(0)
    if (protein_source_available("S2")) {
      s2 <- protein_no_impute_table()
      notes <- c(notes, paste0(
        "Table S2: ", nrow(s2), " proteins; ",
          attr(s2, "matched_precursor"), " precursor headers and ",
          attr(s2, "matched_abundance"), " abundance headers renamed; %CV: ",
          paste(attr(s2, "cv_conditions"), collapse = ", "), ".",
          if (length(attr(s2, "stats_methods")) > 0) paste0(" Statistics: ", paste(attr(s2, "stats_methods"), collapse = "; "), ".") else ""
      ))
    }
    if (protein_source_available("S3")) {
      s3 <- protein_imputed_table()
      notes <- c(notes, paste0(
        "Table S3: ", nrow(s3), " proteins; ",
          attr(s3, "matched_precursor"), " precursor headers and ",
          attr(s3, "matched_abundance"), " abundance headers renamed; %CV: ",
          paste(attr(s3, "cv_conditions"), collapse = ", "), ".",
          if (length(attr(s3, "stats_methods")) > 0) paste0(" Statistics: ", paste(attr(s3, "stats_methods"), collapse = "; "), ".") else ""
      ))
    }
    s3bc <- tryCatch(batch_corrected_s3_table(), error = function(e) NULL)
    if (!is.null(s3bc)) {
      notes <- c(notes, paste0(
        "Batch-corrected S3: ", nrow(s3bc), " proteins; ",
        attr(s3bc, "matched_abundance"), " batch-corrected abundance columns.",
        if (length(attr(s3bc, "stats_methods")) > 0) paste0(" Statistics: ", paste(attr(s3bc, "stats_methods"), collapse = "; "), ".") else ""
      ))
    }
    if (length(notes) == 0) "Select one or both protein report files to add protein tables to the Excel workbook." else paste(notes, collapse = "\n")
  })

  datatable_schema_key <- function(data) {
    column_text <- paste(colnames(data), collapse = "|")
    paste0(nrow(data), "_", ncol(data), "_", sum(nchar(column_text)), "_", sum(utf8ToInt(column_text)))
  }

  output$protein_no_impute_preview_ui <- renderUI({
    req(protein_source_available("S2"))
    table <- protein_no_impute_table()
    tags$div(
      id = paste0("protein_no_impute_preview_container_", datatable_schema_key(table)),
      DTOutput("protein_no_impute_preview")
    )
  })

  output$protein_imputed_preview_ui <- renderUI({
    req(protein_source_available("S3"))
    table <- protein_imputed_table()
    tags$div(
      id = paste0("protein_imputed_preview_container_", datatable_schema_key(table)),
      DTOutput("protein_imputed_preview")
    )
  })

  output$protein_no_impute_preview <- renderDT({
    req(protein_source_available("S2"))
    datatable(
      protein_no_impute_table(),
      rownames = FALSE,
      selection = "single",
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        deferRender = TRUE,
        processing = TRUE,
        stateSave = FALSE,
        searchDelay = 500
      )
    )
  }, server = TRUE)

  output$protein_imputed_preview <- renderDT({
    req(protein_source_available("S3"))
    datatable(
      protein_imputed_table(),
      rownames = FALSE,
      selection = "single",
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        deferRender = TRUE,
        processing = TRUE,
        stateSave = FALSE,
        searchDelay = 500
      )
    )
  }, server = TRUE)

  output$batch_cache_note <- renderText({
    batch_cache_message()
  })

  output$batch_correction_note <- renderText({
    result <- tryCatch(batch_corrected_s3_result(), error = function(e) NULL)
    if (is.null(result)) {
      return("Upload Table S3 and metadata, choose a batch column, then click Run batch correction. The app uses direct sva::ComBat when available; HarmonizR is optional.")
    }
    prepared <- result$prepared
    batch_counts <- table(prepared$sample_map$Batch)
    filter_summary <- prepared$filter_summary
    paste0(
      "Batch correction complete using ", result$correction_method, " mode ", result$combat_mode, ".\n",
      nrow(prepared$data_as_input), " proteins and ", ncol(prepared$data_as_input), " samples corrected on log2 abundance values.\n",
      "Prefilter kept ", filter_summary$kept_rows, " of ", filter_summary$source_rows,
      " proteins present in at least ", filter_summary$min_batches_for_feature, " batches; removed ",
      filter_summary$removed_rows, ".\n",
      "Batches: ", paste(paste(names(batch_counts), batch_counts, sep = "="), collapse = "; "), ".\n",
      prepared$confounding$message
    )
  })

  output$batch_corrected_s3_preview <- renderDT({
    datatable(
      batch_corrected_s3_table(),
      rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 10, deferRender = TRUE, processing = TRUE)
    )
  }, server = TRUE)

  output$cv_plot <- renderPlot({
    cv_plot_obj()
  }, res = 120)

  preview_height <- function(height_inches, minimum = 300) {
    if (is.null(height_inches) || !is.finite(height_inches)) return(paste0(minimum, "px"))
    paste0(max(minimum, round(height_inches * 90)), "px")
  }

  panel_row_count <- function(n_panels, ncol_value) {
    n_panels <- suppressWarnings(as.integer(n_panels))
    if (is.na(n_panels) || n_panels < 1) n_panels <- 1L
    ncol_value <- suppressWarnings(as.integer(ncol_value))
    if (is.na(ncol_value) || ncol_value < 1) ncol_value <- 1L
    ceiling(n_panels / ncol_value)
  }

  output$cv_plot_ui <- renderUI({
    plotOutput("cv_plot", height = preview_height(input$cv_figure_height, 300))
  })

  output$pca_plot_ui <- renderUI({
    plotOutput("pca_plot", height = preview_height(input$pca_figure_height, 400))
  })

  output$clustvis_pca_plot_ui <- renderUI({
    if (isTRUE(input$clustvis_pca_interactive)) {
      plotly::plotlyOutput("clustvis_pca_plot_interactive", height = preview_height(input$clustvis_pca_height, 400))
    } else {
      plotOutput("clustvis_pca_plot", height = preview_height(input$clustvis_pca_height, 400))
    }
  })

  output$gene_set_enrichment_plot_ui <- renderUI({
    plotOutput("gene_set_enrichment_plot", height = preview_height(input$gsea_plot_height, 350))
  })

  output$volcano_plot_ui <- renderUI({
    if (isTRUE(input$volcano_interactive)) {
      plotly::plotlyOutput("volcano_plot_interactive", height = preview_height(input$volcano_figure_height, 400))
    } else {
      plotOutput("volcano_plot", height = preview_height(input$volcano_figure_height, 400))
    }
  })

  output$feature_plot_ui <- renderUI({
    n_features <- length(unique(as.character(input$feature_select)))
    rows <- panel_row_count(max(n_features, 1L), input$feature_plot_ncol)
    min_height <- 260 * rows
    if (isTRUE(input$feature_interactive)) {
      plotly::plotlyOutput("feature_plot_interactive", height = preview_height(input$feature_figure_height, min_height))
    } else {
      plotOutput("feature_plot", height = preview_height(input$feature_figure_height, min_height))
    }
  })

  output$script_box_plot_ui <- renderUI({
    n_features <- length(unique(as.character(input$script_box_features)))
    rows <- panel_row_count(max(n_features, 1L), input$script_box_ncol)
    min_height <- 260 * rows
    plotOutput("script_box_plot", height = preview_height(input$script_box_height, min_height))
  })

  output$identification_overview_plot_ui <- renderUI({
    plotOutput("identification_overview_plot", height = preview_height(input$identification_overview_height, 300))
  })

  output$run_identifications_plot_ui <- renderUI({
    plotOutput("run_identifications_plot", height = preview_height(input$run_identifications_height, 300))
  })

  output$identification_overview_plot <- renderPlot({
    identification_overview_plot_obj()
  }, res = 120)

  output$run_identifications_plot <- renderPlot({
    run_identifications_plot_obj()
  }, res = 120)

  output$run_identifications_note <- renderText({
    overview <- identifications_overview_data()
    stacked <- run_identifications_data()
    metric_label <- if (input$identification_metric == "ProteinGroups") "protein groups" else "precursors"
    paste0(
      nrow(overview),
      " runs loaded. Plotting ",
      metric_label,
      "; stacked table total = ",
      format(sum(stacked$Value, na.rm = TRUE), big.mark = ","),
      " identifications."
    )
  })

  output$cv_plot_note <- renderText({
    plot_data <- cv_plot_data()
    median_text <- paste0(
      plot_data$medians$Condition, ": ",
      format(round(plot_data$medians$MedianCV, 1), nsmall = 1), "%"
    )
    paste0(plot_data$mode, ". Median %CV: ", paste(median_text, collapse = "; "), ".")
  })

  output$pca_data_note <- renderText({
    result <- pca_results()
    paste0(
      "Retained ", result$retained_features, " of ", result$total_features,
      " proteins at the ", input$pca_min_observed_percent,
      "% observed-samples threshold; ", result$missing_values_before_imputation,
      " missing values imputed before PCA. PCA used ncp = ", result$ncp_used,
      " and ", result$variable_features, " variable proteins."
    )
  })

  output$clustvis_pca_note <- renderText({
    result <- clustvis_pca_results()
    source_label <- switch(
      input$clustvis_pca_source,
      "no_impute" = "Table S2 non-imputed protein report",
      "S3_batch_corrected" = "Table S3 batch-corrected protein report (PCA uses corrected log2 values)",
      "Table S3 imputed protein report"
    )
    paste0(
      source_label,
      ". PCA sample subset included ", result$included_samples, " sample(s) and excluded ",
      result$excluded_by_subset, " sample(s)",
      ". Retained ", result$retained_features, " of ", result$total_features,
      " proteins at the ", input$clustvis_pca_min_observed_percent,
      "% observed-samples threshold; ", result$missing_values_before_imputation,
      " missing values handled by ", result$imputation_method,
      "; ", result$variable_features, " variable proteins used for prcomp."
    )
  })

  output$script_box_note <- renderText({
    df <- script_box_plot_df()
    style_note <- if (identical(input$script_box_plot_style, "mean_sd")) {
      "Bars show mean +/- standard deviation; points are individual samples."
    } else {
      "Boxes show median and interquartile range; whiskers use the ggplot2 boxplot rule; points are individual samples."
    }
    source_label <- switch(
      input$script_box_source,
      "no_impute" = "Table S2 non-imputed protein report",
      "S3_batch_corrected" = "Table S3 batch-corrected protein report (corrected log2 abundance)",
      "Table S3 imputed protein report"
    )
    selected_features <- unique(as.character(df$Feature))
    paste0(
      "Faceted boxplot for ", length(selected_features), " protein group", if (length(selected_features) == 1L) "" else "s",
      " from ", source_label,
      " using ", length(unique(df$Sample)), " samples across ",
      length(unique(as.character(df$GroupValue))), " groups. ",
      style_note, " ",
      if (requireNamespace("ggprism", quietly = TRUE)) "Using ggprism theme." else "Using theme_classic because ggprism is not installed."
    )
  })

  output$correlation_note <- renderText({
    data <- correlation_results_data()
    valid <- sum(is.finite(data$PValue), na.rm = TRUE)
    significant <- sum(is.finite(data$QValue) & data$QValue < 0.05, na.rm = TRUE)
    source_label <- switch(
      input$correlation_source,
      "no_impute" = "Table S2 non-imputed protein report",
      "S3_batch_corrected" = "Table S3 batch-corrected protein report",
      "Table S3 imputed protein report"
    )
    scale_label <- if (identical(input$correlation_value_scale, "log2")) "log2 abundance" else "raw abundance"
    method_label <- switch(
      input$correlation_method,
      pearson = "Pearson correlation",
      spearman = "Spearman correlation",
      "linear regression"
    )
    group_filter <- if ("GroupFilter" %in% colnames(data) && nrow(data) > 0) data$GroupFilter[1] else "all samples"
    covariate_text <- if ("Covariates" %in% colnames(data) && nrow(data) > 0 && nzchar(data$Covariates[1])) {
      paste0(" Adjusted for: ", data$Covariates[1], ".")
    } else {
      ""
    }
    paste0(
      source_label,
      " using ",
      scale_label,
      " and ",
      method_label,
      " across ",
      group_filter,
      ".",
      covariate_text,
      " ",
      valid,
      " proteins had enough finite values for comparison against ",
      input$correlation_feature,
      "; ",
      significant,
      " have BH q-value < 0.05."
    )
  })

  output$download_destination_note <- renderText({
    destination_dir <- download_destination_path()
    destination_text <- if (is.na(destination_dir) || !nzchar(destination_dir)) "(no folder selected)" else destination_dir
    paste0(
      "Destination folder: ", destination_text, "\n",
      download_destination_message()
    )
  })

    output$pca_plot <- renderPlot({ pca_plot_obj() }, res = 120)
    output$clustvis_pca_plot <- renderPlot({ clustvis_pca_plot_obj() }, res = 120)
    output$clustvis_pca_plot_interactive <- plotly::renderPlotly({
      validate(need(requireNamespace("plotly", quietly = TRUE), "Install the R package 'plotly' to use interactive PCA."))
      p_obj <- plotly::ggplotly(clustvis_pca_plot_obj(), tooltip = "text")
      color_var <- if (!is.null(input$clustvis_pca_color_by) && input$clustvis_pca_color_by != "None") input$clustvis_pca_color_by else NULL
      shape_var <- if (!is.null(input$clustvis_pca_shape_by) && input$clustvis_pca_shape_by != "None") input$clustvis_pca_shape_by else NULL
      split_plotly_pca_legend(p_obj, color_var, shape_var)
    })
    output$pca_loadings_plot <- renderPlot({ pca_loadings_plot_obj() }, res = 120)
    output$volcano_plot <- renderPlot({ volcano_plot_obj() }, res = 120)
    output$volcano_plot_interactive <- plotly::renderPlotly({ volcano_interactive_obj() })
    output$gene_set_enrichment_plot <- renderPlot({ gene_set_enrichment_plot_obj() }, res = 120)
    output$feature_plot <- renderPlot({ req(input$feature_select); feature_plot_obj() }, res = 120)
    output$feature_plot_interactive <- plotly::renderPlotly({
      req(input$feature_select)
      validate(need(requireNamespace("plotly", quietly = TRUE), "Install the R package 'plotly' to use the interactive feature plot."))
      plotly::ggplotly(feature_plot_obj(), tooltip = "text") %>%
        plotly::layout(dragmode = "zoom")
    })
    output$script_box_plot <- renderPlot({ req(input$script_box_features); script_box_plot_obj() }, res = 120)
    output$correlation_lollipop_plot <- renderPlot({ req(input$correlation_feature); correlation_lollipop_plot_obj() }, res = 120)

    output$correlation_lollipop_plot_ui <- renderUI({
      height_px <- max(250, round(input$correlation_figure_height * 96))
      plotOutput("correlation_lollipop_plot", height = paste0(height_px, "px"))
    })

    output$scores_table <- renderDT({
      datatable(pca_results()$scores, options = list(scrollX = TRUE, pageLength = 10))
    })

    output$clustvis_pca_scores_table <- renderDT({
      datatable(clustvis_pca_results()$scores, options = list(scrollX = TRUE, pageLength = 10))
    })

    output$pca_loadings_table <- renderDT({
      datatable(pca_loadings_data(), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10))
    })

    output$gene_set_enrichment_note <- renderText({
      data <- gene_set_enrichment_data()
      paste0(
        nrow(data),
        " gene sets tested from ",
        attr(data, "gene_sets"),
        " loaded sets using ",
        attr(data, "ranked_genes"),
        " ranked genes."
      )
    })

    output$gene_set_enrichment_table <- renderDT({
      data <- gene_set_enrichment_data()
      preferred <- intersect(c("GeneSet", "Direction", "NES", "pval", "BH_FDR", "size", "leadingEdge"), colnames(data))
      datatable(data[, preferred, drop = FALSE], rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10))
    })

    output$correlation_results_table <- renderDT({
      data <- correlation_results_data()
      preferred <- intersect(
        c("Protein", "DisplayProtein", "Method", "Beta", "R2", "PValue", "QValue", "Correlation", "N",
          "GroupFilter", "Covariates",
          "PG.Genes", "PG.ProteinNames", "PG.ProteinDescriptions", "PG.ProteinGroups", "RawFeatureID"),
        colnames(data)
      )
      datatable(
        data[, preferred, drop = FALSE],
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 15)
      )
    })

    output$volcano_note <- renderText({
      data <- volcano_plot_data()
      source_label <- switch(
        input$volcano_source,
        "S2" = "Table S2. Protein, no impute",
        "S3_batch_corrected" = "Table S3. Protein, imputed, batch-corrected",
        "Table S3. Protein, imputed"
      )
      paste0(
        source_label,
        ": ",
        sum(data$Status == "Increased", na.rm = TRUE),
        " increased and ",
        sum(data$Status == "Decreased", na.rm = TRUE),
        " decreased proteins pass the selected thresholds."
      )
    })

    output$volcano_hits_table <- renderDT({
      data <- volcano_hits_data()
      display_columns <- volcano_hits_display_columns(
        c("Protein", "ProteinName", "ProteinDescription", "Status", "Log2FoldChange", "PValue", "BH_FDR"),
        input$stats_bh_fdr
      )
      datatable(
        data[, display_columns, drop = FALSE],
        rownames = FALSE,
        selection = list(mode = "multiple", selected = NULL, target = "row"),
        options = list(scrollX = TRUE, pageLength = 10)
      )
    })

    output$volcano_counts_table <- renderDT({
      data <- volcano_counts_data()
      datatable(
        data[, setdiff(colnames(data), "ComparisonID"), drop = FALSE],
        rownames = FALSE,
        selection = "single",
        options = list(dom = "t", paging = FALSE, ordering = TRUE, scrollX = TRUE)
      )
    })

  output$download_scores <- downloadHandler(
    filename = function() build_export_filename("pca_scores.csv"),
    content = function(file) {
      out_name <- build_export_filename("pca_scores.csv")
      write.csv(pca_results()$scores, file, row.names = FALSE)
      save_download_copy(file, out_name)
    }
  )

  output$download_built_metadata <- downloadHandler(
    filename = function() build_export_filename("Table_S1_Metadata.csv"),
    content = function(file) {
      out_name <- build_export_filename("Table_S1_Metadata.csv")
      write.csv(exported_metadata(), file, row.names = FALSE, na = "NaN")
      save_download_copy(file, out_name)
    }
  )

  output$download_active_metadata <- downloadHandler(
    filename = function() build_export_filename("active_project_metadata.csv"),
    content = function(file) {
      write.csv(built_metadata(), file, row.names = FALSE, na = "")
      save_download_copy(file, build_export_filename("active_project_metadata.csv"))
    }
  )

  output$download_protein_s2_csv <- downloadHandler(
    filename = function() build_export_filename("Table_S2_Protein_no_impute.csv"),
    content = function(file) {
      out_name <- build_export_filename("Table_S2_Protein_no_impute.csv")
      write.csv(protein_no_impute_table(), file, row.names = FALSE, na = "NaN")
      save_download_copy(file, out_name)
    }
  )

  output$download_protein_s3_csv <- downloadHandler(
    filename = function() build_export_filename("Table_S3_Protein_imputed.csv"),
    content = function(file) {
      out_name <- build_export_filename("Table_S3_Protein_imputed.csv")
      write.csv(protein_imputed_table(), file, row.names = FALSE, na = "NaN")
      save_download_copy(file, out_name)
    }
  )

  output$download_batch_corrected_s3_csv <- downloadHandler(
    filename = function() build_export_filename("Table_S3_Protein_imputed_batch_corrected.csv"),
    content = function(file) {
      out_name <- build_export_filename("Table_S3_Protein_imputed_batch_corrected.csv")
      write.csv(batch_corrected_s3_table(), file, row.names = FALSE, na = "NaN")
      save_download_copy(file, out_name)
    }
  )

  output$download_condition_setup_tsv <- downloadHandler(
    filename = function() {
      template_file <- project_file("condition_setup_template_file")
      source_name <- if (!is.null(template_file) && !is.null(template_file$name) && nzchar(template_file$name)) {
        template_file$name
      } else {
        "ConditionSetup.tsv"
      }
      ext <- tools::file_ext(source_name)
      stem <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", source_name) else source_name
      out_name <- if (nzchar(ext)) paste0(stem, "_app.", ext) else paste0(stem, "_app.tsv")
      build_export_filename(out_name)
    },
    content = function(file) {
      template_file <- project_file("condition_setup_template_file")
      source_name <- if (!is.null(template_file) && !is.null(template_file$name) && nzchar(template_file$name)) {
        template_file$name
      } else {
        "ConditionSetup.tsv"
      }
      ext <- tools::file_ext(source_name)
      stem <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", source_name) else source_name
      out_name <- if (nzchar(ext)) paste0(stem, "_app.", ext) else paste0(stem, "_app.tsv")
      out_name <- build_export_filename(out_name)
      write.table(condition_setup_table(), file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
      save_download_copy(file, out_name)
    }
  )

  output$download_evosep_csl <- downloadHandler(
    filename = function() {
      template_file <- project_file("evosep_template_file")
      source_name <- if (!is.null(template_file) && !is.null(template_file$name) && nzchar(template_file$name)) {
        template_file$name
      } else {
        "evosep_run_queue.csl"
      }
      ext <- tools::file_ext(source_name)
      stem <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", source_name) else source_name
      out_name <- if (nzchar(ext)) paste0(stem, "_app.", ext) else paste0(stem, "_app.csl")
      build_export_filename(out_name)
    },
    content = function(file) {
      template_file <- project_file("evosep_template_file")
      source_name <- if (!is.null(template_file) && !is.null(template_file$name) && nzchar(template_file$name)) {
        template_file$name
      } else {
        "evosep_run_queue.csl"
      }
      ext <- tools::file_ext(source_name)
      stem <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", source_name) else source_name
      out_name <- if (nzchar(ext)) paste0(stem, "_app.", ext) else paste0(stem, "_app.csl")
      out_name <- build_export_filename(out_name)
      csl_text <- evosep_write_csl_text(evosep_template_text(), evosep_queue_table())
      writeLines(csl_text, file, useBytes = TRUE)
      save_download_copy(file, out_name)
    }
  )

  output$download_project_bundle <- downloadHandler(
    filename = function() {
      priority_ids <- c(
        "meta_file", "run_order_file", "protein_no_impute_file", "protein_imputed_file",
        setdiff(project_file_ids, c("meta_file", "run_order_file", "protein_no_impute_file", "protein_imputed_file"))
      )
      source_names <- vapply(priority_ids, function(id) {
        source <- project_file(id)
        if (is.null(source) || is.null(source$name)) "" else as.character(source$name)
      }, character(1))
      project_bundle_filename(source_names, project_number = project_number_for_downloads())
    },
    content = function(file) {
      priority_ids <- c(
        "meta_file", "run_order_file", "protein_no_impute_file", "protein_imputed_file",
        setdiff(project_file_ids, c("meta_file", "run_order_file", "protein_no_impute_file", "protein_imputed_file"))
      )
      source_names <- vapply(priority_ids, function(id) {
        source <- project_file(id)
        if (is.null(source) || is.null(source$name)) "" else as.character(source$name)
      }, character(1))
      out_name <- project_bundle_filename(source_names, project_number = project_number_for_downloads())
      write_project_bundle(file)
      save_download_copy(file, out_name)
    }
  )

  output$download_metadata_workbook <- downloadHandler(
    filename = function() build_export_filename(input$metadata_workbook_filename),
    content = function(file) {
      validate(need(requireNamespace("openxlsx", quietly = TRUE), "Install the R package 'openxlsx' to download an Excel workbook."))

      workbook <- openxlsx::createWorkbook()
      metadata_export <- exported_metadata()
      header_style <- openxlsx::createStyle(
        fgFill = "#1F4E78",
        fontColour = "#FFFFFF",
        textDecoration = "bold",
        halign = "center",
        border = "Bottom",
        wrapText = TRUE
      )

      add_supplementary_sheet <- function(sheet_name, table_data, with_filter = FALSE) {
        table_data[] <- lapply(table_data, function(column) {
          if (is.numeric(column)) column[is.nan(column) | is.infinite(column)] <- NA_real_
          column
        })
        openxlsx::addWorksheet(workbook, sheet_name)
        openxlsx::writeData(
          workbook,
          sheet = sheet_name,
          x = table_data,
          headerStyle = header_style,
          withFilter = with_filter,
          keepNA = TRUE,
          na.string = "NaN"
        )
        openxlsx::freezePane(workbook, sheet = sheet_name, firstRow = TRUE)
        openxlsx::setRowHeights(workbook, sheet = sheet_name, rows = 1, heights = 42)

        derived_cols <- which(grepl("(_percent_CV|_log2_fold_change|_(paired|unpaired)_t_test_p_value|_BH_FDR)$", colnames(table_data)))
        abundance_cols <- which(grepl("_Protein_group_abundance$", colnames(table_data)))
        digits <- max(1, as.integer(input$numeric_sig_figs))
        threshold <- suppressWarnings(as.numeric(input$scientific_threshold))
        if (!is.finite(threshold) || threshold < 0) threshold <- 0

        abundance_style <- openxlsx::createStyle(numFmt = "0.0E+00")
        if (length(abundance_cols) > 0) {
          openxlsx::addStyle(
            workbook, sheet = sheet_name, style = abundance_style,
            rows = seq_len(nrow(table_data)) + 1, cols = abundance_cols,
            gridExpand = TRUE, stack = TRUE
          )
        }

        make_number_format <- function(value) {
          if (isTRUE(input$scientific_small_values) && value != 0 && abs(value) < threshold) {
            if (digits == 1) return("0E+00")
            return(paste0("0.", paste(rep("0", digits - 1), collapse = ""), "E+00"))
          }
          decimal_places <- if (value == 0) digits - 1 else digits - floor(log10(abs(value))) - 1
          if (decimal_places < 0) {
            if (digits == 1) return("0E+00")
            return(paste0("0.", paste(rep("0", digits - 1), collapse = ""), "E+00"))
          }
          if (decimal_places == 0) return("0")
          paste0("0.", paste(rep("0", decimal_places), collapse = ""))
        }

        for (column_number in derived_cols) {
          values <- suppressWarnings(as.numeric(table_data[[column_number]]))
          valid <- which(is.finite(values))
          if (length(valid) > 0) {
            formats <- vapply(values[valid], make_number_format, character(1))
            for (format_text in unique(formats)) {
              rows <- valid[formats == format_text] + 1
              openxlsx::addStyle(
                workbook, sheet = sheet_name,
                style = openxlsx::createStyle(numFmt = format_text),
                rows = rows, cols = column_number, gridExpand = TRUE, stack = TRUE
              )
            }
          }
        }

        compact_widths <- ifelse(
          grepl("_Protein_group_abundance$|_quantified_precursors$|_percent_CV|_log2_fold_change|_(paired|unpaired)_t_test_p_value|_BH_FDR", colnames(table_data)),
          12,
          pmin(pmax(nchar(colnames(table_data)) + 2, 8), 22)
        )
        openxlsx::setColWidths(workbook, sheet = sheet_name, cols = seq_len(ncol(table_data)), widths = compact_widths)
      }

      add_supplementary_sheet("Table S1. Metadata", metadata_export, with_filter = TRUE)
      if (protein_source_available("S2")) {
        add_supplementary_sheet("Table S2. Protein, no impute", protein_no_impute_table(), with_filter = isTRUE(input$enable_protein_filters))
      }
      if (protein_source_available("S3")) {
        add_supplementary_sheet("Table S3. Protein, imputed", protein_imputed_table(), with_filter = isTRUE(input$enable_protein_filters))
      }
      openxlsx::saveWorkbook(workbook, file, overwrite = TRUE)
      save_download_copy(file, build_export_filename(input$metadata_workbook_filename))
    }
  )

  output$download_png <- downloadHandler(
    filename = function() build_export_filename(input$png_filename),
    content = function(file) {
      out_name <- build_export_filename(input$png_filename)
      ggsave(file, plot = pca_plot_obj(), width = input$pca_figure_width, height = input$pca_figure_height, dpi = 600)
      save_download_copy(file, out_name)
    }
  )

  output$download_svg <- downloadHandler(
    filename = function() build_export_filename(input$svg_filename),
    content = function(file) {
      out_name <- build_export_filename(input$svg_filename)
      ggsave(file, plot = pca_plot_obj(), width = input$pca_figure_width, height = input$pca_figure_height, device = svglite::svglite)
      save_download_copy(file, out_name)
    }
  )

  output$download_clustvis_pca_png <- downloadHandler(
    filename = function() build_export_filename("clustvis_like_pca.png"),
    content = function(file) {
      out_name <- build_export_filename("clustvis_like_pca.png")
      ggsave(file, plot = clustvis_pca_plot_obj(), width = input$clustvis_pca_width, height = input$clustvis_pca_height, dpi = 600)
      save_download_copy(file, out_name)
    }
  )

  output$download_clustvis_pca_svg <- downloadHandler(
    filename = function() build_export_filename("clustvis_like_pca.svg"),
    content = function(file) {
      validate(need(requireNamespace("svglite", quietly = TRUE), "Install the R package 'svglite' to download SVG files."))
      out_name <- build_export_filename("clustvis_like_pca.svg")
      ggsave(file, plot = clustvis_pca_plot_obj(), width = input$clustvis_pca_width, height = input$clustvis_pca_height, device = svglite::svglite)
      save_download_copy(file, out_name)
    }
  )

  output$download_pca_loadings_csv <- downloadHandler(
    filename = function() build_export_filename("pca_pc1_pc2_protein_loadings.csv"),
    content = function(file) {
      out_name <- build_export_filename("pca_pc1_pc2_protein_loadings.csv")
      write.csv(pca_loadings_data(), file, row.names = FALSE, na = "NaN")
      save_download_copy(file, out_name)
    }
  )

  output$download_gene_set_enrichment_csv <- downloadHandler(
    filename = function() {
      comparison <- if (is.null(input$gsea_comparison)) "comparison" else stats_comparison_prefix(input$gsea_comparison)
      build_export_filename(paste0("gene_set_enrichment_", comparison, ".csv"))
    },
    content = function(file) {
      comparison <- if (is.null(input$gsea_comparison)) "comparison" else stats_comparison_prefix(input$gsea_comparison)
      out_name <- build_export_filename(paste0("gene_set_enrichment_", comparison, ".csv"))
      write.csv(gene_set_enrichment_data(), file, row.names = FALSE, na = "NaN")
      save_download_copy(file, out_name)
    }
  )

  output$download_gene_set_enrichment_png <- downloadHandler(
    filename = function() {
      comparison <- if (is.null(input$gsea_comparison)) "comparison" else stats_comparison_prefix(input$gsea_comparison)
      build_export_filename(paste0("gene_set_enrichment_", comparison, ".png"))
    },
    content = function(file) {
      comparison <- if (is.null(input$gsea_comparison)) "comparison" else stats_comparison_prefix(input$gsea_comparison)
      out_name <- build_export_filename(paste0("gene_set_enrichment_", comparison, ".png"))
      ggsave(file, plot = gene_set_enrichment_plot_obj(), width = input$gsea_plot_width, height = input$gsea_plot_height, dpi = 300)
      save_download_copy(file, out_name)
    }
  )

  output$download_gene_set_enrichment_svg <- downloadHandler(
    filename = function() {
      comparison <- if (is.null(input$gsea_comparison)) "comparison" else stats_comparison_prefix(input$gsea_comparison)
      build_export_filename(paste0("gene_set_enrichment_", comparison, ".svg"))
    },
    content = function(file) {
      validate(need(requireNamespace("svglite", quietly = TRUE), "Install the R package 'svglite' to download SVG files."))
      comparison <- if (is.null(input$gsea_comparison)) "comparison" else stats_comparison_prefix(input$gsea_comparison)
      out_name <- build_export_filename(paste0("gene_set_enrichment_", comparison, ".svg"))
      ggsave(file, plot = gene_set_enrichment_plot_obj(), width = input$gsea_plot_width, height = input$gsea_plot_height, device = svglite::svglite)
      save_download_copy(file, out_name)
    }
  )

  output$download_feature_png <- downloadHandler(
    filename = function() {
      features <- as.character(input$feature_select)
      stem <- if (length(features) == 1L) gsub("[^A-Za-z0-9_\\-]", "_", features) else paste0(length(features), "_proteins")
      build_export_filename(paste0(stem, "_barplot.png"))
    },
    content = function(file) {
      features <- as.character(input$feature_select)
      stem <- if (length(features) == 1L) gsub("[^A-Za-z0-9_\\-]", "_", features) else paste0(length(features), "_proteins")
      out_name <- build_export_filename(paste0(stem, "_barplot.png"))
      ggsave(file, plot = feature_plot_obj(), width = input$feature_figure_width, height = input$feature_figure_height, dpi = 600)
      save_download_copy(file, out_name)
    }
  )

  output$download_feature_svg <- downloadHandler(
    filename = function() {
      features <- as.character(input$feature_select)
      stem <- if (length(features) == 1L) gsub("[^A-Za-z0-9_\\-]", "_", features) else paste0(length(features), "_proteins")
      build_export_filename(paste0(stem, "_barplot.svg"))
    },
    content = function(file) {
      validate(need(requireNamespace("svglite", quietly = TRUE), "Install the R package 'svglite' to download SVG files."))
      features <- as.character(input$feature_select)
      stem <- if (length(features) == 1L) gsub("[^A-Za-z0-9_\\-]", "_", features) else paste0(length(features), "_proteins")
      out_name <- build_export_filename(paste0(stem, "_barplot.svg"))
      ggsave(file, plot = feature_plot_obj(), width = input$feature_figure_width, height = input$feature_figure_height, device = svglite::svglite)
      save_download_copy(file, out_name)
    }
  )

  output$download_script_box_png <- downloadHandler(
    filename = function() {
      features <- as.character(input$script_box_features)
      stem <- if (length(features) == 1L) gsub("[^A-Za-z0-9_\\-]", "_", features) else paste0(length(features), "_protein_groups")
      build_export_filename(paste0(stem, "_script_boxplot.png"))
    },
    content = function(file) {
      features <- as.character(input$script_box_features)
      stem <- if (length(features) == 1L) gsub("[^A-Za-z0-9_\\-]", "_", features) else paste0(length(features), "_protein_groups")
      out_name <- build_export_filename(paste0(stem, "_script_boxplot.png"))
      ggsave(file, plot = script_box_plot_obj(), width = input$script_box_width, height = input$script_box_height, dpi = 600)
      save_download_copy(file, out_name)
    }
  )

  output$download_script_box_svg <- downloadHandler(
    filename = function() {
      features <- as.character(input$script_box_features)
      stem <- if (length(features) == 1L) gsub("[^A-Za-z0-9_\\-]", "_", features) else paste0(length(features), "_protein_groups")
      build_export_filename(paste0(stem, "_script_boxplot.svg"))
    },
    content = function(file) {
      validate(need(requireNamespace("svglite", quietly = TRUE), "Install the R package 'svglite' to download SVG files."))
      features <- as.character(input$script_box_features)
      stem <- if (length(features) == 1L) gsub("[^A-Za-z0-9_\\-]", "_", features) else paste0(length(features), "_protein_groups")
      out_name <- build_export_filename(paste0(stem, "_script_boxplot.svg"))
      ggsave(file, plot = script_box_plot_obj(), width = input$script_box_width, height = input$script_box_height, device = svglite::svglite)
      save_download_copy(file, out_name)
    }
  )

  output$download_correlation_png <- downloadHandler(
    filename = function() {
      feature <- if (is.null(input$correlation_feature)) "reference" else input$correlation_feature
      build_export_filename(paste0(gsub("[^A-Za-z0-9_\\-]", "_", feature), "_correlation_lollipop.png"))
    },
    content = function(file) {
      feature <- if (is.null(input$correlation_feature)) "reference" else input$correlation_feature
      out_name <- build_export_filename(paste0(gsub("[^A-Za-z0-9_\\-]", "_", feature), "_correlation_lollipop.png"))
      ggsave(file, plot = correlation_lollipop_plot_obj(), width = input$correlation_figure_width, height = input$correlation_figure_height, dpi = 600)
      save_download_copy(file, out_name)
    }
  )

  output$download_correlation_svg <- downloadHandler(
    filename = function() {
      feature <- if (is.null(input$correlation_feature)) "reference" else input$correlation_feature
      build_export_filename(paste0(gsub("[^A-Za-z0-9_\\-]", "_", feature), "_correlation_lollipop.svg"))
    },
    content = function(file) {
      validate(need(requireNamespace("svglite", quietly = TRUE), "Install the R package 'svglite' to download SVG files."))
      feature <- if (is.null(input$correlation_feature)) "reference" else input$correlation_feature
      out_name <- build_export_filename(paste0(gsub("[^A-Za-z0-9_\\-]", "_", feature), "_correlation_lollipop.svg"))
      ggsave(file, plot = correlation_lollipop_plot_obj(), width = input$correlation_figure_width, height = input$correlation_figure_height, device = svglite::svglite)
      save_download_copy(file, out_name)
    }
  )

  output$download_correlation_csv <- downloadHandler(
    filename = function() {
      feature <- if (is.null(input$correlation_feature)) "reference" else input$correlation_feature
      build_export_filename(paste0(gsub("[^A-Za-z0-9_\\-]", "_", feature), "_correlation_results.csv"))
    },
    content = function(file) {
      feature <- if (is.null(input$correlation_feature)) "reference" else input$correlation_feature
      out_name <- build_export_filename(paste0(gsub("[^A-Za-z0-9_\\-]", "_", feature), "_correlation_results.csv"))
      write.csv(correlation_results_data(), file, row.names = FALSE, na = "NaN")
      save_download_copy(file, out_name)
    }
  )

  output$download_cv_png <- downloadHandler(
    filename = function() build_export_filename("protein_group_cv_distribution_per_condition.png"),
    content = function(file) {
      out_name <- build_export_filename("protein_group_cv_distribution_per_condition.png")
      ggsave(file, plot = cv_plot_obj(), width = input$cv_figure_width, height = input$cv_figure_height, dpi = 300)
      save_download_copy(file, out_name)
    }
  )

  output$download_cv_svg <- downloadHandler(
    filename = function() build_export_filename("protein_group_cv_distribution_per_condition.svg"),
    content = function(file) {
      validate(need(requireNamespace("svglite", quietly = TRUE), "Install the R package 'svglite' to download SVG files."))
      out_name <- build_export_filename("protein_group_cv_distribution_per_condition.svg")
      ggsave(file, plot = cv_plot_obj(), width = input$cv_figure_width, height = input$cv_figure_height, device = svglite::svglite)
      save_download_copy(file, out_name)
    }
  )

  output$download_identification_overview_png <- downloadHandler(
    filename = function() {
      metric <- if (input$identification_metric == "ProteinGroups") "protein_groups" else "precursors"
      build_export_filename(paste0("run_identifications_overview_", metric, ".png"))
    },
    content = function(file) {
      metric <- if (input$identification_metric == "ProteinGroups") "protein_groups" else "precursors"
      out_name <- build_export_filename(paste0("run_identifications_overview_", metric, ".png"))
      ggsave(file, plot = identification_overview_plot_obj(), width = input$identification_overview_width, height = input$identification_overview_height, dpi = 300)
      save_download_copy(file, out_name)
    }
  )

  output$download_identification_overview_svg <- downloadHandler(
    filename = function() {
      metric <- if (input$identification_metric == "ProteinGroups") "protein_groups" else "precursors"
      build_export_filename(paste0("run_identifications_overview_", metric, ".svg"))
    },
    content = function(file) {
      validate(need(requireNamespace("svglite", quietly = TRUE), "Install the R package 'svglite' to download SVG files."))
      metric <- if (input$identification_metric == "ProteinGroups") "protein_groups" else "precursors"
      out_name <- build_export_filename(paste0("run_identifications_overview_", metric, ".svg"))
      ggsave(file, plot = identification_overview_plot_obj(), width = input$identification_overview_width, height = input$identification_overview_height, device = svglite::svglite)
      save_download_copy(file, out_name)
    }
  )

  output$download_run_identifications_png <- downloadHandler(
    filename = function() {
      metric <- if (input$identification_metric == "ProteinGroups") "protein_groups" else "precursors"
      build_export_filename(paste0("run_identifications_stacked_", metric, ".png"))
    },
    content = function(file) {
      metric <- if (input$identification_metric == "ProteinGroups") "protein_groups" else "precursors"
      out_name <- build_export_filename(paste0("run_identifications_stacked_", metric, ".png"))
      ggsave(file, plot = run_identifications_plot_obj(), width = input$run_identifications_width, height = input$run_identifications_height, dpi = 300)
      save_download_copy(file, out_name)
    }
  )

  output$download_run_identifications_svg <- downloadHandler(
    filename = function() {
      metric <- if (input$identification_metric == "ProteinGroups") "protein_groups" else "precursors"
      build_export_filename(paste0("run_identifications_stacked_", metric, ".svg"))
    },
    content = function(file) {
      validate(need(requireNamespace("svglite", quietly = TRUE), "Install the R package 'svglite' to download SVG files."))
      metric <- if (input$identification_metric == "ProteinGroups") "protein_groups" else "precursors"
      out_name <- build_export_filename(paste0("run_identifications_stacked_", metric, ".svg"))
      ggsave(file, plot = run_identifications_plot_obj(), width = input$run_identifications_width, height = input$run_identifications_height, device = svglite::svglite)
      save_download_copy(file, out_name)
    }
  )

  output$download_volcano_png <- downloadHandler(
    filename = function() {
      comparison <- if (is.null(input$volcano_comparison)) "comparison" else stats_comparison_prefix(input$volcano_comparison)
      build_export_filename(paste0("volcano_", comparison, ".png"))
    },
    content = function(file) {
      comparison <- if (is.null(input$volcano_comparison)) "comparison" else stats_comparison_prefix(input$volcano_comparison)
      out_name <- build_export_filename(paste0("volcano_", comparison, ".png"))
      ggsave(file, plot = volcano_plot_obj(), width = input$volcano_figure_width, height = input$volcano_figure_height, dpi = 300)
      save_download_copy(file, out_name)
    }
  )

  output$download_volcano_svg <- downloadHandler(
    filename = function() {
      comparison <- if (is.null(input$volcano_comparison)) "comparison" else stats_comparison_prefix(input$volcano_comparison)
      build_export_filename(paste0("volcano_", comparison, ".svg"))
    },
    content = function(file) {
      validate(need(requireNamespace("svglite", quietly = TRUE), "Install the R package 'svglite' to download SVG files."))
      comparison <- if (is.null(input$volcano_comparison)) "comparison" else stats_comparison_prefix(input$volcano_comparison)
      out_name <- build_export_filename(paste0("volcano_", comparison, ".svg"))
      ggsave(file, plot = volcano_plot_obj(), width = input$volcano_figure_width, height = input$volcano_figure_height, device = svglite::svglite)
      save_download_copy(file, out_name)
    }
  )

  output$download_volcano_html <- downloadHandler(
    filename = function() {
      comparison <- if (is.null(input$volcano_comparison)) "comparison" else stats_comparison_prefix(input$volcano_comparison)
      build_export_filename(paste0("volcano_interactive_", comparison, ".zip"))
    },
    content = function(file) {
      validate(need(requireNamespace("htmlwidgets", quietly = TRUE) && requireNamespace("zip", quietly = TRUE), "Install R packages 'htmlwidgets' and 'zip' to export the interactive volcano plot."))
      interactive_dir <- tempfile("volcano_interactive_")
      dir.create(interactive_dir)
      html_file <- file.path(interactive_dir, "volcano_interactive.html")
      htmlwidgets::saveWidget(volcano_interactive_obj(), file = html_file, selfcontained = FALSE)
      zip::zipr(file, list.files(interactive_dir, recursive = TRUE, full.names = TRUE), root = interactive_dir)
      comparison <- if (is.null(input$volcano_comparison)) "comparison" else stats_comparison_prefix(input$volcano_comparison)
      save_download_copy(file, build_export_filename(paste0("volcano_interactive_", comparison, ".zip")))
    }
  )

  output$download_volcano_counts_csv <- downloadHandler(
    filename = function() build_export_filename("significant_protein_counts_by_comparison.csv"),
    content = function(file) {
      data <- volcano_counts_data()
      write.csv(data[, setdiff(colnames(data), "ComparisonID"), drop = FALSE], file, row.names = FALSE, na = "NaN")
      save_download_copy(file, build_export_filename("significant_protein_counts_by_comparison.csv"))
    }
  )

  output$download_all <- downloadHandler(
    filename = function() build_export_filename("pca_exports.zip"),
    content = function(file) {
      tmpdir <- tempfile("pca_export_dir_")
      dir.create(tmpdir)
      scores_file <- file.path(tmpdir, build_export_filename("pca_scores.csv"))
      png_file <- file.path(tmpdir, build_export_filename(input$png_filename))
      svg_file <- file.path(tmpdir, build_export_filename(input$svg_filename))

      write.csv(pca_results()$scores, scores_file, row.names = FALSE)
      ggsave(png_file, plot = pca_plot_obj(), width = input$pca_figure_width, height = input$pca_figure_height, dpi = 600)
      ggsave(svg_file, plot = pca_plot_obj(), width = input$pca_figure_width, height = input$pca_figure_height, device = svglite::svglite)

      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(tmpdir)
      utils::zip(zipfile = file, files = basename(c(scores_file, png_file, svg_file)))
      save_download_copy(file, build_export_filename("pca_exports.zip"))
    }
  )
}

shinyApp(ui, server)
