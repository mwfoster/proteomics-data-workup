empty_proteomics_metadata_state <- function() {
  list(metadata = NULL, columns = character(0), spqc_metadata_edits = empty_spqc_metadata_edits())
}

merge_proteomics_metadata_edits <- function(committed, draft) {
  committed <- if (is.null(committed)) empty_spqc_metadata_edits() else committed
  draft <- if (is.null(draft)) empty_spqc_metadata_edits() else draft
  if (!nrow(draft)) return(committed)
  draft_keys <- paste(as.character(draft$Sample), as.character(draft$Column), sep = "\r")
  if (nrow(committed)) {
    committed_keys <- paste(as.character(committed$Sample), as.character(committed$Column), sep = "\r")
    committed <- committed[!committed_keys %in% draft_keys, , drop = FALSE]
  }
  rbind(committed, draft)
}

build_proteomics_applied_state <- function(metadata, columns, committed_edits, draft_edits) {
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(metadata) || !"Sample" %in% colnames(metadata)) stop("Metadata must contain at least one Sample row.")
  samples <- as.character(metadata$Sample)
  if (any(is.na(samples) | !nzchar(trimws(samples)))) stop("Metadata Sample values must not be blank.")
  if (anyDuplicated(samples)) stop("Metadata Sample values must be unique before applying changes.")
  columns <- unique(as.character(columns))
  columns <- columns[!is.na(columns) & nzchar(columns)]
  if (!length(columns) || any(!columns %in% colnames(metadata))) stop("Every selected Table S1 column must exist in metadata.")
  edits <- merge_proteomics_metadata_edits(committed_edits, draft_edits)
  list(
    metadata = apply_proteomics_metadata_cell_edits(metadata, edits),
    columns = columns,
    spqc_metadata_edits = edits
  )
}
