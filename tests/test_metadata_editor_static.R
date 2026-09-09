text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")

for (pattern in c(
  'h4("Editable metadata")',
  'actionButton("discard_metadata_edits", "Discard all metadata edits")',
  'observeEvent(input$metadata_preview_cell_edit',
  'locked_metadata_columns <- c("Sample", "RunOrder", "Run Order", "Filename", "FileName", "File Name")',
  'selected <- input$metadata_export_columns',
  'view <- prepare_metadata_editor_view(table, selected, locked_metadata_columns)',
  'isolate(spqc_metadata_draft_edits())',
  'stateSave = TRUE'
)) stopifnot(grepl(pattern, text, fixed = TRUE))

for (pattern in c(
  'h4("Detected SPQC metadata rows")',
  'DTOutput("spqc_metadata_preview")',
  'observeEvent(input$spqc_metadata_preview_cell_edit',
  'actionButton("clear_spqc_metadata_edits"'
)) stopifnot(!grepl(pattern, text, fixed = TRUE))

cat("Unified metadata editor locks identifiers and avoids redraws on cell edits.\n")
