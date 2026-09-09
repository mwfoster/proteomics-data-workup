text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")
required_patterns <- c(
  'settings$workflow_tabs <- input$workflow_tabs',
  'updateTabsetPanel(session, "workflow_tabs"',
  'input$metadata_preview_cell_edit',
  'spqc_metadata_draft_edits(draft_edits)',
  'spqc_metadata_edits(candidate$spqc_metadata_edits)',
  'input$apply_metadata_changes',
  'input$discard_metadata_edits',
  '"cv_conditions"', '"stats_comparisons"', '"stats_paired_comparisons"',
  '"protein_header_label_columns"', '"correlation_groups"'
)
required_patterns <- c(
  required_patterns,
  '"feature_group_by"', '"feature_label_by"',
  '"script_box_features"', '"script_box_label_by"'
)
for (pattern in required_patterns) stopifnot(grepl(pattern, text, fixed = TRUE))
message("Project settings static checks passed.")
