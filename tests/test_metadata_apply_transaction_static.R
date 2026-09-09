text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")

required <- c(
  'actionButton("apply_metadata_changes", "Apply metadata changes", class = "btn-primary")',
  'actionButton("apply_metadata_replacement", "Load modified metadata into draft")',
  'candidate <- build_proteomics_applied_state(',
  'applied_metadata_state(candidate)',
  'metadata_apply_revision(isolate(metadata_apply_revision()) + 1L)',
  'autosave_active_project("applied metadata changes", include_derived = FALSE)'
)
for (pattern in required) stopifnot(grepl(pattern, text, fixed = TRUE))
stopifnot(!grepl('actionButton("apply_spqc_metadata_edits"', text, fixed = TRUE))

observer_start <- regexpr("observeEvent(input$apply_metadata_changes", text, fixed = TRUE)[1L]
stopifnot(observer_start > 0L)
observer <- substr(text, observer_start, observer_start + 3500L)
candidate_position <- regexpr("candidate <- build_proteomics_applied_state(", observer, fixed = TRUE)[1L]
commit_position <- regexpr("applied_metadata_state(candidate)", observer, fixed = TRUE)[1L]
stopifnot(candidate_position > 0L, commit_position > candidate_position)

cat("Metadata apply transaction static checks passed.\n")
