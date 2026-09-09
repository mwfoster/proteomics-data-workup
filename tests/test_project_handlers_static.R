text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")
required <- c(
  "load_project_into_session <- function", "save_active_project_state <- function",
  "autosave_active_project <- function", "shiny::debounce", "3000",
  "project_cache.duckdb", "tempfile(\"proteomics_project_\", fileext = extension)",
  "active_project_display_name", "observeEvent(input$project_open_file",
  "showNotification(conditionMessage(e), type = \"error\")"
)
for (pattern in required) stopifnot(grepl(pattern, text, fixed = TRUE))
message("Project handler static checks passed.")
