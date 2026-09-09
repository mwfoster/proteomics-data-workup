source(if (file.exists("selection_helpers.R")) "selection_helpers.R" else "../selection_helpers.R")

metadata <- data.frame(
  Sample = paste0("S", 1:6),
  Cohort = c("A", "A", "A", "B", "B", "B"),
  Sex = c("F", "M", "F", "M", "F", "M"),
  Constant = rep("same", 6),
  stringsAsFactors = FALSE
)

stopifnot(identical(stats_metadata_group_columns(metadata), c("Cohort", "Sex")))

app_text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")
stopifnot(!grepl('if (!"Condition" %in% colnames(md)) return()', app_text, fixed = TRUE))
generic_restore_block <- sub(
  '.*update_selectize <- c\\((.*?)\\)\\n    for \\(id in update_text\\).*',
  '\\1',
  app_text
)
stopifnot(!grepl('"protein_quantity_order_columns"', generic_restore_block, fixed = TRUE))
stopifnot(!grepl('"stats_group_columns"', generic_restore_block, fixed = TRUE))
stopifnot(grepl('observeEvent(input$workflow_tabs', app_text, fixed = TRUE))
stopifnot(grepl('identical(input$workflow_tabs, "Protein tables")', app_text, fixed = TRUE))
stopifnot(grepl('session$onFlushed(function() {\n      isolate(sync_restored_project_inputs(settings))', app_text, fixed = TRUE))
stopifnot(grepl('isolate(update_protein_header_label_columns(downstream_metadata_columns(), settings$protein_header_label_columns))', app_text, fixed = TRUE))
stopifnot(grepl('isolate(sync_condition_dependent_inputs(settings))', app_text, fixed = TRUE))

cat("Statistics metadata column discovery works without a Condition column\n")
