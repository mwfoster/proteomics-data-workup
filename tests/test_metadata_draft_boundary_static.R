text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")

stopifnot(grepl("applied_metadata_state <- reactiveVal(empty_proteomics_metadata_state())", text, fixed = TRUE))
stopifnot(grepl("draft_metadata <- reactive({", text, fixed = TRUE))
stopifnot(grepl("built_metadata <- reactive({\n    metadata_apply_revision()\n    state <- applied_metadata_state()", text, fixed = TRUE))
stopifnot(grepl('validate(need(!is.null(state$metadata), "Apply metadata changes on the Make metadata tab first."))', text, fixed = TRUE))
stopifnot(grepl("draft_exported_metadata <- reactive({", text, fixed = TRUE))
stopifnot(grepl("built <- draft_metadata()", text, fixed = TRUE))
stopifnot(grepl("table <- draft_metadata()", text, fixed = TRUE))
stopifnot(grepl("table <- apply_proteomics_metadata_cell_edits(table, isolate(spqc_metadata_draft_edits()))", text, fixed = TRUE))
stopifnot(grepl("downstream_metadata_columns <- reactive({\n    state <- applied_metadata_state()", text, fixed = TRUE))
stopifnot(!grepl("table_s1_metadata_choices(available, input$metadata_export_columns)", text, fixed = TRUE))

cat("Metadata draft boundary static checks passed.\n")
