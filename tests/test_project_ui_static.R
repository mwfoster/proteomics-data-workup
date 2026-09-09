text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")
required <- c("project_open_file", "create_new_project", "project_path_preview", "confirm_create_new_project", "new_project_filename", "clear_active_project", "confirm_clear_active_project", "project_status", "download_active_metadata", "metadata_replacement_file", "apply_metadata_replacement", "script_box_point_opacity")
for (id in required) stopifnot(grepl(paste0('"', id, '"'), text, fixed = TRUE))
stopifnot(grepl('fileInput("project_open_file", "Open existing project"', text, fixed = TRUE))
for (pattern in c('downloadButton("download_project_bundle"', 'checkboxInput("project_include_readme"', 'checkboxInput("project_include_session_info"', 'utils::choose.files(', 'utils::choose.dir(', 'project_folder_browser_entries', 'project_file_browser_entries')) {
  stopifnot(!grepl(pattern, text, fixed = TRUE))
}
stopifnot(!grepl('actionButton("save_active_project"', text, fixed = TRUE))
stopifnot(!grepl('checkboxInput("autosave_project"', text, fixed = TRUE))
stopifnot(grepl('downloadButton("download_project_duckdb", "Save As / Download project"', text, fixed = TRUE))
stopifnot(grepl('observeEvent(input$project_open_file', text, fixed = TRUE))
stopifnot(grepl('observeEvent(input$create_new_project', text, fixed = TRUE))
stopifnot(grepl('output$project_path_preview <- renderText', text, fixed = TRUE))
stopifnot(grepl('observeEvent(input$clear_active_project', text, fixed = TRUE))
stopifnot(grepl('observeEvent(input$confirm_clear_active_project', text, fixed = TRUE))
stopifnot(grepl('session$reload()', text, fixed = TRUE))
stopifnot(grepl('No file on your computer will be deleted or modified.', text, fixed = TRUE))
stopifnot(grepl('loaded_proteomics_project_status(project_file_status())', text, fixed = TRUE))
stopifnot(!grepl('applied_metadata_state()$metadata', text, fixed = TRUE))
stopifnot(grepl('isolate(draft_metadata())', text, fixed = TRUE))
stopifnot(grepl('isolate(spqc_metadata_edits())', text, fixed = TRUE))
stopifnot(grepl('geom_point(', text, fixed = TRUE) && grepl('alpha = input$script_box_point_opacity', text, fixed = TRUE))

pca_heading <- regexpr('h4("ClustVis-like PCA")', text, fixed = TRUE)[1]
pca_button <- regexpr('actionButton("run_clustvis_pca", "Run PCA")', text, fixed = TRUE)[1]
pca_source <- regexpr('radioButtons(\n                "clustvis_pca_source"', text, fixed = TRUE)[1]
stopifnot(pca_heading > 0, pca_button > pca_heading, pca_button < pca_source)
forbidden <- c("RSQLite::SQLite", "project_db_name", "Prototype: save active projects into a local SQLite database")
for (pattern in forbidden) stopifnot(!grepl(pattern, text, fixed = TRUE))
message("Project UI static checks passed.")
