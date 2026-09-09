text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")

stopifnot(grepl(
  "clustvis_pca_results <- eventReactive(list(input$run_clustvis_pca, metadata_apply_revision()), {",
  text,
  fixed = TRUE
))
stopifnot(grepl("req(input$run_clustvis_pca > 0)", text, fixed = TRUE))

cat("Applied metadata refreshes an already-run PCA.\n")
