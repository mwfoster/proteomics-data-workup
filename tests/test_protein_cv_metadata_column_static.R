text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")

for (pattern in c(
  'selectInput("protein_cv_group_col", "Metadata column for %CV calculation"',
  '"Groups for %CV calculation"',
  'cv_group_col <- input$protein_cv_group_col',
  'md[[cv_group_col]]',
  'updateSelectInput(session, "protein_cv_group_col"'
)) stopifnot(grepl(pattern, text, fixed = TRUE))

cat("Protein-table CV calculation supports a selected metadata column.\n")
