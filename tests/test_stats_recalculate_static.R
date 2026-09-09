text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")

stopifnot(grepl('actionButton("run_protein_stats", "Recalculate statistics"', text, fixed = TRUE))
stopifnot(grepl('verbatimTextOutput("protein_stats_status")', text, fixed = TRUE))
stopifnot(grepl('protein_stats_refresh_revision(0L)', text, fixed = TRUE))
stopifnot(grepl('protein_stats_refresh_revision <- reactiveVal(0L)', text, fixed = TRUE))
stopifnot(grepl('observeEvent(input$run_protein_stats', text, fixed = TRUE))
stopifnot(grepl('protein_stats_refresh_revision()', text, fixed = TRUE))
stopifnot(grepl('actionButton("stop_protein_stats", "Stop statistics run")', text, fixed = TRUE))

cat("Explicit statistics recalculation controls OK\n")
