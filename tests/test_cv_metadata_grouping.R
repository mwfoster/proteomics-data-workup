source(if (file.exists("selection_helpers.R")) "selection_helpers.R" else "../selection_helpers.R")

metadata <- data.frame(
  Sample = c("S1", "S2", "S3", "S4"),
  Condition = c("Case", "Case", "Control", "Control"),
  Batch = c("1", "1", "2", "2"),
  Constant = "same",
  stringsAsFactors = FALSE
)

stopifnot(identical(cv_metadata_group_columns(metadata), c("Condition", "Batch")))
stopifnot(identical(
  cv_sample_groups(metadata, "Batch", c("S4", "S1", "missing")),
  c(S4 = "2", S1 = "1", missing = NA_character_)
))

app_text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")
stopifnot(grepl('selectInput("cv_plot_group_col", "Group CV by metadata field"', app_text, fixed = TRUE))
stopifnot(grepl('"cv_plot_group_col"', app_text, fixed = TRUE))

message("CV metadata grouping tests passed.")
