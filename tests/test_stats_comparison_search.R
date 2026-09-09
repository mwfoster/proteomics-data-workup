source(if (file.exists("selection_helpers.R")) "selection_helpers.R" else "../selection_helpers.R")

options <- stats_comparison_selectize_options()

stopifnot(identical(unlist(options$plugins), c("drag_drop", "remove_button")))
stopifnot(identical(options$searchField, c("text", "value")))
stopifnot(identical(options$closeAfterSelect, FALSE))
stopifnot(is.numeric(options$maxOptions) && options$maxOptions >= 1000)
stopifnot(grepl("Type", options$placeholder, fixed = TRUE))

cat("Statistics comparison selector supports type-to-search for large choice lists.\n")
