app_text <- paste(readLines(if (file.exists("app.R")) "app.R" else "../app.R", warn = FALSE), collapse = "\n")

function_start <- regexpr("protein_sample_names_for_table <- function", app_text, fixed = TRUE)[1L]
stopifnot(function_start > 0L)
function_text <- substr(app_text, function_start, function_start + 1200L)
stopifnot(grepl("resolve_proteomics_processed_sample_ids", function_text, fixed = TRUE))
stopifnot(grepl("processed_sample_map", function_text, fixed = TRUE))

cat("PCA restored-header mapping static test passed.\n")
