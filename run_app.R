script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

if (is.null(script_file) || !nzchar(script_file)) {
  script_file <- "run_app.R"
}

app_dir <- dirname(normalizePath(script_file, mustWork = FALSE))

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Package 'shiny' is required. Run source('install_packages.R') first.", call. = FALSE)
}

shiny::runApp(app_dir, launch.browser = TRUE)
