app_dir <- normalizePath(dirname(sys.frame(1)$ofile %||% getwd()), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(app_dir, "app.R"))) {
  app_dir <- getwd()
}

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Package 'shiny' is required. Run source('install_packages.R') first.", call. = FALSE)
}

options(shiny.maxRequestSize = 1024 * 1024 * 1024)
shiny::runApp(app_dir, host = "0.0.0.0", port = 6875, launch.browser = interactive())
