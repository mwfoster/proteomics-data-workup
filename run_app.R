args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
app_file <- sub(file_arg, "", args[grepl(paste0("^", file_arg), args)])
app_dir <- if (length(app_file) > 0) dirname(normalizePath(app_file[1], winslash = "/", mustWork = FALSE)) else getwd()
if (!file.exists(file.path(app_dir, "app.R"))) app_dir <- getwd()

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Package 'shiny' is required. Run source('install_packages.R') first.", call. = FALSE)
}

options(shiny.maxRequestSize = 1024 * 1024 * 1024)
shiny::runApp(app_dir, host = "0.0.0.0", port = 6875, launch.browser = interactive())
