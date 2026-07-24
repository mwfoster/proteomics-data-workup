r_version_lib <- paste(R.version$major, sub("\\..*", "", R.version$minor), sep = ".")
user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (!nzchar(user_lib) || normalizePath(user_lib, mustWork = FALSE) == normalizePath(file.path(path.expand("~"), "R", "library"), mustWork = FALSE)) {
  user_lib <- file.path(path.expand("~"), "R", "library", r_version_lib)
}
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_LIBS_USER = user_lib, R_LIBS_SITE = "", R_LIBS = user_lib)
try(assign(".Library.site", character(0), envir = baseenv()), silent = TRUE)
.libPaths(unique(c(user_lib, file.path(R.home(), "library"))))

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
