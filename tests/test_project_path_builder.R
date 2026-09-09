source(if (file.exists("project_duckdb.R")) "project_duckdb.R" else "../project_duckdb.R")

folder <- file.path(tempdir(), "proteomics-projects")

path <- proteomics_project_path(folder, "BM project")
stopifnot(identical(path, normalizePath(file.path(folder, "BM project.duckdb"), mustWork = FALSE)))

path_with_extension <- proteomics_project_path(folder, "BM_project.duckdb")
stopifnot(identical(path_with_extension, normalizePath(file.path(folder, "BM_project.duckdb"), mustWork = FALSE)))

tryCatch(
  {
    proteomics_project_path("", "BM_project")
    stop("Expected a missing-folder error")
  },
  error = function(e) stopifnot(grepl("folder", conditionMessage(e), ignore.case = TRUE))
)

tryCatch(
  {
    proteomics_project_path(folder, "")
    stop("Expected a missing-filename error")
  },
  error = function(e) stopifnot(grepl("filename", conditionMessage(e), ignore.case = TRUE))
)

message("Project folder and filename path tests passed.")
