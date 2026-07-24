user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (!nzchar(user_lib)) {
  user_lib <- file.path(path.expand("~"), "R", "library")
  Sys.setenv(R_LIBS_USER = user_lib)
}
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(user_lib, .libPaths())))

if (.Platform$OS.type != "windows") {
  r_dir <- file.path(path.expand("~"), ".R")
  makevars <- file.path(r_dir, "Makevars")
  dir.create(r_dir, recursive = TRUE, showWarnings = FALSE)
  existing <- if (file.exists(makevars)) readLines(makevars, warn = FALSE) else character()
  needed <- c(
    "",
    "# Proteomics Data Workup: required by BH/Boost-backed packages such as fgsea.",
    "CXX11STD = -std=gnu++14",
    "CXX14STD = -std=gnu++14"
  )
  if (!any(grepl("^CXX11STD\\s*=\\s*-std=gnu\\+\\+14", existing))) {
    writeLines(c(existing, needed), makevars)
  }
}

install_if_missing <- function(packages, repos = "https://cloud.r-project.org") {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing, repos = repos, dependencies = TRUE)
  }
}

cran_packages <- c(
  "shiny",
  "ggplot2",
  "DT",
  "dplyr",
  "stringr",
  "missMDA",
  "FactoMineR",
  "svglite",
  "readxl",
  "openxlsx",
  "jsonlite",
  "zip",
  "plotly",
  "htmlwidgets",
  "msigdbr",
  "BiocManager"
)

install_if_missing(cran_packages)

if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea", ask = FALSE, update = FALSE)
}

required_packages <- c(cran_packages, "fgsea")
still_missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing) > 0) {
  stop("These packages are still missing: ", paste(still_missing, collapse = ", "), call. = FALSE)
}

message("Package installation check complete.")
