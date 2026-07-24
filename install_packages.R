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

package_ok <- function(package, min_version = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    return(FALSE)
  }
  if (!is.null(min_version) && utils::packageVersion(package) < package_version(min_version)) {
    return(FALSE)
  }
  TRUE
}

install_if_needed <- function(requirements, repos = "https://cloud.r-project.org") {
  needed <- names(requirements)[!vapply(names(requirements), function(package) {
    min_version <- requirements[[package]]
    package_ok(package, if (nzchar(min_version)) min_version else NULL)
  }, logical(1))]
  if (length(needed) > 0) {
    install.packages(
      needed,
      repos = repos,
      dependencies = c("Depends", "Imports", "LinkingTo")
    )
  }
}

cran_requirements <- c(
  shiny = "",
  ggplot2 = "",
  DT = "",
  dplyr = "",
  stringr = "",
  missMDA = "",
  FactoMineR = "2.16",
  svglite = "",
  readxl = "",
  openxlsx = "",
  jsonlite = "",
  zip = "",
  plotly = "",
  htmlwidgets = "",
  msigdbr = "",
  BiocManager = ""
)

install_if_needed(cran_requirements)

if (!package_ok("fgsea")) {
  bioc_args <- list(pkgs = "fgsea", ask = FALSE, update = FALSE)
  if (getRversion() >= "4.6.0") {
    bioc_args$version <- "3.23"
  }
  do.call(BiocManager::install, bioc_args)
}

required_packages <- c(names(cran_requirements), "fgsea")
still_missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing) > 0) {
  stop("These packages are still missing: ", paste(still_missing, collapse = ", "), call. = FALSE)
}

if (utils::packageVersion("FactoMineR") < package_version("2.16")) {
  stop("FactoMineR >= 2.16 is required; installed version is ", utils::packageVersion("FactoMineR"), call. = FALSE)
}

message("Package installation check complete.")
