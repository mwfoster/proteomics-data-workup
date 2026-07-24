options(repos = c(CRAN = "https://cloud.r-project.org"))

r_version_lib <- paste(R.version$major, sub("\\..*", "", R.version$minor), sep = ".")
default_user_lib <- file.path(path.expand("~"), "R", "library", r_version_lib)
user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (!nzchar(user_lib) || normalizePath(user_lib, mustWork = FALSE) == normalizePath(file.path(path.expand("~"), "R", "library"), mustWork = FALSE)) {
  user_lib <- default_user_lib
}
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
base_lib <- file.path(R.home(), "library")
Sys.setenv(R_LIBS_USER = user_lib, R_LIBS_SITE = "", R_LIBS = user_lib)
try(assign(".Library.site", character(0), envir = baseenv()), silent = TRUE)
.libPaths(unique(c(user_lib, base_lib)))

message("R version: ", R.version.string)
message("R executable: ", file.path(R.home("bin"), "R"))
message("User library: ", user_lib)
message("Library paths:\n  ", paste(.libPaths(), collapse = "\n  "))

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

install_one <- function(package, min_version = NULL, repos = getOption("repos")) {
  if (package_ok(package, min_version)) {
    message("OK: ", package, " ", as.character(utils::packageVersion(package)))
    return(TRUE)
  }

  version_note <- if (!is.null(min_version)) paste0(" >= ", min_version) else ""
  message("Installing: ", package, version_note)
  tryCatch(
    {
      install.packages(
        package,
        lib = user_lib,
        repos = repos,
        dependencies = c("Depends", "Imports", "LinkingTo")
      )
      if (package_ok(package, min_version)) {
        message("Installed: ", package, " ", as.character(utils::packageVersion(package)))
        TRUE
      } else {
        warning("Package installed command returned, but package is still unavailable or too old: ", package, call. = FALSE)
        FALSE
      }
    },
    error = function(e) {
      warning("Failed to install ", package, ": ", conditionMessage(e), call. = FALSE)
      FALSE
    }
  )
}

cran_requirements <- c(
  Rcpp = "1.1.0",
  cli = "3.6.5",
  fansi = "1.0.6",
  utf8 = "1.2.6",
  rlang = "1.1.6",
  vctrs = "0.6.5",
  lifecycle = "1.0.4",
  glue = "1.8.0",
  pillar = "1.11.0",
  tibble = "3.3.0",
  pkgconfig = "2.0.3",
  data.table = "1.17.8",
  htmltools = "0.5.8",
  textshaping = "1.0.3",
  systemfonts = "1.2.3",
  cpp11 = "0.5.2",
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

for (package in names(cran_requirements)) {
  min_version <- cran_requirements[[package]]
  install_one(package, if (nzchar(min_version)) min_version else NULL)
}

if (package_ok("BiocManager") && !package_ok("fgsea")) {
  tryCatch(
    {
      message("Installing optional Bioconductor package: fgsea")
      BiocManager::install("fgsea", ask = FALSE, update = FALSE)
    },
    error = function(e) {
      warning(
        "Optional package 'fgsea' could not be installed. The app will still run, ",
        "but the GSEA tab will require fgsea. Details: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

required_packages <- setdiff(names(cran_requirements), c("Rcpp", "cli", "fansi", "utf8", "rlang", "vctrs", "lifecycle", "glue", "pillar", "tibble", "pkgconfig", "data.table", "htmltools", "textshaping", "systemfonts", "cpp11"))
still_missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing) > 0) {
  stop("These packages are still missing: ", paste(still_missing, collapse = ", "), call. = FALSE)
}

if (utils::packageVersion("FactoMineR") < package_version("2.16")) {
  stop("FactoMineR >= 2.16 is required; installed version is ", utils::packageVersion("FactoMineR"), call. = FALSE)
}

message("Package installation check complete.")
