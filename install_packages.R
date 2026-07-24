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

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea", ask = FALSE, update = FALSE)
}

message("Package installation check complete.")
