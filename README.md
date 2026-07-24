# Proteomics Data Workup

A Shiny app for proteomics data workup, including metadata construction, protein supplementary table export, PCA, CV plots, volcano plots, gene set enrichment, feature plots, boxplots, run-identification summaries, and Evosep run queue generation.

## Repository Layout

- `app.R` - main Shiny app.
- `run_app.R` - local launcher for development/testing.
- `install_packages.R` - installs required R packages.
- `.gitignore` - excludes uploaded study data, generated exports, R session state, and local caches.
- `docs/ubuntu-shiny-server.md` - Ubuntu VM deployment notes.

## Local Quick Start

Install R, then from the repository root run:

```r
source("install_packages.R")
source("run_app.R")
```

Or directly:

```r
shiny::runApp(".", launch.browser = TRUE)
```

## Main Features

- Build condition setup files from sample-detail workbooks.
- Build metadata from condition setup, run order, and sample-detail files.
- Export supplementary Excel workbooks with metadata and protein tables.
- Rename protein report headers using metadata-derived labels.
- Calculate CVs and basic protein statistics, including paired comparisons when replicate pairing is valid.
- Create PCA plots, PCA loading summaries, volcano plots, GSEA summaries, feature plots, and boxplots.
- Build Evosep `.csl` run queue files from manifest/sample-detail workbooks and a template `.csl`.
- Save and reopen project ZIP bundles containing uploaded data and app settings.

## Data Handling

Do not commit uploaded study files, exported workbooks, project ZIP bundles, raw proteomics reports, or PHI/sensitive data to this repository. Use the app's project ZIP export only for approved sharing.

## Ubuntu VM Hosting

Recommended VM baseline for larger proteomics tables:

- Ubuntu 24.04 LTS
- 4 CPU
- 40 GB RAM
- 200 GB disk

See `docs/ubuntu-shiny-server.md` for a basic deployment path using Shiny Server.

## Notes

- Excel export requires `openxlsx`.
- Interactive volcano export requires `plotly`, `htmlwidgets`, and `zip`.
- Gene set enrichment uses `msigdbr` and `fgsea`.
- The app sets Shiny's max upload size to 1 GB in `app.R`.
