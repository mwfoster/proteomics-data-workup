# Proteomics Data Workup

## Active DuckDB/RDS projects

The current app uses DuckDB project files (`.duckdb` or `.db`) as its primary save format. Enter an active project path in the **Make metadata** tab to create or reopen a project, save immediately, or enable debounced autosave. A `.rds` path uses the matching RDS fallback when DuckDB is unavailable.

Projects retain metadata, editable SPQC assignments, processed protein tables, batch-corrected results, statistics, and selections across tabs. Replacing one metadata source preserves metadata-independent protein processing while invalidating sample mapping, batch correction, and statistics that need to be recalculated.

Project ZIP exports include uploaded source files and, when DuckDB is installed, an embedded `project_cache.duckdb`. Legacy SQLite project catalogs are left untouched and are not silently converted; use the earlier app version if one must be reopened.

This folder contains a portable Shiny app for building proteomics supplementary tables, metadata, PCA plots, CV plots, volcano plots, feature-level views, box plots, and run-identification summaries.

## Files in this folder

- `app.R` - the Shiny app.
- `run_app.R` - launches the app from this folder.
- `install_packages.R` - installs the R packages used by the app.
- `README.md` - this guide.

## Quick Start

1. Install R if needed.
2. Open R or RStudio.
3. Set the working directory to this folder.
4. Run:

```r
source("install_packages.R")
shiny::runApp(".", launch.browser = TRUE)
```

You can also run:

```r
source("run_app.R")
```

## Required R Packages

The app uses:

- `shiny`
- `ggplot2`
- `DT`
- `dplyr`
- `stringr`
- `missMDA`
- `FactoMineR`
- `svglite`
- `readxl`
- `openxlsx`
- `jsonlite`
- `zip`
- `plotly`
- `htmlwidgets`
- `msigdbr`
- `BiocManager`

If packages are missing, run:

```r
source("install_packages.R")
```

## Typical Inputs

The app lets you choose files from any location on your computer. Common inputs include:

- Condition setup TSV
- Run-order TSV
- Order sample-details workbook (`.xlsx`)
- Protein group report without imputation
- Protein group report with imputation
- CV distribution table, or protein reports for calculating CVs
- Identification overview TSV
- Run identifications TSV for precursors
- Run identifications TSV for protein groups

## Sharing Data With Another User

Inside the app, use the project ZIP export option on the `Make metadata` tab if you want to share the uploaded data files and app settings with another person. The recipient can open that ZIP from the same tab.

For privacy, this share folder does not include project data files by default. Add data files manually only if they are approved for sharing.

## Notes

- Excel workbook export requires `openxlsx`.
- Interactive volcano export requires `plotly`, `htmlwidgets`, and `zip`.
- Missing values are exported as `NaN` where supported.
- Protein table sample measurement headers can be renamed from one or more metadata variables. Selected variables are joined with underscores, such as `Condition_SampleName` or `Condition_Replicate`.
- Protein tables can be exported directly from their preview area as CSV files. CV columns can be placed before or after statistics columns.
- The `Batch correction` tab uses optional `HarmonizR` / ComBat correction and defaults to the imputed protein table (`Table S3`).
- The `Make condition setup` tab can build a condition setup file from a sample-details workbook and an optional existing condition setup template.
- If sample detail names include embedded replicate numbers such as `CKD_01`, `CKD-02`, or `Control 3`, the app strips the replicate number from the condition/label and writes the number into the `Replicate` column.
