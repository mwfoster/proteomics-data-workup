# Proteomics Data Workup

Portable R/Shiny application for processing, analyzing, and visualizing proteomics data.

## Install

1. Install R if needed.
2. Open R or RStudio.
3. Set the working directory to this folder.
4. Install the required R packages:

```r
source("install_packages.R")
```

## Run

Launch the application with:

```r
shiny::runApp(".", launch.browser = TRUE)
```

Alternatively, run:

```r
source("run_app.R")
```

DuckDB project files can be opened from the app. Use the app's project download control to save a new or updated project to your computer.
