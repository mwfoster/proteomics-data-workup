# Proteomics Data Workup

Portable R/Shiny application for processing, analyzing, and visualizing proteomics data.

## Use the hosted app

The current hosted application is available at [https://rapid-2905.vm.duke.edu/](https://rapid-2905.vm.duke.edu/). No local R installation is needed to use the hosted app. Upload your data or DuckDB project through the app, and download updated projects and results to your computer.

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

DuckDB is the primary project format and supports project autosave. RDS projects remain available as a fallback when DuckDB is unavailable. Legacy SQLite project catalogs are not converted automatically.

## Apply metadata changes

Metadata edits are staged on the **Make metadata** tab. This includes edits to the metadata table, SPQC assignments, sample exclusions, and changes uploaded through a modified metadata file.

1. Make or upload your metadata changes on the **Make metadata** tab.
2. Click **Apply metadata changes** to commit the tab's draft metadata for analysis.
3. Check the relevant downstream tables and plots. Rerun PCA, statistics, or batch correction as needed if the changes affect sample inclusion, grouping, pairing, batch assignments, or model covariates.
4. Save or download the updated DuckDB project to retain the applied changes.

Until you click **Apply metadata changes**, subsequent tabs continue to use the previously applied metadata. Switching tabs does not apply draft edits. Applying metadata changes is not the same as saving a project or rerunning an analysis.

## Run on an Ubuntu VM with Docker

The repository includes a `Dockerfile` and `compose.yaml`. On an Ubuntu VM with Docker and the Docker Compose plugin installed, run:

```bash
git clone https://github.com/mwfoster/proteomics-data-workup.git
cd proteomics-data-workup
docker compose up -d --build
docker compose ps
```

The container listens only on `127.0.0.1:6875` on the VM so it can be published securely through a reverse proxy. See [VM-DEPLOYMENT.md](VM-DEPLOYMENT.md) for nginx, HTTPS, large DuckDB uploads, logs, and update instructions.

## License

This project is available under the [MIT License](LICENSE).
