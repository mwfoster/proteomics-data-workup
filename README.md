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

DuckDB is the primary project format and supports project autosave. RDS projects remain available as a fallback when DuckDB is unavailable. Legacy SQLite project catalogs are not converted automatically.

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
