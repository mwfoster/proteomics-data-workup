FROM rocker/shiny:4.4.2

ENV DEBIAN_FRONTEND=noninteractive \
    R_LIBS_USER=/usr/local/lib/R/site-library \
    R_LIBS_SITE=/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libjpeg-dev \
    libpng-dev \
    libssl-dev \
    libtiff5-dev \
    libxml2-dev \
    pandoc \
    r-cran-dplyr \
    r-cran-factominer \
    r-cran-ggplot2 \
    r-cran-plotly \
    zip \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "print(.libPaths()); stopifnot(requireNamespace('ggplot2', quietly=TRUE), requireNamespace('dplyr', quietly=TRUE), requireNamespace('FactoMineR', quietly=TRUE), requireNamespace('plotly', quietly=TRUE))"

WORKDIR /srv/shiny-server/proteomics

COPY install_packages.R ./install_packages.R
RUN Rscript install_packages.R

COPY . .

RUN chown -R shiny:shiny /srv/shiny-server/proteomics \
    && mkdir -p /var/lib/shiny-server/bookmarks \
    && chown -R shiny:shiny /var/lib/shiny-server

USER shiny

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/proteomics', host='0.0.0.0', port=3838, launch.browser=FALSE)"]
