FROM rocker/shiny:4.4.2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    cmake \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libmagick++-dev \
    pkg-config \
    libharfbuzz-dev \
    libfribidi-dev \
    libjpeg-dev \
    libtiff-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /srv/shiny-server/

# Copy local Shiny app and renv.lock
COPY . .

RUN R -e "install.packages(c('renv', 'remotes'), repos = 'https://cloud.r-project.org')"
RUN R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org'); BiocManager::install(version = '3.20', ask = FALSE)"
RUN R -e "renv::restore(lockfile = '/srv/shiny-server/renv.lock', prompt = FALSE)"

RUN chown -R shiny:shiny /srv/shiny-server

EXPOSE 3838

CMD ["Rscript", "-e", "shiny::runApp('/srv/shiny-server', host = '0.0.0.0', port = 3838)"]
