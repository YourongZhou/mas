FROM rocker/r-ver:4.4.1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e "options(repos = c(CRAN='https://cloud.r-project.org')); install.packages(c('BiocManager','ggplot2','ggprism','ggrepel','survival','survminer','scales')); BiocManager::install(c('DESeq2','apeglm','pasilla','airway','RTCGA.clinical'), ask = FALSE, update = FALSE)"

WORKDIR /app
CMD ["sleep", "infinity"]
