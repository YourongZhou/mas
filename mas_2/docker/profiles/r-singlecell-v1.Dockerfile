FROM rocker/r-ver:4.4.1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        git \
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
        libglpk40 \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e "options(repos = c(CRAN='https://cloud.r-project.org')); install.packages(c('BiocManager','remotes','Seurat','ggplot2','ggprism','dplyr','patchwork','harmony','NMF','circlize','ggalluvial','rmarkdown','DoubletFinder','SoupX')); BiocManager::install(c('ComplexHeatmap','DESeq2','muscat','SingleR','celldex'), ask = FALSE, update = FALSE); remotes::install_github('satijalab/seurat-data'); remotes::install_github('jinworks/CellChat'); remotes::install_github('immunogenomics/presto')"

WORKDIR /app
CMD ["sleep", "infinity"]
