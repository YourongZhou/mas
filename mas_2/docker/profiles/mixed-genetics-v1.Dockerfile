FROM rocker/r-ver:4.4.1

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONUNBUFFERED=1 \
    FUSION_PATH=/opt/fusion_twas

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-venv \
        build-essential \
        git \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --no-cache-dir --upgrade pip setuptools wheel \
    && python3 -m pip install --no-cache-dir \
        "pandas>=2.1" \
        "numpy>=1.26" \
        "scipy>=1.11" \
        "statsmodels>=0.14" \
        "requests>=2.31" \
        "seaborn>=0.13" \
        "plotnine>=0.13" \
        "adjustText>=0.8" \
        metaxcan

RUN R -q -e "options(repos = c(CRAN='https://cloud.r-project.org')); install.packages(c('BiocManager','optparse','RColorBrewer','data.table','glmnet','reshape','ggplot2','Rcpp','RcppArmadillo'))" \
    && git clone --depth 1 https://github.com/gusevlab/fusion_twas.git "$FUSION_PATH"

ENV PATH="$FUSION_PATH:${PATH}"

WORKDIR /app
CMD ["sleep", "infinity"]
