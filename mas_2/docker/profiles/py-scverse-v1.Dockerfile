FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential git \
    && rm -rf /var/lib/apt/lists/*

RUN python -m pip install --no-cache-dir --upgrade pip setuptools wheel \
    && python -m pip install --no-cache-dir \
        --index-url https://download.pytorch.org/whl/cpu \
        "torch==2.5.1" \
    && python -m pip install --no-cache-dir \
        "numpy>=1.26" \
        "pandas>=2.1" \
        "scipy>=1.11" \
        "scikit-learn>=1.4" \
        "scikit-misc>=0.1" \
        "statsmodels>=0.14" \
        "matplotlib>=3.8" \
        "seaborn>=0.13" \
        "plotnine>=0.13" \
        "reportlab>=4.0" \
        "scanpy>=1.10" \
        "anndata>=0.10" \
        "adjustText>=0.8" \
        "scrublet>=0.2.3" \
        "squidpy>=1.4" \
        "harmonypy>=0.0.10" \
        "celltypist>=1.6" \
        "pydeseq2>=0.4" \
        "scvelo>=0.3" \
        "cellrank>=2.0" \
    && python -m pip install --no-cache-dir \
        "scvi-tools>=1.1"

WORKDIR /app
CMD ["sleep", "infinity"]
