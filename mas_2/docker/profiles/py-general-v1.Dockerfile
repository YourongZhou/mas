FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential git \
    && rm -rf /var/lib/apt/lists/*

RUN python -m pip install --no-cache-dir --upgrade pip setuptools wheel \
    && python -m pip install --no-cache-dir \
        "numpy>=1.26" \
        "pandas>=2.1" \
        "scipy>=1.11" \
        "scikit-learn>=1.4" \
        "matplotlib>=3.8" \
        "seaborn>=0.13" \
        "statsmodels>=0.14" \
        "plotnine>=0.13" \
        "requests>=2.31" \
        "reportlab>=4.0" \
        "adjustText>=0.8"

WORKDIR /app
CMD ["sleep", "infinity"]
