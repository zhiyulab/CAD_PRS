# ==============================================================================
# Dockerfile — SMuRFless CAD Primary Analysis Pipeline
#
# Builds a self-contained image with:
#   - R 4.4.1 + pinned package versions
#   - Pre-generated synthetic data (~255 MB)
#   - run_pipeline.R (the unified 7-step analysis)
#
# Quick start (from project root):
#   docker build -f docker/Dockerfile -t smurfless-cad-demo .
#   docker run --rm -v $(pwd)/output_docker:/output smurfless-cad-demo
#
# Mount points:
#   /output  — bind-mount a host directory to collect results
#   /data    — (optional) bind-mount real UKB data to override synthetic data
# ==============================================================================

FROM rocker/r-ver:4.4.1

LABEL maintainer="Broad Institute"
LABEL description="SMuRFless CAD primary analysis pipeline (Docker demo)"
LABEL org.opencontainers.image.source="https://github.com/broad-sliang/smurfless-cad"
LABEL org.opencontainers.image.licenses="MIT"

# --- System dependencies ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# --- Install remotes first, then pin package versions ---
RUN R -e "install.packages('remotes', repos = 'https://cloud.r-project.org')"

RUN R -e " \
  remotes::install_version('data.table',  version = '1.16.0',  repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('readr',       version = '2.1.5',   repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('tableone',    version = '0.13.2',  repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('speedglm',    version = '0.3-5',   repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('boot',        version = '1.3-30',  repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('purrr',       version = '1.0.2',   repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('tibble',      version = '3.2.1',   repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('dplyr',       version = '1.1.4',   repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('tidyr',       version = '1.3.1',   repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('ggplot2',     version = '3.5.1',   repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('ggrepel',     version = '0.9.5',   repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('patchwork',   version = '1.2.0',   repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('ggh4x',       version = '0.2.8',   repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('forcats',     version = '1.0.0',   repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('openxlsx',    version = '4.2.5.2', repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('R.utils',     version = '2.12.3',  repos = 'https://cloud.r-project.org', upgrade = 'never'); \
"

# --- Set up working directory ---
WORKDIR /pipeline

# --- Copy pipeline scripts ---
COPY docker/run_pipeline.R          /pipeline/run_pipeline.R
COPY docker/run_demo_docker.sh      /pipeline/run_demo_docker.sh
RUN chmod +x /pipeline/run_demo_docker.sh

# --- Bundle pre-generated synthetic data ---
COPY data_synthetic/ /pipeline/data_synthetic/

# --- Default entrypoint: run the demo ---
CMD ["/pipeline/run_demo_docker.sh"]
