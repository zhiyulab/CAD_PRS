#!/bin/bash
# ==============================================================================
# run.sh — Unified runner for the SMuRFless CAD pipeline
#
# Usage:
#   bash run.sh real      # Run on real UKB data → output_real/
#   bash run.sh demo      # Run on synthetic data → output_demo/
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE="$SCRIPT_DIR/docker/run_pipeline.R"
R_LIB="$HOME/R/library"
mkdir -p "$R_LIB"
export R_LIBS_USER="$R_LIB"

MODE="${1:-}"

if [ "$MODE" != "real" ] && [ "$MODE" != "demo" ]; then
  echo "Usage: bash run.sh <real|demo>"
  echo ""
  echo "  real  — Run on original UKB data (data/) → output_real/"
  echo "  demo  — Run on synthetic data (data_synthetic/) → output_demo/"
  exit 1
fi

if [ "$MODE" = "real" ]; then
  DATA_DIR="$SCRIPT_DIR/data"
  OUTPUT_DIR="$SCRIPT_DIR/output_real"
  LABEL="Real UKB Data"
else
  DATA_DIR="$SCRIPT_DIR/data_synthetic"
  OUTPUT_DIR="$SCRIPT_DIR/output_demo"
  LABEL="Synthetic Demo Data"
fi

if [ ! -d "$DATA_DIR" ]; then
  echo "ERROR: Data directory not found: $DATA_DIR"
  exit 1
fi

echo "============================================"
echo "  SMuRFless CAD Pipeline"
echo "  Mode:   $LABEL"
echo "  Data:   $DATA_DIR"
echo "  Output: $OUTPUT_DIR"
echo "============================================"

# Install missing packages
Rscript -e "
pkgs <- c('data.table','readr','tableone','boot','speedglm','purrr','tibble',
          'dplyr','tidyr','ggplot2','ggrepel','patchwork','ggh4x','forcats','openxlsx','R.utils')
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
if (length(missing) > 0) {
  cat('Installing:', paste(missing, collapse=', '), '\n')
  install.packages(missing, repos='https://cloud.r-project.org', lib='$R_LIB')
} else {
  cat('All R packages available.\n')
}
"

echo ""
mkdir -p "$OUTPUT_DIR"

Rscript "$PIPELINE" \
  --data_dir="$DATA_DIR" \
  --output_dir="$OUTPUT_DIR"

echo ""
echo "============================================"
echo "  Done! [$LABEL]"
echo "  Outputs: $OUTPUT_DIR"
echo "============================================"
