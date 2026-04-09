#!/bin/bash
# ==============================================================================
# run_docker_real.sh — Run the pipeline on real UKB data via Docker
#
# Usage:
#   bash run_docker_real.sh              # expects data/ in current directory
#   bash run_docker_real.sh /path/to/data
#
# Your data directory must follow the structure documented in README.md.
# The data is mounted read-only — nothing is copied into the image.
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
IMAGE_NAME="smurfless-cad-pipeline"
DATA_DIR="${1:-$SCRIPT_DIR/data}"
OUTPUT_DIR="$SCRIPT_DIR/output_real"

if [ ! -d "$DATA_DIR/ukbb_pheno" ] || [ ! -d "$DATA_DIR/prs_58" ]; then
  echo "ERROR: Data directory not found or incomplete: $DATA_DIR"
  echo ""
  echo "Expected structure:"
  echo "  $DATA_DIR/ukbb_pheno/ukbb_PhenoFile.ALL_500k.*.txt"
  echo "  $DATA_DIR/intermediate_cad/20250615_pheno_noGP.tsv.gz"
  echo "  $DATA_DIR/prs_58/PC4.*.adjNormPRS.txt.gz"
  echo "  $DATA_DIR/le8_scores/updated_LE8score_wComposite.txt"
  echo "  $DATA_DIR/sliang/  (pce_zy_2.tsv, family_history_final.csv, etc.)"
  echo ""
  echo "See README.md for full list of required files."
  exit 1
fi

echo "============================================"
echo "  SMuRFless CAD — Real Data Pipeline"
echo "============================================"
echo ""
echo "  Data (read-only): $DATA_DIR"
echo "  Output:           $OUTPUT_DIR"
echo ""

# Build image (reuses cached layers from demo build)
echo ">>> Building Docker image..."
docker build -f "$SCRIPT_DIR/docker/Dockerfile" -t "$IMAGE_NAME" "$SCRIPT_DIR"
echo ""

# Run with real data mounted, overriding the bundled synthetic data
echo ">>> Running pipeline on real data..."
mkdir -p "$OUTPUT_DIR"
docker run --rm \
  -v "$DATA_DIR:/data:ro" \
  -v "$OUTPUT_DIR:/output" \
  --memory=16g \
  --cpus=4 \
  "$IMAGE_NAME" \
  Rscript /pipeline/run_pipeline.R --data_dir=/data --output_dir=/output

echo ""
echo "============================================"
echo "  Done! Outputs in: $OUTPUT_DIR"
echo "============================================"
