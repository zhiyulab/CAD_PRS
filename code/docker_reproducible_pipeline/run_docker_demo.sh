#!/bin/bash
# ==============================================================================
# run_demo.sh — One-click demo. Just needs Docker.
#
# Usage:
#   bash run_demo.sh
#
# This builds a Docker container with all dependencies + synthetic data,
# runs the full 7-step analysis, and saves outputs to output_demo/.
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
IMAGE_NAME="smurfless-cad-demo"
OUTPUT_DIR="$SCRIPT_DIR/output_demo"

echo "============================================"
echo "  SMuRFless CAD — Demo Pipeline"
echo "============================================"
echo ""
echo "  Requirements: Docker"
echo "  Runtime:      ~1 minute"
echo "  Output:       $OUTPUT_DIR"
echo ""

# Build
echo ">>> Building Docker image (first time takes ~10 min, cached after)..."
docker build -f "$SCRIPT_DIR/docker/Dockerfile" -t "$IMAGE_NAME" "$SCRIPT_DIR"
echo ""

# Run
echo ">>> Running pipeline on synthetic data..."
mkdir -p "$OUTPUT_DIR"
docker run --rm \
  -v "$OUTPUT_DIR:/output" \
  --memory=4g \
  --cpus=2 \
  "$IMAGE_NAME"

echo ""
echo "============================================"
echo "  Done! Outputs in: $OUTPUT_DIR"
echo "============================================"
echo ""
echo "Output files:"
find "$OUTPUT_DIR" -type f ! -path "*/.cache/*" | sort | while read -r f; do
  echo "  ${f#$OUTPUT_DIR/}"
done
