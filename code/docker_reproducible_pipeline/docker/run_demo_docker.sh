#!/bin/bash
# ==============================================================================
# run_demo_docker.sh — Runs INSIDE the Docker container
#
# Executes the full analysis pipeline on bundled synthetic data.
# Outputs are written to /output (bind-mount a host directory there).
# ==============================================================================

set -euo pipefail

echo "============================================"
echo "  SMuRFless CAD Analysis — Docker Demo"
echo "============================================"
echo ""
echo "R version: $(R --version | head -1)"
echo "Synthetic data: /pipeline/data_synthetic"
echo "Output dir:     /output"
echo ""

echo "Running analysis pipeline on synthetic data..."
echo ""

Rscript /pipeline/run_pipeline.R \
  --data_dir=/pipeline/data_synthetic \
  --output_dir=/output

echo ""
echo "============================================"
echo "  Demo complete. Outputs in /output/"
echo "============================================"
echo ""
echo "Output files:"
find /output -type f | sort
