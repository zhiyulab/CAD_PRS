# Dockerized Reproducible Pipeline

## Quick Start (Demo)

First, download the synthetic data from Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19476250.svg)](https://doi.org/10.5281/zenodo.19476250)

```bash
curl -L -o data_synthetic.zip https://zenodo.org/records/19476250/files/data_synthetic.zip
unzip data_synthetic.zip -d data_synthetic/
```

Then run:

```bash
bash run_docker_demo.sh
```

This builds a Docker image with all dependencies, runs the full pipeline on the synthetic data, and writes results to `output_demo/`. Only requires Docker. Takes about a minute after the initial build (~10 min for R package compilation, cached after that). Uses ~600 MB RAM.

---

## Running on Real UK Biobank Data

Put your UKB data in `data/` (see [Required Input Data](#required-input-data-for-real-ukb-analysis) for the expected layout), then:

```bash
# Via Docker (no R install needed)
bash run_docker_real.sh                    # reads from data/, writes to output_real/
bash run_docker_real.sh /path/to/data      # custom data location

# Locally (needs R 4.4+)
bash run.sh real
bash run.sh demo      # runs on synthetic data instead
```

Data is mounted read-only -- nothing gets copied into the image. Expect ~10 minutes and ~16 GB RAM on the full UKB cohort.

---

## Repository Layout

```
├── run_docker_demo.sh     # One-click demo via Docker
├── run_docker_real.sh     # Run on real UKB data via Docker
├── run.sh                 # Local runner: bash run.sh <real|demo>
│
├── docker/
│   ├── Dockerfile         # R 4.4.1 + pinned packages
│   ├── run_pipeline.R     # Main analysis script
│   └── run_demo_docker.sh # Container entrypoint
│
├── data_synthetic/        # Synthetic data for demo
├── data/                  # Real UKB data goes here (not included)
├── output_demo/           # Demo output
└── output_real/           # Real data output
```

---

## Reproducibility

- Random seed: 42
- R version: 4.4.1 (pinned in Dockerfile)
- 1000 bootstrap resamples
- All package versions pinned via `remotes::install_version()`

Bootstrap CIs are reproducible given the same seed and R version.

---

## Computational Requirements

| | Time | Peak RAM |
|---|---|---|
| Real UKB (~502k) | ~10 min | ~12.5 GB |
| Synthetic demo (~25k) | ~30 sec | ~570 MB |

Recommended: 16 GB RAM for real data, 4 GB for demo.

---

## R Packages

| Package | Version |
|---------|---------|
| data.table | 1.16.0 |
| readr | 2.1.5 |
| tableone | 0.13.2 |
| boot | 1.3-30 |
| speedglm | 0.3-5 |
| purrr | 1.0.2 |
| tibble | 3.2.1 |
| dplyr | 1.1.4 |
| tidyr | 1.3.1 |
| ggplot2 | 3.5.1 |
| ggrepel | 0.9.5 |
| patchwork | 1.2.0 |
| ggh4x | 0.2.8 |
| forcats | 1.0.0 |
| openxlsx | 4.2.5.2 |
| R.utils | 2.12.3 |

All versions pinned in the Dockerfile.

