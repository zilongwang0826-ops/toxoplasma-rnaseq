#!/usr/bin/env bash
# Run from the repository root in the documented R/Python plotting environment.
set -euo pipefail
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
Rscript scripts/figures/redesign_figures.R
python scripts/figures/export_vector_assets.py
python scripts/build_research_report.py
python scripts/figures/validate_redesign.py
