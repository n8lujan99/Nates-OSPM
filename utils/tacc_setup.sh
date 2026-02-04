#!/usr/bin/env bash
set -euo pipefail

# ---- locate project root ----
ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$ROOT_DIR"

# ---- activate virtual environment ----
if [[ ! -d ".venv" ]]; then
    echo "[ERROR] .venv not found. Create it first."
    exit 1
fi

source .venv/bin/activate
echo "[OK] venv activated"

# ---- thread & process control ----
# TACC sets SLURM_CPUS_PER_TASK automatically inside jobs
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

# Disable accidental oversubscription
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

echo "[OK] JULIA_NUM_THREADS=${JULIA_NUM_THREADS}"

# ---- optional: sanity print ----
julia -e 'println("Julia threads = ", Threads.nthreads())'

echo "[READY] Environment initialized"
