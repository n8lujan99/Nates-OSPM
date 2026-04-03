"""
LS6 requires:
- python/3.12.11
- gcc/12.2.0
- Julia 1.10.x in scratch
- XDG_CACHE_HOME and SLURM_TMPDIR redirected to scratch

Run bootstrap_ospm_ls6.sh once.
Never run utils/setup during normal execution.
"""

#!/usr/bin/env bash
set -euo pipefail

# ---------- hygiene (prevents quota disasters) ----------
export XDG_CACHE_HOME=/scratch/$USER/.cache
export XDG_DATA_HOME=/scratch/$USER/.local
export SLURM_TMPDIR=/scratch/$USER/.slurm
mkdir -p $XDG_CACHE_HOME $XDG_DATA_HOME $SLURM_TMPDIR

# ---------- modules ----------
module purge
module load python/3.12.11
module load gcc/12.2.0

# ---------- paths ----------
ROOT=/scratch/$USER/Nates-OSPM
JULIA_DIR=/scratch/$USER/julia-1.10.10
JULIA_TAR=julia-1.10.10-linux-x86_64.tar.gz

cd /scratch/$USER

# ---------- clone repo ----------
if [ ! -d "$ROOT" ]; then
    git clone git@github.com:n8lujan99/Nates-OSPM.git
fi
cd "$ROOT"

# ---------- python venv ----------
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r utils/_Setup/requirements.txt

# ---------- julia install ----------
if [ ! -d "$JULIA_DIR" ]; then
    cd /scratch/$USER
    wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/$JULIA_TAR
    tar -xzf $JULIA_TAR
fi

# ---------- julia environment ----------
export JULIA_DEPOT_PATH=/scratch/$USER/.julia
export OSPM_JULIA_EXE=$JULIA_DIR/bin/julia
mkdir -p $JULIA_DEPOT_PATH

# ---------- rebuild PyCall (one time) ----------
export PYTHON=$(which python)

env -u LD_LIBRARY_PATH -u LD_PRELOAD \
PYTHON="$PYTHON" JULIA_DEPOT_PATH="$JULIA_DEPOT_PATH" \
$OSPM_JULIA_EXE -e '
using Pkg
Pkg.add("PyCall")
Pkg.build("PyCall")
'

# ---------- verification ----------
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
JULIA_DEPOT_PATH="$JULIA_DEPOT_PATH" \
$OSPM_JULIA_EXE -e '
using PyCall
println("PyCall bound to: ", PyCall.python)
'

echo "Bootstrap complete."
echo "Next time, just:"
echo "  module load python/3.12.11 gcc/12.2.0"
echo "  cd $ROOT"
echo "  source .venv/bin/activate"
echo "  export OSPM_USE_JULIA=1"
echo "  bash utils/start"
