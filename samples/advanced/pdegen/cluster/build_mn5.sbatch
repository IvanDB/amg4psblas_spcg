#!/bin/bash
#SBATCH --job-name=build_sstep
#SBATCH --output=build_sstep_%j.out
#SBATCH --error=build_sstep_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16         # cores for make -j
#SBATCH --time=01:00:00
#SBATCH --qos=gp_debug
#SBATCH --account=ehpc859

# ============================================================================
# Builds PSBLAS, AMG4PSBLAS and the pdegen samples inside a job, so the login
# node never runs a long compile.
#
#   sbatch build_mn5.sbatch
#
# Override defaults if the tree lives elsewhere:
#   sbatch --export=ALL,PSBLAS_SRC=/path/to/psblas3_spcg,AMG_SRC=/path/to/amg4psblas_spcg build_mn5.sbatch
#
# With no overrides both repositories are located automatically, starting from
# the submission directory and falling back to a search under $HOME.
#
# Extrae needs no instrumented binaries, it interposes on MPI
# at run time. This is a plain optimised build.
# ============================================================================

set -e

# ---------------------------------------------------------------- parameters
# The two repositories do not have to be siblings, and do not have to sit in
# any particular place: each one is located independently. Set them explicitly
# to skip the search:
#   sbatch --export=ALL,PSBLAS_SRC=...,AMG_SRC=... build_mn5.sbatch
PSBLAS_SRC=${PSBLAS_SRC:-${SRC:+$SRC/psblas3_spcg}}
AMG_SRC=${AMG_SRC:-${SRC:+$SRC/amg4psblas_spcg}}
PREFIX=${PREFIX:-$HOME/opt/psblas3-sstep}      # PSBLAS install prefix, must be writable
MAKE_J=${MAKE_J:-16}
SKIP_CONFIGURE=${SKIP_CONFIGURE:-0}            # set to 1 to reuse an existing Make.inc

# Locate a repository: honour an explicit setting, then look next to the
# submission directory, then walk down from $HOME.
locate_repo () {                               # $1 = given path, $2 = directory name
  local given=$1 name=$2 c
  if [ -n "$given" ] && [ -d "$given" ]; then echo "$(cd "$given" && pwd)"; return 0; fi
  for c in "$SLURM_SUBMIT_DIR" "$SLURM_SUBMIT_DIR/.." "$SLURM_SUBMIT_DIR/../.." \
           "$SLURM_SUBMIT_DIR/../../.." "$SLURM_SUBMIT_DIR/../../../.." \
           "$SLURM_SUBMIT_DIR/../../../../.."; do
    [ -d "$c/$name" ] && { echo "$(cd "$c/$name" && pwd)"; return 0; }
    [ "$(basename "$(cd "$c" 2>/dev/null && pwd)" 2>/dev/null)" = "$name" ] && \
      { echo "$(cd "$c" && pwd)"; return 0; }
  done
  c=$(find "$HOME" -maxdepth 5 -type d -name "$name" -print -quit 2>/dev/null)
  [ -n "$c" ] && { echo "$c"; return 0; }
  return 1
}

# ------------------------------------------------------------------- modules
module purge
module load bsc/1.0
module load gcc/12.3.0
module load ucx/1.16.0-gcc
module load openmpi/5.0.5-gcc
module load openblas/0.3.27-gcc
module list

PSBLAS_SRC=$(locate_repo "$PSBLAS_SRC" psblas3_spcg) || {
  echo "ERROR: psblas3_spcg not found. Re-submit with:"
  echo "  sbatch --export=ALL,PSBLAS_SRC=/path/to/psblas3_spcg,AMG_SRC=/path/to/amg4psblas_spcg $0"
  exit 1; }
AMG_SRC=$(locate_repo "$AMG_SRC" amg4psblas_spcg) || {
  echo "ERROR: amg4psblas_spcg not found. Re-submit with:"
  echo "  sbatch --export=ALL,PSBLAS_SRC=/path/to/psblas3_spcg,AMG_SRC=/path/to/amg4psblas_spcg $0"
  exit 1; }

echo "=== PSBLAS src  : $PSBLAS_SRC"
echo "=== AMG src     : $AMG_SRC"
echo "=== install to  : $PREFIX"
echo "=== make -j$MAKE_J"
echo "=== compilers   : $(which mpifort) / $(which mpicc)"
mpifort --version | head -1

# -------------------------------------------------------------------- PSBLAS
echo
echo "########## 1/3  PSBLAS"
cd "$PSBLAS_SRC"
if [ "$SKIP_CONFIGURE" != "1" ]; then
  # No --enable-cuda: GPP has no GPUs and the measurements are CPU-side.
  ./configure --prefix="$PREFIX" FC=gfortran CC=gcc CXX=g++
fi
make -j"$MAKE_J"
make install
echo "--- installed:"; ls -1 "$PREFIX/lib"

# ---------------------------------------------------------------- AMG4PSBLAS
echo
echo "########## 2/3  AMG4PSBLAS"
cd "$AMG_SRC"
if [ "$SKIP_CONFIGURE" != "1" ]; then
  ./configure --with-psblas="$PREFIX" FC=gfortran CC=gcc CXX=g++
fi
make -j"$MAKE_J"
echo "--- built:"; ls -1 lib

# ------------------------------------------------------------------- samples
echo
echo "########## 3/3  pdegen samples"
cd "$AMG_SRC/samples/advanced/pdegen"
make -j"$MAKE_J"

echo
echo "########## result"
for b in amg_d_pde3d amg_d_pde3d_sstepbatch; do
  [ -x "runs/$b" ] && echo "  OK       runs/$b" || { echo "  MISSING  runs/$b"; exit 1; }
done
echo
echo "=== BUILD COMPLETE ==="
echo "Next:  cd $AMG_SRC/samples/advanced/pdegen/cluster"
echo "       sbatch --nodes=1 job_mn5.sbatch"
