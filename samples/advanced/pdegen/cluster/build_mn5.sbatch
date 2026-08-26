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
#   sbatch --export=ALL,SRC=$HOME/sstep,PREFIX=$HOME/opt/psblas build_mn5.sbatch
#
# No Score-P here: Extrae needs no instrumented binaries, it interposes on MPI
# at run time. This is a plain optimised build.
# ============================================================================

set -e

# ---------------------------------------------------------------- parameters
SRC=${SRC:-$HOME/sstep}                        # holds psblas3_spcg/ and amg4psblas_spcg/
PREFIX=${PREFIX:-$HOME/opt/psblas3-sstep}      # PSBLAS install prefix, must be writable
MAKE_J=${MAKE_J:-16}
SKIP_CONFIGURE=${SKIP_CONFIGURE:-0}            # set to 1 to reuse an existing Make.inc

# ------------------------------------------------------------------- modules
module purge
module load bsc/1.0
module load gcc/12.3.0
module load ucx/1.16.0-gcc
module load openmpi/5.0.5-gcc
module load openblas/0.3.27-gcc
module list

echo "=== source tree : $SRC"
echo "=== install to  : $PREFIX"
echo "=== make -j$MAKE_J"
echo "=== compilers   : $(which mpifort) / $(which mpicc)"
mpifort --version | head -1

test -d "$SRC/psblas3_spcg"    || { echo "psblas3_spcg not found under $SRC";    exit 1; }
test -d "$SRC/amg4psblas_spcg" || { echo "amg4psblas_spcg not found under $SRC"; exit 1; }

# -------------------------------------------------------------------- PSBLAS
echo
echo "########## 1/3  PSBLAS"
cd "$SRC/psblas3_spcg"
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
cd "$SRC/amg4psblas_spcg"
if [ "$SKIP_CONFIGURE" != "1" ]; then
  ./configure --with-psblas="$PREFIX" FC=gfortran CC=gcc CXX=g++
fi
make -j"$MAKE_J"
echo "--- built:"; ls -1 lib

# ------------------------------------------------------------------- samples
echo
echo "########## 3/3  pdegen samples"
cd "$SRC/amg4psblas_spcg/samples/advanced/pdegen"
make -j"$MAKE_J"

echo
echo "########## result"
for b in amg_d_pde3d amg_d_pde3d_sstepbatch; do
  [ -x "runs/$b" ] && echo "  OK       runs/$b" || { echo "  MISSING  runs/$b"; exit 1; }
done
echo
echo "=== BUILD COMPLETE ==="
echo "Next:  cd $SRC/amg4psblas_spcg/samples/advanced/pdegen/cluster"
echo "       sbatch --nodes=1 job_mn5.sbatch"
