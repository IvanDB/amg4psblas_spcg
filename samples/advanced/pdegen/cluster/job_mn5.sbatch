#!/bin/bash
#SBATCH --job-name=sstepcg
#SBATCH --output=sstepcg_%j.out
#SBATCH --error=sstepcg_%j.err
#SBATCH --ntasks-per-node=112          # MN5-GPP: 2 x Xeon 8480+ = 112 cores/node
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --partition=gpp
#SBATCH --qos=gp_debug
#SBATCH --account=ehpc859
#
# Node count is given at submit time:   sbatch --nodes=4 job_mn5.sbatch
#
#   MODE=converge   run to convergence, gives the time to solution (default)
#   MODE=fixed      equal budget of 200 matvecs, gives the cost per iteration
#   PROFILER=extrae produce a Paraver trace (no rebuild needed)
#   PROFILER=none   plain timing run (default)
#   IDIM=100        override the problem size, for validating the pipeline
#                   cheaply before committing to the production size
#
# gp_debug is capped in nodes and wall time. For the production sizes override
# on the command line, which beats the directives in this file:
#   sbatch --nodes=16 --qos=gp_bsccs --time=02:00:00 job_mn5.sbatch

set -e

BIN=../runs/amg_d_pde3d                # ONE method per run, chosen by KMETHD in the .inp
MODE=${MODE:-converge}
PROFILER=${PROFILER:-none}
IDIM=${IDIM:-}                         # empty: use the size written in the .inp
SRUN_MPI=${SRUN_MPI:-pmix}             # OpenMPI 5 needs PMIx. Without it srun starts
                                       # N independent singletons instead of one job.
                                       # Valid values: srun --mpi=list

module purge
module load bsc/1.0
module load gcc/12.3.0
module load ucx/1.16.0-gcc
module load openmpi/5.0.5-gcc
module load openblas/0.3.27-gcc
if [ "$PROFILER" = "extrae" ]; then
  module load extrae || module load EXTRAE || {
    echo "ERROR: no extrae module. Look for the right name with: module avail extrae"; exit 1; }
fi
module list

if [ "$PROFILER" = "extrae" ] && [ -z "${EXTRAE_HOME:-}" ]; then
  echo "ERROR: EXTRAE_HOME is not set by the module, cannot locate libmpitracef.so"
  exit 1
fi

# The binary is linked against the OpenMPI the build job used. Loading a
# different MPI here produces a PMI mismatch: srun starts N singletons, the
# run looks normal and then dies of memory. Check before spending a queue slot.
echo "=== MPI the binary is linked against:"
if ! ldd "$BIN" | grep -E 'libmpi|libmpi_mpifh|libmpifort'; then
  echo "ERROR: $BIN does not resolve its MPI library. Loaded modules do not"
  echo "       match the ones build_mn5.sbatch used."
  exit 1
fi
if ldd "$BIN" | grep -q "not found"; then
  echo "ERROR: unresolved libraries:"; ldd "$BIN" | grep "not found"
  exit 1
fi

mkdir -p results inputs
NP=$(( SLURM_NNODES * SLURM_NTASKS_PER_NODE ))
TAG="n${SLURM_NNODES}_p${NP}_${MODE}"
[ -n "$IDIM" ] && TAG="${TAG}_d${IDIM}"

case "$MODE" in
  converge) SUF="" ;;
  fixed)    SUF="_fixed" ;;
  *) echo "invalid MODE: $MODE"; exit 1 ;;
esac

run_one () {                           # $1 = cg | sstep
  local m=$1
  local out="results/${m}_${TAG}"
  echo
  echo "=== $m   nodes=$SLURM_NNODES  ranks=$NP  mode=$MODE  profiler=$PROFILER"

  if [ "$PROFILER" = "extrae" ]; then
    export EXTRAE_CONFIG_FILE="$SLURM_SUBMIT_DIR/extrae.xml"
    export EXTRAE_DIR="${out}.extrae"
    mkdir -p "$EXTRAE_DIR"
    # Fortran MPI interposition. No instrumented binary required.
    export LD_PRELOAD="${EXTRAE_HOME}/lib/libmpitracef.so"
  fi

  # An IDIM override writes a copy of the input rather than editing it in place.
  local inp="pde3d_${m}${SUF}.inp"
  if [ -n "$IDIM" ]; then
    local gen="inputs/${m}${SUF}_d${IDIM}.inp"
    awk -v d="$IDIM" 'NR==3{sub(/^[0-9]+/, d)} {print}' "$inp" > "$gen"
    inp="$gen"
  fi

  srun --mpi=$SRUN_MPI --ntasks=$NP $BIN "$inp" > "${out}.txt" 2>&1

  # A PMI mismatch does not fail: every rank runs as a standalone MPI job,
  # each one allocating the whole problem. It shows up as a warning and then
  # as an out-of-memory kill, so catch it explicitly.
  if grep -q "singleton" "${out}.txt"; then
    echo "ERROR: ranks did not form a communicator, they ran as singletons."
    echo "       Check the available PMI plugins with: srun --mpi=list"
    echo "       then re-submit with SRUN_MPI=<name>."
    exit 1
  fi
  unset LD_PRELOAD

  if [ "$PROFILER" = "extrae" ]; then
    # Extrae merges into trace.prv on exit when <merge enabled="yes">.
    mv -f trace.prv trace.pcf trace.row "$EXTRAE_DIR"/ 2>/dev/null || true
    echo "--- trace in $EXTRAE_DIR"
  fi

  grep -E "Krylov method|step size|Iterations to convergence|Time to solve|Time per iteration|Residual 2-norm|Linear system size" \
       "${out}.txt" || true
}

run_one cg
run_one sstep

echo
echo "=== DONE ==="
