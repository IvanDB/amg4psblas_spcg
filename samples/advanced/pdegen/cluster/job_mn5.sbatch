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
ITMAX=${ITMAX:-}                       # cap the run in MATRIX-VECTOR PRODUCTS, divided
                                       # by s for the s-step input so both methods do
                                       # the same work. The main lever on trace size.
TRACE_DIR=${TRACE_DIR:-}               # where Extrae writes. Point it at scratch:
                                       # traces are large and GPFS home has a quota.
SWEEPS=${SWEEPS:-}                     # empty: use the value in the .inp. Each extra
                                       # Jacobi sweep is one more spmm per iteration
                                       # and dilutes the share of all-reduces that
                                       # s-step exists to remove.
MPI_DEBUG=${MPI_DEBUG:-0}              # 1: report which PML OpenMPI selected
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
[ -n "$IDIM" ]   && TAG="${TAG}_d${IDIM}"
[ -n "$SWEEPS" ] && TAG="${TAG}_j${SWEEPS}"
[ -n "$ITMAX" ]  && TAG="${TAG}_m${ITMAX}"

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

  # Which point-to-point layer got selected decides whether the run went over
  # the fabric or fell back to a slow path. Worth knowing before trusting any
  # communication number.
  if [ "$MPI_DEBUG" = "1" ]; then
    export OMPI_MCA_pml_base_verbose=10
  fi

  if [ "$PROFILER" = "extrae" ]; then
    export EXTRAE_CONFIG_FILE="$SLURM_SUBMIT_DIR/extrae.xml"
    export EXTRAE_DIR="${TRACE_DIR:+$TRACE_DIR/}${out##*/}.extrae"
    [ -n "$TRACE_DIR" ] || EXTRAE_DIR="${out}.extrae"
    mkdir -p "$EXTRAE_DIR"
    # Fortran MPI interposition. No instrumented binary required.
    export LD_PRELOAD="${EXTRAE_HOME}/lib/libmpitracef.so"
  fi

  # An IDIM override writes a copy of the input rather than editing it in place.
  local inp="pde3d_${m}${SUF}.inp"
  if [ -n "$IDIM" ] || [ -n "$SWEEPS" ] || [ -n "$ITMAX" ]; then
    local gen="inputs/${m}${SUF}${IDIM:+_d$IDIM}${SWEEPS:+_j$SWEEPS}${ITMAX:+_m$ITMAX}.inp"
    # ITMAX counts outer iterations, each worth s matvecs for the s-step method,
    # so divide there to keep the two runs comparable.
    local sval; sval=$(awk 'NR==11{print $1}' "$inp")
    local cap="$ITMAX"
    if [ -n "$ITMAX" ] && [ "$m" = "sstep" ]; then
      cap=$(( ITMAX / sval )); [ "$cap" -lt 1 ] && cap=1
    fi
    awk -v d="$IDIM" -v j="$SWEEPS" -v c="$cap" '
      NR==3  && d != "" { sub(/^[0-9]+/, d) }
      NR==8  && c != "" { sub(/^[0-9]+/, c) }
      NR==21 && j != "" { sub(/^[0-9]+/, j) }
      { print }' "$inp" > "$gen"
    inp="$gen"
    echo "    input: $inp   (itmax $(awk 'NR==8{print $1}' "$gen"), sweeps $(awk 'NR==21{print $1}' "$gen"))"
  fi

  srun --mpi=$SRUN_MPI --ntasks=$NP $BIN "$inp" > "${out}.txt" 2>&1

  # A PMI mismatch does not fail: every rank runs as a standalone MPI job,
  # each one allocating the whole problem. It shows up as a warning and then
  # as an out-of-memory kill, so catch it explicitly.
  if [ "$MPI_DEBUG" = "1" ]; then
    echo "--- point-to-point layer selected:"
    grep -m2 -iE "select.*pml|pml.*select" "${out}.txt" || echo "    (nothing reported)"
    unset OMPI_MCA_pml_base_verbose
  fi

  # A UCX version mismatch degrades the transport silently. On a benchmark whose
  # whole point is communication cost, that invalidates the numbers.
  if grep -q "UCP API version is incompatible" "${out}.txt"; then
    echo "WARNING: UCX version mismatch, the transport may have fallen back."
    grep -m1 "UCP API version is incompatible" "${out}.txt"
    echo "         Look for a newer module with: module avail ucx"
  fi

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
    du -sh "$EXTRAE_DIR" 2>/dev/null
    if [ ! -s "$EXTRAE_DIR/trace.prv" ]; then
      echo "    WARNING: no trace.prv. The merge may have failed; the .mpits are kept,"
      echo "             so it can be redone offline with mpi2prv."
    fi
  fi

  grep -E "Krylov method|step size|Iterations to convergence|Time to solve|Time per iteration|Residual 2-norm|Linear system size" \
       "${out}.txt" || true
}

run_one cg
run_one sstep

echo
echo "=== DONE ==="
