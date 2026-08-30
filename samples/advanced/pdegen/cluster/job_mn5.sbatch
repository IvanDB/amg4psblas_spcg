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
WEAK=${WEAK:-}                         # unknowns per rank: derives IDIM from the rank
                                       # count, keeping the local problem fixed. Strong
                                       # scaling shrinks it into cache, which favours
                                       # whichever method holds fewer vectors.
ITMAX=${ITMAX:-}                       # cap the run in MATRIX-VECTOR PRODUCTS, divided
                                       # by s for the s-step input so both methods do
                                       # the same work. The main lever on trace size.
TRACE_DIR=${TRACE_DIR:-}               # where Extrae writes. Left empty it looks for
                                       # scratch, since traces are large and GPFS home
                                       # is quota limited and slow for many small files.
STEPS=${STEPS:-}                       # s, the block size. Reductions fall as 1/s but
                                       # the Gram work grows as s, so there is an
                                       # optimum rather than "bigger is better".
SWEEPS=${SWEEPS:-}                     # empty: use the value in the .inp. Each extra
                                       # Jacobi sweep is one more spmm per iteration
                                       # and dilutes the share of all-reduces that
                                       # s-step exists to remove.
REPS=${REPS:-1}                        # repetitions of each method inside the SAME
                                       # allocation. Past 8 nodes a single run varies
                                       # by 15-20% with network contention, so one
                                       # measurement is not a measurement.
PSB_SPMM_TIMERS=${PSB_SPMM_TIMERS:-}   # 1: break the sparse product into halo
                                       # exchange and local computation. Their
                                       # ratio decides whether overlapping or
                                       # deepening the halo is worth anything.
PSB_SPMM_BARRIER=${PSB_SPMM_BARRIER:-} # 1: barrier in front of each phase, so
                                       # imbalance is charged to the barrier and
                                       # not to the exchange behind it. Perturbs
                                       # the run: use it to attribute, not to time.
export PSB_SPMM_TIMERS PSB_SPMM_BARRIER
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

# Locate the Fortran MPI tracing library. Not every extrae module exports
# EXTRAE_HOME, so fall back to searching LD_LIBRARY_PATH for the library
# itself, which is what LD_PRELOAD actually needs.
# EXTRAE_LIB can be set explicitly. The module ships one library per MPI
# flavour and picks the wrong one when it cannot recognise OpenMPI: preloading
# the Intel MPI build into an OpenMPI program dies with "Attempting to use an
# MPI routine before initializing MPICH".
EXTRAE_LIB=${EXTRAE_LIB:-}
if [ "$PROFILER" = "extrae" ] && [ -n "$EXTRAE_LIB" ]; then
  [ -f "$EXTRAE_LIB" ] || { echo "ERROR: EXTRAE_LIB=$EXTRAE_LIB does not exist"; exit 1; }
  echo "=== Extrae library (given): $EXTRAE_LIB"
elif [ "$PROFILER" = "extrae" ]; then
  for c in "${EXTRAE_HOME:-}/lib" "${EXTRAE_ROOT:-}/lib" ${LD_LIBRARY_PATH//:/ }; do
    [ -n "$c" ] && [ -f "$c/libmpitracef.so" ] && { EXTRAE_LIB="$c/libmpitracef.so"; break; }
  done
  if [ -z "$EXTRAE_LIB" ]; then
    echo "ERROR: libmpitracef.so not found. EXTRAE_HOME=[${EXTRAE_HOME:-unset}]"
    echo "       Searched EXTRAE_HOME/lib, EXTRAE_ROOT/lib and LD_LIBRARY_PATH."
    echo "       Locate it with: find /apps -name libmpitracef.so 2>/dev/null | head"
    exit 1
  fi
  echo "=== Extrae library (found): $EXTRAE_LIB"
  case "$EXTRAE_LIB" in
    *impi*|*intel*|*mpich*)
      echo "WARNING: that looks like the Intel MPI build, but this binary uses"
      echo "         OpenMPI. Point EXTRAE_LIB at the openmpi one instead:"
      echo "         find /apps -name libmpitracef.so | grep -i openmpi" ;;
  esac
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

# Locate a writable scratch directory for the traces. The group name from
# id -gn is not the scratch path: it is the project group that owns it.
if [ "$PROFILER" = "extrae" ] && [ -z "$TRACE_DIR" ]; then
  # The user subdirectory under the group scratch may not exist yet: the group
  # directory is created for the project, the per-user one is not.
  for c in "${SCRATCH:-}" /gpfs/scratch/*/"$USER" /gpfs/scratch/*; do
    [ -n "$c" ] || continue
    [ -d "$c" ] || continue
    case "$c" in */"$USER") d="$c" ;; *) d="$c/$USER" ;; esac
    mkdir -p "$d/traces" 2>/dev/null && [ -w "$d/traces" ] && { TRACE_DIR="$d/traces"; break; }
  done
  if [ -n "$TRACE_DIR" ]; then
    echo "=== traces to $TRACE_DIR"
  else
    echo "WARNING: no writable scratch found, traces go under results/ in GPFS home."
    echo "         Find yours with: ls -d /gpfs/scratch/*/\$USER   then set TRACE_DIR."
  fi
fi
[ -n "$TRACE_DIR" ] && mkdir -p "$TRACE_DIR"

mkdir -p results inputs
NP=$(( SLURM_NNODES * SLURM_NTASKS_PER_NODE ))

# Weak scaling: hold unknowns per rank constant, so IDIM follows the rank count.
if [ -n "$WEAK" ]; then
  IDIM=$(awk -v w="$WEAK" -v p="$NP" 'BEGIN{printf "%d", int((w*p)^(1/3) + 0.5)}')
  echo "=== weak scaling: $WEAK unknowns/rank x $NP ranks -> IDIM $IDIM ($(( IDIM**3 )) unknowns)"
fi

TAG="n${SLURM_NNODES}_p${NP}_${MODE}"
[ -n "$IDIM" ]   && TAG="${TAG}_d${IDIM}"
[ -n "$SWEEPS" ] && TAG="${TAG}_j${SWEEPS}"
[ -n "$ITMAX" ]  && TAG="${TAG}_m${ITMAX}"
# SSTEP_SYNC is read inside the solver, not here, but it changes what the run
# measures, so it belongs in the name: otherwise a run with barriers silently
# overwrites the one without, which is exactly its own baseline.
[ -n "$STEPS" ]  && TAG="${TAG}_s${STEPS}"
[ -n "${SSTEP_SYNC:-}" ] && TAG="${TAG}_sync"
# The halo/compute split of the sparse product. It changes nothing about the
# arithmetic, but it does add timing calls inside every matvec, so a timed run
# is not the run to quote a solve time from: keep it under its own name.
[ -n "${PSB_SPMM_TIMERS:-}" ] && TAG="${TAG}_spmm"
[ -n "${PSB_SPMM_BARRIER:-}" ] && TAG="${TAG}_bar"

case "$MODE" in
  converge) SUF="" ;;
  fixed)    SUF="_fixed" ;;
  *) echo "invalid MODE: $MODE"; exit 1 ;;
esac

run_one () {                           # $1 = cg | sstep, $2 = repetition index
  local m=$1 rep=${2:-}
  local out="results/${m}_${TAG}${rep:+_r$rep}"
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
    # Extrae needs an ABSOLUTE path here: given a relative one it ignores it and
    # drops set-0/ in the working directory instead.
    if [ -n "$TRACE_DIR" ]; then
      EXTRAE_DIR="$TRACE_DIR/${out##*/}.extrae"
    else
      EXTRAE_DIR="$PWD/${out}.extrae"
    fi
    mkdir -p "$EXTRAE_DIR"
    EXTRAE_DIR=$(cd "$EXTRAE_DIR" && pwd)
    export EXTRAE_DIR EXTRAE_FINAL_DIR="$EXTRAE_DIR"
    # Fortran MPI interposition. No instrumented binary required.
    export LD_PRELOAD="$EXTRAE_LIB"
  fi

  # An IDIM override writes a copy of the input rather than editing it in place.
  local inp="pde3d_${m}${SUF}.inp"
  if [ -n "$IDIM" ] || [ -n "$SWEEPS" ] || [ -n "$ITMAX" ] || [ -n "$STEPS" ]; then
    local gen="inputs/${m}${SUF}${IDIM:+_d$IDIM}${SWEEPS:+_j$SWEEPS}${ITMAX:+_m$ITMAX}${STEPS:+_s$STEPS}.inp"
    # ITMAX counts outer iterations, each worth s matvecs for the s-step method,
    # so divide there to keep the two runs comparable.
    local sval; sval=${STEPS:-$(awk 'NR==11{print $1}' "$inp")}
    local cap="$ITMAX"
    if [ -n "$ITMAX" ] && [ "$m" = "sstep" ]; then
      cap=$(( ITMAX / sval )); [ "$cap" -lt 1 ] && cap=1
    fi
    awk -v d="$IDIM" -v j="$SWEEPS" -v c="$cap" -v p="$STEPS" '
      NR==3  && d != "" { sub(/^[0-9]+/, d) }
      NR==8  && c != "" { sub(/^[0-9]+/, c) }
      NR==11 && p != "" { sub(/^[0-9]+/, p) }
      NR==21 && j != "" { sub(/^[0-9]+/, j) }
      { print }' "$inp" > "$gen"
    inp="$gen"
    echo "    input: $inp   (itmax $(awk 'NR==8{print $1}' "$gen"), s $(awk 'NR==11{print $1}' "$gen"), sweeps $(awk 'NR==21{print $1}' "$gen"))"
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

# Alternate the two methods rather than looping each in turn, so a slow patch of
# the network hits both equally instead of biasing one.
if [ "$REPS" -le 1 ]; then
  run_one cg
  run_one sstep
else
  for i in $(seq 1 "$REPS"); do
    echo
    echo "########## repetition $i of $REPS"
    run_one cg    "$i"
    run_one sstep "$i"
  done
fi

echo
echo "=== DONE ==="
