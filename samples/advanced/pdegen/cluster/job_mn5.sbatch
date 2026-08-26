#!/bin/bash
#SBATCH --job-name=sstepcg
#SBATCH --output=sstepcg_%j.out
#SBATCH --error=sstepcg_%j.err
#SBATCH --ntasks-per-node=112          # MN5-GPP: 2 x Xeon 8480+ = 112 cores/node
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --qos=gp_debug
#SBATCH --account=ehpc859
#
# Node count is given at submit time:   sbatch --nodes=4 job_mn5.sbatch
#
#   MODE=converge   run to convergence, gives the time to solution (default)
#   MODE=fixed      equal budget of 200 matvecs, gives the cost per iteration
#   PROFILER=extrae produce a Paraver trace (no rebuild needed)
#   PROFILER=none   plain timing run (default)

set -e

BIN=../runs/amg_d_pde3d                # ONE method per run, chosen by KMETHD in the .inp
MODE=${MODE:-converge}
PROFILER=${PROFILER:-none}

module purge
module load bsc/1.0
module load gcc/12.3.0
module load ucx/1.16.0-gcc
module load openmpi/5.0.5-gcc
module load openblas/0.3.27-gcc
[ "$PROFILER" = "extrae" ] && module load extrae
module list

mkdir -p results
NP=$(( SLURM_NNODES * SLURM_NTASKS_PER_NODE ))
TAG="n${SLURM_NNODES}_p${NP}_${MODE}"

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

  srun --ntasks=$NP $BIN "pde3d_${m}${SUF}.inp" > "${out}.txt" 2>&1
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
