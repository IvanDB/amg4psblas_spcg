#!/bin/bash
#SBATCH --job-name=sstepcg
#SBATCH --output=logs/%x_%j_n%N.out
#SBATCH --error=logs/%x_%j_n%N.err
#SBATCH --ntasks-per-node=112          # MN5-GPP: 2 x Xeon 8480+ = 112 core/nodo
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --account=CHANGE_ME            # verifica con:  bsc_acct
#SBATCH --qos=CHANGE_ME                # verifica con:  bsc_queues
# I nodi si scelgono al lancio:  sbatch --nodes=4 job_mn5.sbatch
set -euo pipefail

BIN=../amg_d_pde3d                     # UN metodo per run, scelto dal .inp
MODE=${MODE:-converge}                 # converge | fixed
PROFILER=${PROFILER:-none}             # none | scorep | extrae | mpip

module purge
module load intel impi                 # oppure: module load gcc openmpi
module list 2>&1

mkdir -p logs results
NP=$(( SLURM_NNODES * SLURM_NTASKS_PER_NODE ))
TAG="n${SLURM_NNODES}_p${NP}_${MODE}"

case "$MODE" in
  converge) SUF="" ;;                  # a convergenza: misura il tempo a soluzione
  fixed)    SUF="_fixed" ;;            # budget fisso di 200 matvec: costo per iterazione
  *) echo "MODE non valido"; exit 1 ;;
esac

run_one () {                           # $1 = cg | sstep
  local m=$1
  local out="results/${m}_${TAG}"
  echo "=== $m  nodi=$SLURM_NNODES  rank=$NP  mode=$MODE"

  case "$PROFILER" in
    scorep)
      # richiede binari ricompilati con i wrapper scorep (vedi note)
      export SCOREP_EXPERIMENT_DIRECTORY="${out}.scorep"
      export SCOREP_ENABLE_PROFILING=true
      export SCOREP_ENABLE_TRACING=false      # metti true SOLO dopo scorep-score
      export SCOREP_FILTERING_FILE="$SLURM_SUBMIT_DIR/scorep.filter"
      export SCOREP_TOTAL_MEMORY=64M
      export SCOREP_MPI_ENABLE_GROUPS=cg,coll,p2p,xnonblock
      ;;
    extrae)
      export EXTRAE_CONFIG_FILE="$SLURM_SUBMIT_DIR/extrae.xml"
      export LD_PRELOAD=${EXTRAE_HOME}/lib/libmpitracef.so
      ;;
    mpip)
      export LD_PRELOAD=${MPIP_HOME}/lib/libmpiP.so
      export MPIP="-t 5.0 -k 4"
      ;;
  esac

  srun --ntasks=$NP $BIN "pde3d_${m}${SUF}.inp" > "${out}.txt" 2>&1
  unset LD_PRELOAD || true

  grep -E "Krylov method|step size|Iterations to convergence|Time to solve|Time per iteration|Residual 2-norm|Linear system size" \
       "${out}.txt" || true
}

run_one cg
run_one sstep
