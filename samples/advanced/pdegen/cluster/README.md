# s-step CG measurements on MareNostrum 5

Account `ehpc859`, QOS `gp_debug`. Adjust in both scripts if you move to a
production QOS with longer wall time.

## 1. Build

    sbatch build_mn5.sbatch

Compiles PSBLAS, AMG4PSBLAS and the pdegen samples on a compute node, so the
login node never runs a long compile. Can be submitted from anywhere: both
repositories are located automatically, starting from the submission directory
and falling back to a search under `$HOME`. They need not be siblings.

PSBLAS is installed to `$HOME/opt/psblas3-sstep` by default. To pin any of it:

    sbatch --export=ALL,PSBLAS_SRC=/path/to/psblas3_spcg,AMG_SRC=/path/to/amg4psblas_spcg build_mn5.sbatch
    sbatch --export=ALL,PREFIX=$HOME/opt/psblas build_mn5.sbatch

PSBLAS *must* be installed before AMG4PSBLAS is configured, otherwise the
samples link against a stale library and the fixed bugs come back silently.
Add `SKIP_CONFIGURE=1` to reuse an existing `Make.inc` on a rebuild.

No Score-P: Extrae interposes on MPI at run time and needs no instrumented
binary, so this is a plain optimised build.

## 2. Inputs

The driver is `amg_d_pde3d`, which runs a SINGLE method per job, chosen by
KMETHD. The batch driver `amg_d_pde3d_sstepbatch` runs seven solvers back to
back: good for correctness sweeps, useless for profiling because the trace
would blend all of them.

| file                    | KMETHD  | ITMAX | purpose                        |
|-------------------------|---------|-------|--------------------------------|
| `pde3d_cg.inp`          | CG      | 2000  | time to solution               |
| `pde3d_sstep.inp`       | SSTEPCG | 2000  | time to solution               |
| `pde3d_cg_fixed.inp`    | CG      | 200   | fixed budget, 200 matvecs      |
| `pde3d_sstep_fixed.inp` | SSTEPCG | 40    | fixed budget, 40 x s=5 matvecs |

All at IDIM=400 (64M unknowns), L1-JACOBI, s=5, BASE_TYPE=M.

For SSTEPCG, ITMAX counts OUTER iterations, each worth s matrix-vector
products. That is why the fixed-budget pair is 200 against 40: same 200
matvecs for both methods, which isolates the cost per iteration from the
number of iterations needed.

BASE_TYPE=M with s=5 is the only configuration currently verified correct: it
reproduces CG iteration for iteration. The Chebyshev basis needs roughly four
times as many iterations and breaks down for s >= 10, so measuring it would
profile a numerical defect rather than the communication pattern.

## 3. Runs

    for N in 1 2 4 8 16 32; do sbatch --nodes=$N job_mn5.sbatch; done   # strong scaling
    MODE=fixed       sbatch --nodes=8 job_mn5.sbatch                    # cost per iteration
    PROFILER=extrae  sbatch --nodes=8 job_mn5.sbatch                    # Paraver trace

Each job runs CG and SSTEPCG in the same allocation, so the comparison is not
polluted by landing on different nodes. Output goes to `results/`.

## 4. Profiling with Extrae

`extrae.xml` traces MPI plus PAPI counters, with MPI call stacks six levels
deep. That depth is the point: in Paraver it separates an MPI_Allreduce issued
by `psb_dscg` from one issued by the preconditioner or the convergence check,
which is exactly the breakdown needed here.

What to look for: CG issues 3 all-reduces per matvec, s-step v1 issues 3 per
OUTER iteration, so 3 per s matvecs. At s=5 that is five times fewer. The
question the trace has to answer is whether that saving shows up in wall time,
and at how many nodes it starts to.

To profile inside the kernel rather than just at MPI boundaries, two options:

- rebuild with `-finstrument-functions` and set `<user-functions enabled="yes">`
  with a function list. Cheap to try, but the list needs pruning or the trace
  explodes on the inner loops;
- add `Extrae_event` calls at the phase boundaries in `psb_dscg` (matrix power
  kernel, Gram factorisation, Gram solves, block update). More work, far more
  readable in Paraver, and it survives compiler inlining.

Merge happens automatically on exit. Open `trace.prv` in Paraver.
