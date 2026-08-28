# s-step CG measurements on MareNostrum 5

Account `ehpc859`. Available QOS:

| QOS              | max time   | max procs | nodes | note                  |
|------------------|------------|-----------|-------|-----------------------|
| `gp_debug`       | 02:00:00   | 3584      | 32    | set in both scripts   |
| `gp_ehpc`        | 3-00:00:00 | 89600     | 800   | default, production   |
| `gp_interactive` | 02:00:00   | 32        | -     | login nodes only      |

gp_debug reaches 32 nodes, which covers the whole scaling grid below, and
usually schedules sooner. gp_ehpc is the default and is there for anything
longer or wider.

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

Validate the pipeline on a small case first. IDIM=100 is 1M unknowns, a couple
of minutes, and it exercises exactly the same code path as the production size:

    IDIM=100 sbatch --nodes=1 job_mn5.sbatch
    IDIM=100 PROFILER=extrae sbatch --nodes=1 job_mn5.sbatch

Then the real measurements:

    for N in 1 2 4 8 16 32; do sbatch --nodes=$N job_mn5.sbatch; done   # strong scaling
    MODE=fixed       sbatch --nodes=8 job_mn5.sbatch                    # cost per iteration
    PROFILER=extrae  sbatch --nodes=8 job_mn5.sbatch                    # Paraver trace

An IDIM override writes a copy under `inputs/` instead of editing the input in
place, and tags the result files with the size.

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

## 5. Scaling campaign

Strong scaling, fixed 64M unknowns, growing node count. The local problem
shrinks as nodes grow, which is the regime where s-step is supposed to win:

| nodes | ranks | unknowns/rank |
|-------|-------|---------------|
| 1     | 112   | 571k          |
| 2     | 224   | 286k          |
| 4     | 448   | 143k          |
| 8     | 896   | 71k           |
| 16    | 1792  | 36k           |

    for N in 1 2 4 8 16; do sbatch --nodes=$N job_mn5.sbatch; done

gp_debug caps nodes and wall time, so anything past a couple of nodes needs a
production QOS: `sbatch --nodes=16 --qos=<production> --time=01:00:00 ...`

Collect everything into one table when the jobs land:

    ./collect.sh            # markdown
    ./collect.sh --tsv      # for gnuplot or a spreadsheet

The number to plot is time per iteration against node count, CG and SSTEPCG on
the same axes. Where the s-step curve crosses below CG is the result. At one
node it sits about 13% above: all-reduces are intra-node and cost almost
nothing, so the method pays its Gram-matrix arithmetic for no saving.

## 6. Tracing strategy

Do not trace every point of the scaling curve. Traces at 1792 ranks are large
and mostly redundant. Two runs answer the question:

- the smallest node count, where CG wins, and
- the largest one available, where s-step should.

Comparing the two shows what changed, which is more informative than either
one alone.

    PROFILER=extrae sbatch --nodes=1  job_mn5.sbatch
    PROFILER=extrae sbatch --nodes=16 job_mn5.sbatch

If the traces get unwieldy, the knob is `<trace-control>` in extrae.xml, which
restricts tracing to a window instead of the whole run, including the matrix
generation that is of no interest here.

What the trace has to confirm, before anything else is read into it: CG issues
about 3 all-reduces per matrix-vector product, s-step 3 per OUTER iteration,
so 3 per s of them. At s=5 that is a factor of five. If the counts do not show
that ratio, the model is wrong and the timings mean nothing.

## 7. The heavy traced run

Two things decide whether this succeeds: how big the trace gets, and where it
lands.

Trace size grows with ranks times MPI calls. At 1792 ranks a full solve is
around 300 iterations, each with roughly six halo exchanges and three
all-reduces, and every event carries six levels of call stack. That is tens of
millions of events. Nothing about the per-iteration pattern needs 300
iterations to be visible, so cap the run:

    ITMAX=30 counts MATVECS, and is divided by s for the s-step input, so both
    methods do the same 30 and stay comparable.

Traces must not go to GPFS home, which is quota limited and slow for many
small writes:

    TRACE_DIR=<your scratch>/traces

Find your scratch with `bsc_quota`, which lists the filesystems and their
paths. Left unset, the job looks for one itself.

The run worth tracing, given what the scaling curve showed:

    mkdir -p /gpfs/scratch/<group>/$USER/traces
    SWEEPS=1 ITMAX=30 MODE=fixed PROFILER=extrae \
      TRACE_DIR=/gpfs/scratch/<group>/$USER/traces \
      sbatch --nodes=16 --qos=<production> --time=00:30:00 job_mn5.sbatch

SWEEPS=1 because six Jacobi sweeps bury the all-reduces under six halo
exchanges; MODE=fixed because equal work makes the two traces comparable
side by side.

Take the same run at one node as a baseline. The comparison between the two is
what shows how the communication share moved, which neither trace tells alone.

If the merge fails at high rank count, keep-mpits is on in extrae.xml, so it
can be redone afterwards on a login node with mpi2prv.
