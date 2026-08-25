# Misure s-step CG su MN5

Binario: `../amg_d_pde3d` — UN metodo per run, scelto da KMETHD nel .inp.
(Il batch `amg_d_pde3d_sstepbatch` esegue 7 solver di seguito: va bene per
verificare la correttezza, NON per profilare, perché mescola i profili.)

## Input

| file | KMETHD | ITMAX | scopo |
|---|---|---|---|
| pde3d_cg.inp          | CG      | 2000 | tempo a soluzione |
| pde3d_sstep.inp       | SSTEPCG | 2000 | tempo a soluzione |
| pde3d_cg_fixed.inp    | CG      | 200  | budget fisso 200 matvec |
| pde3d_sstep_fixed.inp | SSTEPCG | 40   | budget fisso 40x5 = 200 matvec |

Tutti a IDIM=400 (64M incognite), L1-JACOBI, s=5, BASE_TYPE=M.
ATTENZIONE: s=5 con base monomiale e' l'unica configurazione oggi
verificata corretta. Chebyshev e s>=10 vanno in breakdown.

ITMAX per s-step conta le iterazioni ESTERNE: 40 esterne x s=5 = 200 matvec.

## Lancio

    mkdir -p logs results
    for N in 1 2 4 8 16 32; do sbatch --nodes=$N job_mn5.sbatch; done

    MODE=fixed     sbatch --nodes=8 job_mn5.sbatch     # costo per iterazione
    PROFILER=extrae sbatch --nodes=8 job_mn5.sbatch    # profilo MPI

Compila prima i campi --account e --qos in job_mn5.sbatch
(controlla con `bsc_acct` e `bsc_queues`).
