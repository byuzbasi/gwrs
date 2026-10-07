# Core benchmark

## SCAD/MCP path benchmark

`benchmark-nonconvex.R` benchmarks the common RcppArmadillo GWR-SCAD/MCP
path core under four controlled cases: serial/parallel crossed with
screening on/off. Each penalty uses one frozen simulated dataset and its own
exact all-zero anchor. The script verifies convergence, stationarity,
coordinatewise objective gaps, and numerical agreement with the serial
unscreened reference before finalizing a run.

The default is a deliberately small protocol smoke test:

```sh
R_LIBS=/absolute/path/to/project-library \
GWRS_NONCONVEX_BENCH_OUTPUT=/absolute/path/nonconvex-benchmark-smoke-v0.4.0 \
GWRS_NONCONVEX_BENCH_MODE=smoke \
Rscript inst/benchmarks/benchmark-nonconvex.R
```

The same command resumes an interrupted run. Read-only status is available
with:

```sh
GWRS_NONCONVEX_BENCH_ACTION=status \
GWRS_NONCONVEX_BENCH_OUTPUT=/absolute/path/nonconvex-benchmark-smoke-v0.4.0 \
Rscript inst/benchmarks/benchmark-nonconvex.R
```

A moderate TRUBA production template is:

```sh
RCPP_PARALLEL_NUM_THREADS=10 \
R_LIBS=/absolute/path/to/project-library \
GWRS_NONCONVEX_BENCH_OUTPUT=/absolute/path/nonconvex-benchmark-production-v0.4.0-n25000-p30-k120 \
GWRS_NONCONVEX_BENCH_MODE=production \
GWRS_NONCONVEX_BENCH_CONFIRM_LARGE=YES \
GWRS_NONCONVEX_BENCH_N=25000 \
GWRS_NONCONVEX_BENCH_P=30 \
GWRS_NONCONVEX_BENCH_K=120 \
GWRS_NONCONVEX_BENCH_N_LAMBDA=25 \
Rscript inst/benchmarks/benchmark-nonconvex.R
```

Runtime is machine-dependent because the benchmark intentionally includes an
unscreened baseline. The final `task-summary.csv` contains timings, speedups,
convergence audits, numerical differences, retained-object sizes, and the
dense spatial allocation avoided. `manifest-sha256.csv` is written only after
both penalty tasks pass. A completed status must report `state = complete`,
`finalized = TRUE`, and `manifest_valid = TRUE`.

Changing any dimension, seed, tolerance, installed binary, or script requires
a new output directory. `GWRS_NONCONVEX_BENCH_MAX_TASKS=1` can be used for a
controlled resume test; `GWRS_NONCONVEX_BENCH_RECOVER_LOCK=YES` is only for a
confirmed stale lock.

