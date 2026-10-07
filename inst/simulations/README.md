# Resumable simulation workflow

This directory contains a local-first simulation runner for future `gwrs`
studies. It is separate from the statistical estimators: R manages immutable
configuration, deterministic task seeds, checkpoints, status files, and
manifests, while model fitting continues to use the package's compiled cores.

The supplied study is a software-workflow example, not a scientific Monte
Carlo design or a performance claim. Its defaults are intentionally moderate.

| Mode | n | p | k | Replications |
|:--|--:|--:|--:|--:|
| `smoke` | 300 | 6 | 40 | 3 |
| `standard` | 2,500 | 12 | 60 | 20 |

## Paper-level penalized selection study

`penalized-selection-study.R` is the separately versioned scientific runner.
Unlike the fixed-anchor validation below, it compares GWR, GWR-Ridge,
GWR-Lasso, GWR-Elastic Net, GWR-SCAD, GWR-MCP, and the simulation-only
Oracle-union and Oracle-local GWR references after spatial cross-validation.
One checkpoint represents one generated dataset and contains all eight
methods, so data, spatial folds, the selected common bandwidth, same-location
test response, and fixed-grid target data are exactly paired.

The frozen heterogeneous `penalized-selection-v2` production design contains
39 cells, 1,900 dataset shards, and 15,200 method fits:

| Section | Cells | Replications | Dataset shards | Purpose |
|:--|--:|--:|--:|:--|
| main | 32 | 50 | 1,600 | two heterogeneous supports, n=2k/10k, rx=0/.7, rs=0/.7, SNR=1/3 |
| reference | 2 | 100 | 200 | stationary sparse and null-slope references, rx=rs=.7 |
| stress | 5 | 20 | 100 | two p=100, local-collinearity, clustered-sampling, and L-domain cells |
| total | 39 | -- | 1,900 | 15,200 method fits |

The main scenarios are `smooth_multiscale_local_sparse`, with broad and
compact plateau-Wendland coefficient surfaces, and
`regional_sparse_discontinuous`, with abrupt regional supports. The two
reference scenarios are `stationary_sparse` and `null_slopes`. All scenarios
retain the common spatially varying intercept, so the latter name means null
slopes rather than a constant-response model. The five stress cells use
`n=10,000`, SNR 3, and `rs=.7`: smooth and regional `p=100` designs; a `p=30`
local-collinearity design with target correlation .9 inside its declared
region; clustered sampling; and an L-shaped-domain boundary design.

For predictors, `rx` is the AR(1) cross-predictor correlation and `rs` is the
spatial autoregressive parameter. After the columnwise AR(1) construction, a
sparse row-standardized 8-nearest-neighbor graph defines
`(I - rs W) X = innovation`; no dense spatial weight matrix is formed. The
spatial predictor process is constructed jointly over the sampled training
locations and fixed target grid, and target predictors use only the centers
and scales learned from the training rows. The design records realized global,
local-region, outside-region, and spatial-neighbor correlation checks.

Coefficient support is stored as an explicit location-by-predictor logical
mask rather than recovered from a floating-point coefficient tolerance.
Separate explicit masks identify support-boundary bands. Independent train,
same-location test, and fixed-grid target errors use
`sigma = sd(training signal) / sqrt(SNR)`. The target design is a 32-by-32
midpoint grid, restricted to the observed L-shaped domain when applicable.
Every estimator receives the same spatial fold assignment. The common `k` is
selected by GWR spatial-CV RMSE from:

```text
max(2(p+1), ceiling(.025 n))
max(4(p+1), ceiling(.050 n))
max(8(p+1), ceiling(.100 n))
```

Penalized models then tune lambda on the same folds. Ridge uses 30 log-spaced
values from `1e2` to `1e-4`; Lasso, Elastic Net, SCAD, and MCP use 30-value
penalty-specific paths. Elastic Net fixes alpha at 0.5, SCAD gamma at 3.7, and
MCP gamma at 3. The Oracle-union model knows the union of predictors active
anywhere. Oracle-local knows the true support separately at every target and
uses grouped cross-location GWR fits. Both are infeasible simulation
benchmarks, not candidate methods.

The primary prediction, estimation, and selection metrics are fixed-grid
noise-free signal RMSE, active-coefficient RMSE, and MCC. Summaries also retain
same-location response/signal RMSE; fixed-grid response, signal, and
total/active/inactive coefficient RMSE; intercept and inactive-leakage RMSE;
active attenuation bias; support-boundary versus core-active coefficient RMSE;
and domain-edge versus interior coefficient and signal RMSE. Sparse estimators
report TP/FP/FN/TN, TPR, FPR, specificity, precision, false discovery rate,
F1, MCC, support IoU, selection accuracy, and mean selected predictors. GWR,
Ridge, Oracle-union, and Oracle-local have `NA` selection metrics. A constant
prediction against a truth containing both classes has MCC zero, so a complete
support miss stays in the Monte Carlo summary. MCC is undefined (`NA`) only
for a one-class truth, including `null_slopes`; FPR, specificity, and false
discoveries remain directly auditable as applicable. Nominal `design_snr` is
recorded separately from realized train, same-location test, and target SNR.
Convergence and timing, including fixed-grid prediction, are retained without
silently deleting failed or nonconverged fits.

A local smoke test is:

```sh
GWRS_SELECTION_SIM_OUTPUT=/absolute/path/penalized-selection-smoke-v2 \
GWRS_SELECTION_SIM_MODE=smoke \
GWRS_SELECTION_SIM_THREADS=2 \
GWRS_SELECTION_SIM_SEED=20260826 \
Rscript --vanilla inst/simulations/penalized-selection-study.R
```

The exact same command resumes. `smoke` has two tiny dataset shards so setting
`GWRS_SELECTION_SIM_MAX_TASKS=1` on the first call proves an interrupted
resume. Production is intentionally not shown as a local command; use the
selection-specific TRUBA handoff documented in `truba/README.md` after the
debug, worker-scaling pilot, and separate max-cell resource calibration pass.

The v2 smoke task IDs are
`smoke-smooth_multiscale_local_sparse-baseline-n00120-p006-rx020-rs030-s30-rep-00001`
and
`smoke-regional_sparse_discontinuous-baseline-n00120-p006-rx020-rs030-s30-rep-00001`.
The current TRUBA example identifies the debug run as
`penalized-selection-v2-debug-v0.4.0-seed20260826-wf4`.

`max_cell_calibration` is an operational mode, not a scientific design cell.
It reproduces the exact first `regional_sparse_discontinuous`, `n=10,000`,
`p=100`, `rx=rs=.7` production stress shard with five folds, the 32-by-32
target grid, and the full 30-value paths. Its exact task ID is
`max-cell-calibration-stress-regional_sparse_discontinuous-p100-n10000-p100-rx070-rs070-s30-rep-00001`
and the current TRUBA example run ID is
`penalized-selection-v2-max-cell-calibration-v0.4.0-seed20260826-wf4`.
Bandwidth CV is retained for audit, while all eight resource fits force the
maximum candidate `k=1,000`. Its output has a distinct scope and must not be
combined with production replications.

After a study has a valid `COMPLETED` marker and all source manifests pass,
create a new, versioned analysis directory with:

```sh
Rscript --vanilla \
  inst/simulations/summarize-penalized-selection-study.R \
  /absolute/path/completed-study \
  /absolute/path/analysis-v2
```

The paper summarizer accepts only the exact 1,900-shard, 15,200-fit production
contract; smoke, pilot, and calibration outputs cannot produce paper tables.
It refuses partial or altered inputs and an existing output directory. It
writes the primary all-completed analysis and a pre-specified all-converged
sensitivity analysis (paired sensitivity rows require both fits to converge).
Paired error contrasts include every non-GWR method versus GWR and direct
SCAD/Lasso and MCP/Lasso comparisons. Its `favors` label names a method only
when the paired 95% interval excludes zero. Outputs also include convergence,
bandwidth, and design tables, a session record, SHA-256 manifest, and
completion marker.

## Joint GWR/Lasso/SCAD/MCP recovery study

`nonconvex-study.R` provides a separate, tall-format comparison workflow.
Every scenario-replication dataset is regenerated from the same recorded
`data_seed` for all requested methods, so GWR, GWR-Lasso, GWR-SCAD, and
GWR-MCP are evaluated on identical observations. Predictor standardization,
kernel, neighbor count, response noise, predictor correlation, solver
tolerance, selection threshold, and penalty concavity parameters are all
recorded in the immutable configuration.

The current protocol evaluates each penalized method at an explicit fraction
of its own all-zero anchor. This makes it useful for software and recovery
validation, but it is not yet the final Monte Carlo design for a paper. A
paper-level comparison should separately approve the tuning design (for
example nested spatial cross-validation), scenarios, estimands, and reporting
rules before production execution.

The smoke run is:

```sh
R_LIBS=/absolute/path/to/project-library \
GWRS_NONCONVEX_SIM_OUTPUT=/absolute/path/nonconvex-study-smoke-v0.4.0 \
GWRS_NONCONVEX_SIM_MODE=smoke \
Rscript inst/simulations/nonconvex-study.R
```

The same command resumes. To inspect progress without fitting:

```sh
GWRS_NONCONVEX_SIM_ACTION=status \
GWRS_NONCONVEX_SIM_OUTPUT=/absolute/path/nonconvex-study-smoke-v0.4.0 \
Rscript inst/simulations/nonconvex-study.R
```

The production defaults are intentionally moderate for a cluster: two
scenarios, four methods, 50 replications, `n = 10,000`, `p = 30`, and
`k = 100` (400 checkpointed tasks). A directly resumable TRUBA template is:

```sh
RCPP_PARALLEL_NUM_THREADS=10 \
R_LIBS=/absolute/path/to/project-library \
GWRS_NONCONVEX_SIM_OUTPUT=/absolute/path/nonconvex-study-production-v0.4.0-n10000-p30-k100 \
GWRS_NONCONVEX_SIM_MODE=production \
GWRS_NONCONVEX_SIM_CONFIRM_LARGE=YES \
GWRS_NONCONVEX_SIM_THREADS=10 \
GWRS_NONCONVEX_SIM_N=10000 \
GWRS_NONCONVEX_SIM_P=30 \
GWRS_NONCONVEX_SIM_K=100 \
GWRS_NONCONVEX_SIM_REPLICATIONS=50 \
Rscript inst/simulations/nonconvex-study.R
```

With `GWRS_NONCONVEX_SIM_KEEP_TASK_DATA=NO` (the default), checkpoints retain
compact one-row method summaries rather than fitted objects or simulated
matrices. `task-summary.csv` reports response, signal, intercept, active and
inactive coefficient RMSE; auditable TP/FP/FN/TN counts; TPR, FPR, precision,
F1, Matthews correlation coefficient (`mcc`), support IoU, and selection
accuracy; solver audits; anchor and fit timing; and memory. Selection metrics
are evaluated only for estimators that produce exact zeros. They are `NA` for
unpenalized GWR (and must likewise remain `NA` for Ridge when it is added to a
paper-level comparison). A constant prediction with a two-class truth has
MCC zero; MCC is `NA` for a one-class truth such as the all-null case, where
FPR remains available. Exact runtime depends on the node and convergence
behavior. Reissue the identical command to resume and use a new output
directory whenever any recorded setting changes.

A guard refuses `n > 20,000`, `p > 100`, `k > 200`, or more than 500 tasks
unless `GWRS_SIM_CONFIRM_LARGE=YES` is supplied explicitly. Full fitted objects
are not retained by default; each task stores compact recovery, selection,
convergence, timing, and memory summaries.

## First run and resume

Use a new versioned output directory. Reissuing the identical command is the
resume command.

```sh
R_LIBS=/absolute/path/to/project-library \
GWRS_SIM_OUTPUT=/absolute/path/simulation-smoke-v1 \
GWRS_SIM_MODE=smoke \
Rscript inst/simulations/template-study.R
```

For a controlled interruption test, limit the first invocation and then omit
the limit on the second invocation:

```sh
R_LIBS=/absolute/path/to/project-library \
GWRS_SIM_OUTPUT=/absolute/path/simulation-smoke-v1 \
GWRS_SIM_MODE=smoke \
GWRS_SIM_MAX_TASKS=1 \
Rscript inst/simulations/template-study.R
```

The exact first command, without `GWRS_SIM_MAX_TASKS`, completes the remaining
tasks. Completed task files are validated and skipped rather than overwritten.

## Status without fitting

```sh
GWRS_SIM_ACTION=status \
GWRS_SIM_OUTPUT=/absolute/path/simulation-smoke-v1 \
Rscript inst/simulations/template-study.R
```

`progress.csv` gives the current phase, active task, counts, percentage, and a
simple remaining-time estimate. `task-status.csv` gives task-level status,
attempt count, elapsed time, checkpoint name, and SHA-256. Status mode is
read-only and can be used while another process owns the run lock.

## Output contract

```text
simulation-smoke-v1/
  run-config.rds             immutable exact configuration
  run-config.csv             human-readable configuration
  runtime-signature.csv      package, DLL, scripts, index width, R, RNG
  session-info.txt           initial R session, libraries, thread environment
  scientific-signature.csv   size/SHA256 signature of immutable run inputs
  preflight.csv              dimensions and rough input allocations
  task-plan.csv              immutable task IDs and seeds
  progress.csv               replaceable, reconstructible run summary
  task-status.csv            replaceable, reconstructible task ledger
  tasks/                     immutable completed checkpoints
  errors/                    immutable error records by attempt
  quarantine/                preserved invalid checkpoints or recovered locks
  task-summary.csv            compact combined scientific summaries
  task-manifest.csv           task seeds, sizes, SHA256, timing, completion UTC
  manifest-sha256.csv         immutable output size/SHA256 manifest
  COMPLETED                   written only after every final validation passes
```

The configuration, task plan, runtime signature, task checkpoints, error
records, and final manifest are never silently overwritten. Mutable progress
files contain only derived tracking information and can be reconstructed from
the immutable files.

## Failure and lock handling

Failed tasks retain an error record and are retried on the next identical run
by default. Set `GWRS_SIM_RETRY_FAILED=NO` to inspect without retrying. A run
lock prevents two writers from using the same output directory. After a
confirmed process crash, `GWRS_SIM_RECOVER_LOCK=YES` archives the old lock in
`quarantine/`; it should not be used while a process is still active.

An invalid or unreadable checkpoint is moved to `quarantine/` and recomputed.
A changed configuration, task plan, workflow script, installed shared library,
package version, platform, or RNG kind is refused. A finalized run is accepted
only when every file listed in `manifest-sha256.csv` still matches its size and
checksum and `COMPLETED` binds that manifest to the scientific signature.

The complete, single-node TRUBA/ARF installation, debug, pilot, production,
and exact-resume handoff is documented in the project-level `truba/README.md`.

## Scope of resume

The checkpoint unit is one scenario-replication task. Interruption can lose at
most the currently running task. Continuing from an internal C++ optimization
iteration would require a separate solver-state API and is deliberately not
part of this workflow layer.
