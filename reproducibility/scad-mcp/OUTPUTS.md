# Output schema and acceptance

## Simulation

- `study/run-config.rds`, `run-config.csv`, `task-plan.csv`: exact design and seeds.
- `study/runtime-signature.csv`, `scientific-signature.csv`, `session-info.txt`:
  software, compiled engine, scripts and scientific identity.
- `study/tasks/`: atomic per-dataset checkpoints and checksums; all eight
  methods share the same generated observations, folds and evaluation targets.
- `study/task-summary.csv`: method-level prediction, coefficient and selection
  metrics, convergence diagnostics, tuning choices and timings.
- `study/task-status.csv`: validated completed/failed/pending task records.
- `study/manifest-sha256.csv` and `COMPLETED`: final output sizes/checksums.
- `summary-v1/`: cell-level metrics, paired contrasts, Monte Carlo summaries and
  separate all-completed/all-converged analysis sets.
- `report/v1/`: CSV/LaTeX tables and publication figures.

Full completion requires 1,900 valid datasets and 15,200 method rows. Smoke
completion requires two datasets and 16 method rows. Undefined null-support
selection metrics remain missing; convergence warnings are retained and
reported, not silently filtered from the primary estimand.

## ACS application

- `numeric/task-plan.*`, `identity.rds`, `scientific-signature.*` and
  `runtime-*`: complete plan and binding to data/configuration/software.
- `shards/`: checksummed per-path and refit results.
- `numeric/selection-*.rds`: validated tuning selections for outer folds and
  full-data refits.
- `numeric/oof-*.rds`, `numeric/oof-*.csv`: pooled/fold prediction metrics and
  common-graph residual Moran/LISA diagnostics.
- `numeric/full-*.rds`, `numeric/gwr-*.csv`: full-data selections, coefficients
  and unpenalized GWR diagnostics, including F1/F2/F3.
- `manifest-sha256.csv`, `COMPLETED`: final numerical completion record.
- `invocations/`: retained preflight, task, failure and verification logs.

Full verification requires 2,272 tasks and 3,107 out-of-fold predictions per
method, plus the original convergence/stationarity and output consistency
gates. Smoke verification requires its 210 tasks. F-tests are not classical
post-selection tests for penalized estimates. Selection maps use the recorded
coefficient threshold, not a significance decision.

Both studies expose `progress.json` and `progress.tsv` alongside immutable
numerical outputs. These operational records and invocation logs may change
on verification/resume; completed numerical shards must not change. An
interrupted run is incomplete even if some outputs or a stale heartbeat exist.
The verify command, not the progress percentage alone, establishes acceptance.
