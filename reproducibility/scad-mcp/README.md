# SCAD/MCP paper: independent reproduction

Reproduction code for **Nonconvex Geographically Weighted Regression for Local
Variable Selection**, Bahadir Yuzbasi. The package and this bundle are supplied
under GPL (>= 3). Public input data are credited to the U.S. Census Bureau.

## Requirements

- macOS or Linux; Python >= 3.9; R >= 4.1; a C++17 compiler and GNU make.
- R packages listed in `setup.R` and each source archive's DESCRIPTION.
  The real-data workflow requires `sf` and its GDAL/GEOS/PROJ dependencies;
  rebuilding geometry also requires `spdep`. Binary R packages may be used
  where available. Installation of system libraries is platform-specific.
- R and Python must be on PATH. Paths containing spaces are supported.
- Write access to the extracted directory. Use a new copy for a different
  compiler, R version or dependency stack. Do not share a writable runtime
  library between simultaneous setup operations.

No operating-system packages are installed by the launcher. By default setup
only uses already installed R dependencies and builds the included `gwrs`
source into `runtime/`. The optional `--install-dependencies` flag explicitly
allows downloading missing declared R dependencies into `runtime/dependencies`;
it does not replace packages in a user's library. The included source archives
must be used for reproducing the paper, even if a newer `gwrs` is installed.

All commands below start from the extracted `scad-mcp` directory. Replace
`/path/to/scad-mcp` with that directory's absolute path once. The numerical
settings and seeds are read from the included code/configuration. No fitting
occurs during setup, check, preflight or data reconstruction.

## Setup and integrity checks

```sh
cd /path/to/scad-mcp
python3 reproduce.py simulation setup
python3 reproduce.py real-data setup
python3 reproduce.py simulation check
python3 reproduce.py real-data check
python3 reproduce.py simulation preflight --mode full --run-id paper-v1
python3 reproduce.py real-data preflight --mode full --run-id paper-v1
```

If setup reports missing R dependencies, install the declared packages yourself
or rerun that setup command with `--install-dependencies`. Setup verifies the
loaded package version, project library and source-archive receipt. Each run
records actual R, platform and dependency versions and compiled-library hashes.
A changed runtime or scientific signature is refused on resume.

## Small lifecycle tests

These are executable checks, not paper results. The simulation smoke contains
two small datasets (n=120, p=6); the application smoke uses synthetic n=40, p=5
data and 210 small model-path tasks. It does not refit the 3,107-county study.

```sh
cd /path/to/scad-mcp
python3 reproduce.py simulation smoke --run-id smoke-v1 --max-new-tasks 1
python3 reproduce.py simulation resume --run-id smoke-v1
python3 reproduce.py simulation verify --run-id smoke-v1
python3 reproduce.py real-data smoke --run-id smoke-v1 --max-new-tasks 4 --workers 2
python3 reproduce.py real-data resume --run-id smoke-v1 --workers 2
python3 reproduce.py real-data verify --run-id smoke-v1
```

Use the same mode and simulation thread count when resuming. Run/smoke refuses
an existing output directory; resume validates existing shards and computes
only missing work. Invalid existing shards cause an error; retain their logs
and investigate instead of deleting completion markers. Do not change a
library, source, grid or input during a run.

## Full paper simulation (run explicitly)

The experiment has **39 design cells, 1,900 dataset replications and 15,200
method fits**, including two simulation-only oracle benchmarks. The exact
configuration and dataset-level seeds are in `simulation/reference/`.
Dimensions vary across cells; the largest cell has n=10,000 and p=100.
The study uses GWR, Ridge, Lasso, Elastic Net, SCAD and MCP. The two oracle
fits use known simulated support and are not empirical candidate methods.

```sh
cd /path/to/scad-mcp
python3 reproduce.py simulation run --mode full --run-id paper-v1 --threads 1 --confirm-full
```

Resume, progress, completion verification and reporting:

```sh
cd /path/to/scad-mcp
python3 reproduce.py simulation status --mode full --run-id paper-v1
python3 reproduce.py simulation resume --mode full --run-id paper-v1 --threads 1 --confirm-full
python3 reproduce.py simulation verify --mode full --run-id paper-v1 --threads 1
python3 reproduce.py simulation report --mode full --run-id paper-v1 --threads 1
```

Output: `simulation/results/paper-v1/`. The original execution used 56 native
threads. `--threads` is a resource setting; choose it before starting, and keep
it unchanged on resume. The default of one avoids assuming a particular
machine. BLAS remains single-threaded. There is no dataset-level parallelism
in this launcher. The paper summarizer retains valid fits with convergence
warnings in the primary analysis and reports a separate converged sensitivity
analysis; a completed run is not a claim that every fit converged.

## Full county application (run explicitly)

The complete application has **3,107 counties, 26 predictors and 2,272 tasks**:
1,760 nested inner tuning paths, 60 outer refits, 440 full-sample tuning paths
and 12 final refits. Local and global versions of GWR/OLS, Ridge, Lasso,
Elastic Net, SCAD and MCP are compared using the saved spatial folds and
training-only preprocessing. Exact grids and numerical gates are in
`real-data/study-v3/config/numerics.json` and `policy.json`.

```sh
cd /path/to/scad-mcp
python3 reproduce.py real-data run --mode full --run-id paper-v1 --workers 1 --confirm-full
```

Resume, progress, completion verification and reporting:

```sh
cd /path/to/scad-mcp
python3 reproduce.py real-data status --mode full --run-id paper-v1
python3 reproduce.py real-data resume --mode full --run-id paper-v1 --workers 1 --confirm-full
python3 reproduce.py real-data verify --mode full --run-id paper-v1
python3 reproduce.py real-data report --mode full --run-id paper-v1
```

Output: `real-data/study-v3/results/paper-v1/`, with publication figures and
tables in sibling `paper-v1-publication/` and `paper-v1-distance/` directories.
Workers run independent R processes, each with one native/BLAS thread. Increase
`--workers` only to match available RAM and CPUs; no nested thread pool is used.
The three internal `study-v*` directories are self-contained implementation
layers of this one application, not separate studies to run manually.

Reporting reads completed validated outputs only. Simulation reporting refuses
an existing report directory; a valid existing summary may be reused. ACS
reporting verifies an existing completed report. Neither command refits models.

## Input reconstruction and dataset use

The included public Census snapshot is sufficient to rebuild the input table,
projected geometry, spatial folds and common diagnostic graph offline:

```sh
cd /path/to/scad-mcp
python3 reproduce.py data setup
python3 reproduce.py data prepare --run-id rebuild-v1
```

The derived copy goes to `work/rebuild-v1/`. Repeating the command validates and
reuses valid existing stages. Original input files are never overwritten.
The final check compares every model-data field with the supplied table
(numeric tolerance 1e-12), and requires identical saved fold and graph tables.
External geometry-library differences can fail this strict reconstruction
check; preserve the evidence and use the supplied verified analysis inputs
rather than silently changing the paper's sample or folds.

For optional source acquisition, `data download --run-id download-v1` creates a
separate working copy and verifies included downloads. The original snapshot
is retained. This action requires network access if any source file is absent.
It is not a request to refresh the analysis to a different ACS vintage.

The current `gwrs` package additionally supplies the same 3,107 x 91 table as
`data(acs2024_counties)`. See `?acs2024_counties` for all fields, exclusions,
units, MOE handling and citation guidance. The historical paper libraries
predate this dataset export; they read the included CSV directly.

## Monitoring, resources and reproducibility

Progress is written to stdout, `progress.json` and `progress.tsv` every
30 seconds and at stage transitions (`--heartbeat` changes the interval).
Completion percentages use validated tasks. ETA uses measured current-run
costs only after pending task classes are represented; otherwise it is unknown.
ETA covers the current computational stage, excluding later validation and
optional rendering. Any SLURM wall-time value is an allocation deadline,
not a runtime prediction. The status command flags stale running heartbeats.

Full-run elapsed time is not estimated from the tiny smoke tests. Runtime
and memory depend strongly on hardware and candidates. For ACS, the inherited
preflight requires at least 6 GiB free and recommends 12 GiB. Simulation storage
and peak RAM must be budgeted for the maximum n/p cell; its compiled engine
uses sparse neighborhoods. Reserve additional space for rendered figures.

Checkpoints bind the study configuration, runtime and source identity. A
completion marker alone is insufficient: use `verify`, which checks task
counts, dimensions, numerical validity and SHA-256 manifests. Never mix a
partial run from another machine/runtime with the same output directory.
Full production fitting is not run automatically by setup or tests.

See [SOURCES.md](SOURCES.md) and [OUTPUTS.md](OUTPUTS.md). The [validation report](VALIDATION.md) distinguishes executed smoke checks from unexecuted full
production runs. Bitwise equality across different R/compiler/BLAS/GDAL
versions is not promised; record and compare the numerical diagnostics.
