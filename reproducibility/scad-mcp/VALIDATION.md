# Distribution validation — 7 October 2026

## Executed locally

- R 4.6.0, macOS arm64. A relocated directory containing spaces was used,
  with commands launched from a different working directory.
- Both included historical source packages installed into their own project
  libraries. No author-installed `gwrs` library was required.
- Full-design preflight: 1,900 simulation datasets and 2,272 ACS tasks; zero
  full-study models fitted during preflight.
- Simulation lifecycle: two small datasets, interrupted after one, resumed to
  completion, independently verified. Existing shard bytes were preserved.
- ACS lifecycle: synthetic n=40/p=5, interrupted after four tasks, resumed to
  all 210 tasks, independently verified. Existing shard bytes were preserved.
- Completed resume performed no new fitting. Attempts to overwrite existing
  runs, verify a smoke run as full, use a corrupted checkpoint or use modified
  source code were rejected in both appropriate workflows.
- Progress JSON/TSV, completed counts and zero failed tasks were checked.
- Forty R source files and the distributed Python sources parsed successfully.
- The current R package check completed with **Status: OK**: 1,203 testthat
  expectations passed; zero failures, warnings or skipped tests. Examples and
  vignette code passed. PDF manuals were not built.
- ACS data export: 3,107 x 91, all fields and row order identical to the reference
  CSV. Leading-zero geographic identifiers and the accented county name are
  retained, with UTF-8 encoding declared.
- Input reconstruction from the public Census snapshot passed: all data fields,
  saved fold assignments and the common graph matched the supplied references.
- The simulation summary and report were regenerated from existing completed
  results, without fitting: 1,872 primary rows and 468 direct paired rows were
  reconciled, and nine figures were created.
- Ninety-eight provenance/equality checks passed. Current estimator source and
  historical R/native source files are unchanged. The simulation's top-level
  numerical function expressions match the original source; ACS model,
  preprocessing and diagnostic implementations are unchanged.

## Scope and limits

No full simulation or empirical estimation run was repeated for this release.
The complete ACS publication renderer was retained with a subprocess path
quoting correction; it was not rerun against a newly estimated full application.
The small lifecycle test exercises estimation, selection, aggregation and
numerical verification, not the full-size publication maps.

This validation does not certify bitwise equality on untested operating systems,
R versions, compilers or external spatial libraries. Linux users should run
setup, preflight and the documented smoke tests before committing substantial
resources. Exact scientific grids, seeds and thresholds remain unchanged.

`manifest-sha256.json` identifies every distributed workflow file. The numerical
implementation signature excludes Markdown documentation; adding this report
does not change the tested computational sources. The current package archive
has only the local build username removed from DESCRIPTION; all other checked
archive members are unchanged. Personal paths, host identifiers and internal
operator notes are excluded from the publication tree and nested archives.
