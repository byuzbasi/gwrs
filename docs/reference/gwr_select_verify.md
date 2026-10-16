# Verify a joint GWR selection checkpoint

Verify a joint GWR selection checkpoint

## Usage

``` r
gwr_select_verify(checkpoint)
```

## Arguments

- checkpoint:

  Directory created by \[gwr_select()\].

## Value

Verification summary with candidate counts, search status and
provenance. All recorded packs/results are checked for size and SHA-256;
errors stop verification. Uncommitted orphan files are not counted.

## Details

This is read-only and does not fit models. Completion and numerical
eligibility are separate: a complete search can contain failed
candidates. A stale run.lock is never removed automatically. After
confirming its owner is no longer running, archive that lock directory
before resuming.
