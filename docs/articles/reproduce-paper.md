# Reproduce the SCAD/MCP paper

The accompanying manuscript is *Nonconvex Geographically Weighted
Regression for Local Variable Selection*. Its independent reproduction
bundle contains the simulation and ACS application scripts, public
Census input snapshot, historical package source archives,
configuration, seeds and SHA-256 manifest.

## Download the materials

- [Independent reproduction
  ZIP](https://github.com/byuzbasi/gwrs/blob/main/reproducibility/gwrs-scad-mcp-reproduction-v2.zip)
- [Current package source, including ACS
  data](https://github.com/byuzbasi/gwrs/blob/main/reproducibility/gwrs_0.4.0.9013.tar.gz)
- [Installation and run
  guide](https://github.com/byuzbasi/gwrs/blob/main/reproducibility/scad-mcp/README.md)
- [Executed validation and
  limits](https://github.com/byuzbasi/gwrs/blob/main/reproducibility/scad-mcp/VALIDATION.md)
- [Source
  provenance](https://github.com/byuzbasi/gwrs/blob/main/reproducibility/scad-mcp/SOURCES.md)

Repository downloads require repository access. The independently
supplied ZIP can also be used without a GitHub checkout.

## Check before a full run

After extracting the reproduction ZIP, open a terminal in its `scad-mcp`
directory. With the required R/Python toolchains and declared
dependencies available, the launcher provides setup and source checks:

``` sh
python3 reproduce.py simulation setup
python3 reproduce.py real-data setup
python3 reproduce.py simulation check
python3 reproduce.py real-data check
```

Follow the standalone guide for preflight, small smoke tests, full
execution, progress, checkpoint/resume, verification and reporting. A
completion marker or progress percentage alone is insufficient: the
verification command checks task completeness, scientific signatures and
output integrity.

## Software versions and validation

The current package is **0.4.0.9013**. The paper’s simulation used
**0.4.0** and its county analysis used **0.4.0.9006**. The launcher
installs those historical sources into separate project libraries.
Retain them to reproduce the reported analyses.

The distribution validation records a local package check with **Status:
OK** and **1,203 passing testthat expectations**, together with small
interrupted and resumed study workflows. These are recorded checks for
the distributed sources, not a claim of continuous-integration coverage
on every platform. No full study was rerun for the documentation
preview.

## Cite the software and data

``` r

citation("gwrs")
```

The package citation identifies the software and accompanying
manuscript. The manuscript remains an unpublished work until a
publication record exists. Please also acknowledge the Census Bureau
when using the prepared ACS table.
