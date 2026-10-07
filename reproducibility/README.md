# Reproducing the SCAD/MCP paper

The [standalone reproduction guide](scad-mcp/README.md) contains setup,
preflight, smoke, run, resume, verify and report commands for both studies.
It can be used from this checkout or from the separately supplied ZIP.

The source snapshot, model inputs, historical estimation engines and SHA-256
manifest are included. No author home directory, installed package library,
cluster account or existing project checkout is required.

The repository remains private. A repository URL alone does not provide
reviewer access; the separate ZIP can be uploaded as supplementary code.

The package at the repository root is the current implementation. The paper's
simulation and ACS application use the versions recorded in their original
analyses (0.4.0 and 0.4.0.9006 respectively), installed into separate project
libraries by the reproduction launcher. See [source provenance](scad-mcp/SOURCES.md).

## Downloadable files

- [Independent reproduction ZIP](gwrs-scad-mcp-reproduction-v2.zip)
- [Current package source with the ACS dataset](gwrs_0.4.0.9013.tar.gz)
- [SHA-256 checksums](SHA256SUMS) and [file sizes](distribution-manifest.json)
- [Executed validation and limits](scad-mcp/VALIDATION.md)

The ZIP contains both complete study workflows, their required historical
package sources, the Census input snapshot and the current package source.
Compiled libraries, full run outputs and author build logs are not included.
