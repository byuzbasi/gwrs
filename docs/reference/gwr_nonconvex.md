# Formula interfaces for nonconvex GWR

Formula interfaces for nonconvex GWR

## Usage

``` r
gwr_nonconvex(
  formula,
  data,
  coords,
  ...,
  penalty = c("scad", "mcp"),
  gamma = NULL
)

gwr_scad(formula, data, coords, ..., gamma = 3.7)

gwr_mcp(formula, data, coords, ..., gamma = 3)
```

## Arguments

- formula:

  Model formula with an intercept.

- data:

  Data frame containing model variables.

- coords:

  Coordinate matrix or two column names in \`data\`.

- ...:

  Additional arguments passed to \[gwr_nonconvex_fit()\].

- penalty:

  Nonconvex penalty, either \`"scad"\` or \`"mcp"\`.

- gamma:

  Concavity parameter. The default is \`3.7\` for SCAD and \`3\` for
  MCP. SCAD requires \`gamma \> 2\`; MCP requires \`gamma \> 1\`.

## Value

An object inheriting from \`gwr_nonconvex_fit\` and \`gwrs_fit\`.

## Examples

``` r
quake <- datasets::quakes[1:32, ]
fit <- gwr_scad(
  mag ~ depth + stations, quake, coords = c("long", "lat"),
  k = 20, lambda = 0.05,
  control = gwrs_control(n_threads = 1, diagnostics = "none")
)
fit
```
