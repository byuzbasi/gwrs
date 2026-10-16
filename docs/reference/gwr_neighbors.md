# Build a reusable nearest-neighbor structure

\`gwr_neighbors()\` computes nearest neighbors once and stores them in
the target-contiguous layout used by the C++ fitting kernels. Reusing
this object avoids repeating spatial searches across models and tuning
parameters.

## Usage

``` r
gwr_neighbors(
  coords,
  query_coords = NULL,
  k,
  include_self = is.null(query_coords),
  algorithm = c("kd_tree", "cover_tree", "CR", "brute")
)
```

## Arguments

- coords:

  Numeric coordinate matrix with observations in rows, or an \`sf\`
  point object. Optional matrix \`"crs"\` metadata is retained.

- query_coords:

  Optional matrix or \`sf\` point coordinates at new target locations.
  When supplied, neighbors are selected from \`coords\` for every row of
  \`query_coords\`.

- k:

  Number of retained neighbors, including the target itself when
  \`include_self = TRUE\`.

- include_self:

  Whether each observation is inserted as its own first neighbor.

- algorithm:

  Nearest-neighbor algorithm passed to \[FNN::get.knn()\].

## Value

An object of class \`gwrs_neighbors\`.

## Examples

``` r
quake <- datasets::quakes[1:30, ]
coords <- as.matrix(quake[, c("long", "lat")])

neighbors <- gwr_neighbors(coords, k = 18)
neighbors
head(neighbors$index[, 1:3])

query <- coords[1:3, , drop = FALSE] + 0.01
cross_neighbors <- gwr_neighbors(
  coords, query_coords = query, k = 12, include_self = FALSE
)
dim(cross_neighbors$index)
```
