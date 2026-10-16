# Predict from a fitted geographically weighted model

For Ridge, Lasso, Elastic Net, SCAD and MCP, new target locations are
fitted from retained training data using only their training neighbors.
\`newdata\` supplies predictor values at targets, while \`newcoords\`
controls the local training samples.

## Usage

``` r
# S3 method for class 'gwrs_fit'
predict(
  object,
  newdata = NULL,
  newcoords = NULL,
  neighbors = NULL,
  k = NULL,
  type = c("response", "coefficients", "both"),
  ...
)
```

## Arguments

- object:

  A fitted \`gwrs_fit\`.

- newdata:

  Optional target design matrix or formula data frame.

- newcoords:

  Coordinates corresponding to \`newdata\`.

- neighbors:

  Optional cross-neighbor object constructed with
  \`gwr_neighbors(training_coords, query_coords = newcoords, ...)\`.

- k:

  Number of training neighbors when \`neighbors\` is not supplied.

- type:

  Return predictions, local coefficients, or both.

- ...:

  Unused.

## Value

A numeric prediction vector, coefficient matrix, or a list of both.

## Examples

``` r
train <- datasets::quakes[1:30, ]
target <- datasets::quakes[31:33, ]
train_coords <- as.matrix(train[, c("long", "lat")])
train_x <- as.matrix(train[, c("depth", "stations")])
fit <- gwr_en_fit(
  train_x, train$mag, coords = train_coords, k = 18,
  lambda = 0.05, alpha = 0.7,
  control = gwrs_control(n_threads = 1)
)
prediction <- predict(
  fit,
  newdata = as.matrix(target[, c("depth", "stations")]),
  newcoords = as.matrix(target[, c("long", "lat")]),
  k = 18,
  type = "both"
)
prediction$fit
prediction$coefficients
```
