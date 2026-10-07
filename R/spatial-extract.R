gwrs_spatial_coordinates <- function(object) {
  if (!inherits(object, "gwrs_fit")) {
    stop("`object` must inherit from 'gwrs_fit'.", call. = FALSE)
  }
  if (is.null(object$neighbors)) {
    stop("The fitted object does not retain its spatial neighbor structure.",
         call. = FALSE)
  }
  coordinates <- object$neighbors$query_coords %||%
    object$neighbors$coords
  coordinates <- as_coordinate_matrix(coordinates)
  if (ncol(coordinates) < 2L) {
    stop("Heat-map coordinates must contain at least two columns.",
         call. = FALSE)
  }
  if (nrow(coordinates) != nrow(object$coefficients)) {
    stop("Stored coordinates and local coefficients have different row counts.",
         call. = FALSE)
  }
  coordinates[, 1:2, drop = FALSE]
}

gwrs_coefficient_names <- function(object) {
  names <- colnames(object$coefficients)
  if (!is.null(names) && length(names) == ncol(object$coefficients)) {
    return(names)
  }
  predictors <- object$predictor_names
  if (is.null(predictors) ||
      length(predictors) != ncol(object$coefficients) - 1L) {
    predictors <- paste0("x", seq_len(ncol(object$coefficients) - 1L))
  }
  c("(Intercept)", predictors)
}

gwrs_spatial_variable_options <- function(object) {
  coefficient_names <- gwrs_coefficient_names(object)
  unique(c(
    coefficient_names,
    if (length(coefficient_names) > 1L) {
      paste0("beta_", seq_len(length(coefficient_names) - 1L))
    },
    "beta_0", "local_r2", "residuals", "fitted", "response",
    "abs_residuals",
    if (!is.null(object$process)) "process"
  ))
}

resolve_gwrs_spatial_variable <- function(object, variable) {
  if (!is.character(variable) || length(variable) != 1L ||
      is.na(variable) || !nzchar(trimws(variable))) {
    stop("`variable` must be one non-empty character value.", call. = FALSE)
  }
  variable <- trimws(variable)
  coefficient_names <- gwrs_coefficient_names(object)

  exact <- match(variable, coefficient_names)
  if (!is.na(exact)) {
    return(list(
      type = "coefficient", column = exact,
      name = coefficient_names[exact], label = coefficient_names[exact]
    ))
  }
  folded <- which(tolower(coefficient_names) == tolower(variable))
  if (length(folded) == 1L) {
    return(list(
      type = "coefficient", column = folded,
      name = coefficient_names[folded], label = coefficient_names[folded]
    ))
  }
  if (length(folded) > 1L) {
    stop(
      sprintf(
        "`variable = %s` is ambiguous. Matching coefficients: %s.",
        sQuote(variable), paste(coefficient_names[folded], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  normalized <- tolower(gsub("[ .-]+", "_", variable))
  if (normalized %in% c("intercept", "(intercept)", "beta_0", "beta0")) {
    return(list(
      type = "coefficient", column = 1L,
      name = coefficient_names[1L], label = coefficient_names[1L]
    ))
  }
  beta_match <- regexec("^beta_?([0-9]+)$", normalized)
  beta_parts <- regmatches(normalized, beta_match)[[1L]]
  if (length(beta_parts) == 2L) {
    slope <- as.integer(beta_parts[2L])
    column <- slope + 1L
    if (!is.na(column) && column >= 1L && column <= length(coefficient_names)) {
      return(list(
        type = "coefficient", column = column,
        name = coefficient_names[column], label = coefficient_names[column]
      ))
    }
  }

  special <- switch(
    normalized,
    local_r2 = list(type = "local_r2", name = "local_r2",
                    label = "Local R-squared"),
    local_r_squared = list(type = "local_r2", name = "local_r2",
                           label = "Local R-squared"),
    residual = list(type = "residuals", name = "residuals",
                    label = "Residuals"),
    residuals = list(type = "residuals", name = "residuals",
                     label = "Residuals"),
    fitted = list(type = "fitted", name = "fitted", label = "Fitted"),
    fitted_values = list(type = "fitted", name = "fitted",
                         label = "Fitted"),
    response = list(type = "response", name = "response", label = "Response"),
    y = list(type = "response", name = "response", label = "Response"),
    abs_residual = list(type = "abs_residuals", name = "abs_residuals",
                        label = "Absolute residuals"),
    abs_residuals = list(type = "abs_residuals", name = "abs_residuals",
                         label = "Absolute residuals"),
    process = list(type = "process", name = "process", label = "Process"),
    NULL
  )
  if (!is.null(special)) return(special)

  choices <- gwrs_spatial_variable_options(object)
  stop(
    sprintf(
      "Unknown spatial variable %s. Valid choices include: %s.",
      sQuote(variable), paste(choices, collapse = ", ")
    ),
    call. = FALSE
  )
}

gwrs_local_r2_values <- function(object) {
  value <- object$model_diagnostics$local$local_r_squared %||%
    object$local_r2
  if (is.null(value)) {
    stop(
      "Local R-squared is unavailable; refit with statistical diagnostics enabled.",
      call. = FALSE
    )
  }
  as.double(value)
}

gwrs_response_values <- function(object) {
  if (!is.null(object$y)) return(as.double(object$y))
  if (!is.null(object$fitted.values) && !is.null(object$residuals)) {
    return(as.double(object$fitted.values + object$residuals))
  }
  stop("The response cannot be reconstructed from this fitted object.",
       call. = FALSE)
}

gwrs_is_sparse_selection_fit <- function(object) {
  regularized <- !is.null(object$lambda) && length(object$lambda) == 1L &&
    is.finite(object$lambda) && object$lambda > 0
  if (!regularized) return(FALSE)
  if ((object$penalty %||% "") %in% c("scad", "mcp")) return(TRUE)
  !is.null(object$alpha) && length(object$alpha) == 1L &&
    is.finite(object$alpha) && object$alpha > 0
}

extract_gwrs_spatial_values <- function(object,
                                        variable,
                                        selection_tolerance = 1e-10,
                                        na.rm = TRUE) {
  coordinates <- gwrs_spatial_coordinates(object)
  resolved <- resolve_gwrs_spatial_variable(object, variable)
  values <- switch(
    resolved$type,
    coefficient = as.double(object$coefficients[, resolved$column]),
    local_r2 = gwrs_local_r2_values(object),
    residuals = as.double(object$residuals),
    abs_residuals = abs(as.double(object$residuals)),
    fitted = as.double(object$fitted.values),
    response = gwrs_response_values(object),
    process = {
      if (is.null(object$process)) {
        stop("A process surface is available only for GWR-KTDD fits.",
             call. = FALSE)
      }
      as.double(object$process)
    }
  )
  if (length(values) != nrow(coordinates)) {
    stop("The selected spatial value has an incompatible length.",
         call. = FALSE)
  }
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm)) {
    stop("`na.rm` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(selection_tolerance) ||
      length(selection_tolerance) != 1L ||
      !is.finite(selection_tolerance) || selection_tolerance < 0) {
    stop("`selection_tolerance` must be one nonnegative finite number.",
         call. = FALSE)
  }

  selected <- rep(TRUE, length(values))
  if (identical(resolved$type, "coefficient") &&
      resolved$column > 1L && gwrs_is_sparse_selection_fit(object)) {
    selected <- abs(values) > selection_tolerance
  }
  points <- data.frame(
    x = coordinates[, 1L],
    y = coordinates[, 2L],
    value = values,
    observation_id = seq_along(values),
    selected = selected,
    stringsAsFactors = FALSE
  )
  invalid <- !is.finite(points$value)
  if (any(invalid) && !na.rm) {
    stop("The requested spatial variable contains non-finite values.",
         call. = FALSE)
  }
  if (any(invalid)) points <- points[!invalid, , drop = FALSE]
  if (nrow(points) < 3L) {
    stop("At least three finite spatial values are required.", call. = FALSE)
  }
  unique_location <- !duplicated(points[c("x", "y")])
  if (sum(unique_location) < 3L) {
    stop("At least three unique spatial locations are required.",
         call. = FALSE)
  }
  list(
    points = points,
    variable = resolved$name,
    label = resolved$label,
    value_type = resolved$type,
    coefficient = identical(resolved$type, "coefficient"),
    coefficient_column = resolved$column %||% NA_integer_,
    original_geometry = NULL
  )
}

add_gwrs_spatial_significance <- function(extracted,
                                          object,
                                          alpha,
                                          p_adjust,
                                          estimator) {
  if (!isTRUE(extracted$coefficient)) {
    stop("Significance masking is available only for local coefficient surfaces.",
         call. = FALSE)
  }
  inference <- gwr_local_inference(
    object,
    alpha = alpha,
    adjust = p_adjust,
    estimator = estimator
  )
  local <- inference$local[
    inference$local$coefficient == extracted$variable,
    , drop = FALSE
  ]
  if (nrow(local) != nrow(object$coefficients)) {
    stop("Local inference could not be matched to the requested coefficient.",
         call. = FALSE)
  }
  match_index <- match(extracted$points$observation_id, local$location)
  extracted$points$p_value <- local$p_value[match_index]
  extracted$points$p_adjusted <- local$p_adjusted[match_index]
  extracted$points$significant <- local$significant[match_index]
  extracted$inference <- inference$settings
  extracted
}
