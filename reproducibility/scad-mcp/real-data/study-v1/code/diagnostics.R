# Independent convex diagnostics. These functions inspect a completed native
# iterate and never call a fitting routine or alter the optimization state.

usa_study_v1_kernel_weights <- function(distance, kernel, adaptive,
                                    fixed_bandwidth) {
  bandwidth <- if (adaptive) max(distance) * (1 + 1e-10) else fixed_bandwidth
  if (!is.finite(bandwidth) || bandwidth <= 1e-12) bandwidth <- 1
  ratio <- distance / bandwidth
  weight <- switch(
    as.character(kernel),
    `0` = exp(-0.5 * ratio^2),
    `1` = ifelse(ratio < 1, (1 - ratio^2)^2, 0),
    `2` = exp(-ratio),
    `3` = ifelse(ratio < 1, (1 - ratio^3)^3, 0),
    `4` = ifelse(ratio <= 1, 1, 0),
    stop("Unknown kernel code in independent diagnostics.")
  )
  stopifnot(all(is.finite(weight)), sum(weight) > 0)
  weight / sum(weight)
}

usa_study_v1_duality <- function(x, y, weight, beta, lambda1, lambda2) {
  x_centered <- sweep(x, 2, colSums(x * weight))
  y_centered <- y - sum(y * weight)
  residual <- y_centered - drop(x_centered %*% beta)
  theta <- sqrt(weight) * residual
  weighted_y <- sqrt(weight) * y_centered
  score <- drop(crossprod(x_centered, weight * residual))
  primal <- sum(weight * residual^2) / 2 +
    lambda1 * sum(abs(beta)) + lambda2 * sum(beta^2) / 2
  if (lambda2 == 0) {
    multiplier <- if (max(abs(score)) == 0) 1 else
      min(1, lambda1 / max(abs(score)))
    theta <- multiplier * theta
    dual <- sum(weighted_y * theta) - sum(theta^2) / 2
  } else {
    conjugate <- sum(pmax(abs(score) - lambda1, 0)^2) / (2 * lambda2)
    dual <- sum(weighted_y * theta) - sum(theta^2) / 2 - conjugate
  }
  gap <- primal - dual
  stopifnot(is.finite(gap), gap >= -1e-10)
  c(duality_gap = gap,
    relative_duality_gap = gap / (1 + abs(primal)))
}

usa_study_v1_working_metrics <- function(record) {
  args <- record$args
  raw <- record$raw
  diagnostics <- raw$working_set_diagnostics
  rows <- vector("list", nrow(raw$converged) * length(args$lambda))
  position <- 0L
  for (target in seq_len(nrow(raw$converged))) {
    index <- args$neighbor_index[, target]
    distance <- args$neighbor_distance[, target]
    weight <- usa_study_v1_kernel_weights(
      distance, args$kernel, args$adaptive, args$fixed_bandwidth
    )
    x <- args$x_train[index, , drop = FALSE]
    y <- args$y_train[index]
    for (lambda_index in seq_along(args$lambda)) {
      position <- position + 1L
      beta <- raw$coefficients[-1, target, lambda_index]
      intercept <- raw$coefficients[1, target, lambda_index]
      residual <- y - intercept - drop(x %*% beta)
      lambda1 <- args$lambda[lambda_index] * args$alpha
      lambda2 <- args$lambda[lambda_index] * (1 - args$alpha)
      score <- drop(crossprod(x, weight * residual)) - lambda2 * beta
      kkt <- max(ifelse(
        abs(beta) <= 1e-12,
        pmax(abs(score) - lambda1, 0),
        abs(score - lambda1 * sign(beta))
      ))
      objective <- sum(weight * residual^2) / 2 +
        lambda1 * sum(abs(beta)) + lambda2 * sum(beta^2) / 2
      predicted <- intercept + sum(args$x_target[target, ] * beta)
      dual <- usa_study_v1_duality(x, y, weight, beta, lambda1, lambda2)
      stopifnot(
        abs(kkt - raw$kkt[target, lambda_index]) <=
          1e-9 * (1 + abs(kkt)),
        abs(predicted - raw$predictions[target, lambda_index]) <=
          1e-10 * (1 + abs(predicted)),
        abs(sum(weight * residual)) <= 1e-10 * (1 + max(abs(y))),
        raw$nonzero[target, lambda_index] == sum(abs(beta) > 1e-12),
        raw$iterations[target, lambda_index] <= args$max_iterations
      )
      at <- c(target, lambda_index)
      value <- data.frame(
        cap = args$max_iterations,
        target = target,
        index = lambda_index,
        lambda = args$lambda[lambda_index],
        converged = raw$converged[target, lambda_index],
        iterations = raw$iterations[target, lambda_index],
        independent_kkt = kkt,
        native_kkt = raw$kkt[target, lambda_index],
        objective = objective,
        nonzero = sum(abs(beta) > 1e-12),
        working_attempts = diagnostics$attempts[target, lambda_index],
        working_accepted = diagnostics$accepted[target, lambda_index],
        working_gain = diagnostics$gain[target, lambda_index],
        rank_rejections = diagnostics$rank_rejected[target, lambda_index],
        selected_size = diagnostics$max_selected_size[target, lambda_index],
        observed_support = diagnostics$max_observed_support[target, lambda_index],
        qr_factorizations = diagnostics$qr_factorizations[target, lambda_index],
        inner_updates = diagnostics$inner_updates[target, lambda_index],
        boundary_hits = diagnostics$boundary_hits[target, lambda_index],
        zero_entries = diagnostics$zero_entries[target, lambda_index],
        inner_limit_hits = diagnostics$inner_limit_hits[target, lambda_index],
        certificate_checks = diagnostics$certificate_checks[target, lambda_index],
        certificate_failures = diagnostics$certificate_failures[target, lambda_index],
        final_polish_attempts = diagnostics$final_polish_attempts[target, lambda_index],
        final_polish_accepted = diagnostics$final_polish_accepted[target, lambda_index],
        max_certificate_kkt = diagnostics$max_certificate_kkt[target, lambda_index],
        max_certificate_ratio = diagnostics$max_certificate_ratio[target, lambda_index],
        kkt_bound = args$tolerance * (1 + lambda1),
        duality_gap = unname(dual[["duality_gap"]]),
        relative_duality_gap = unname(dual[["relative_duality_gap"]]),
        screen_expansions = diagnostics$screen_expansions[target, lambda_index],
        screen_coordinates_added = diagnostics$screen_coordinates_added[target, lambda_index],
        duality_certificate_checks = diagnostics$duality_certificate_checks[target, lambda_index],
        duality_certificate_failures = diagnostics$duality_certificate_failures[target, lambda_index],
        max_relative_duality_gap = diagnostics$max_relative_duality_gap[target, lambda_index],
        max_duality_ratio = diagnostics$max_duality_ratio[target, lambda_index],
        screen_fit_events = diagnostics$screen_fit_events[target, lambda_index],
        screen_fit_batches = diagnostics$screen_fit_batches[target, lambda_index],
        screen_fit_budget_hits = diagnostics$screen_fit_budget_hits[target, lambda_index],
        max_screen_fit_batches = diagnostics$max_screen_fit_batches[target, lambda_index]
      )
      value$kkt_certificate_pass <-
        value$independent_kkt <= value$kkt_bound * (1 + 1e-10)
      value$native_duality_bound <- max(
        10 * args$tolerance, 256 * .Machine$double.eps
      )
      value$native_duality_certificate_pass <-
        value$relative_duality_gap <=
          value$native_duality_bound * (1 + 1e-10)
      value$native_joint_certificate_pass <-
        value$kkt_certificate_pass & value$native_duality_certificate_pass
      rows[[position]] <- value
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  stopifnot(
    all(result$final_polish_attempts == 0),
    all(result$final_polish_accepted == 0),
    all(result$screen_expansions >= 0),
    all(result$screen_coordinates_added >= result$screen_expansions),
    all(result$screen_coordinates_added <=
          ncol(args$x_train) * result$screen_expansions),
    all(result$screen_fit_batches == result$working_attempts),
    all(result$screen_fit_events <= result$screen_fit_batches),
    all(result$screen_fit_batches <= 16 * result$screen_fit_events),
    all(result$screen_fit_budget_hits <= result$screen_fit_events),
    all(result$max_screen_fit_batches <= 16),
    all(result$max_screen_fit_batches <= result$screen_fit_batches),
    all(result$duality_certificate_failures <=
          result$duality_certificate_checks),
    all(result$duality_certificate_checks >= result$converged),
    all(result$native_joint_certificate_pass[result$converged == 1])
  )
  result
}
