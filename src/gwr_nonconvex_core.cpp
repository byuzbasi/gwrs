#include "gwrs_nonconvex.h"
#include "gwrs_nonconvex_block.h"
#include "gwrs_nonconvex_active_block.h"

#include <Rcpp.h>

#include <atomic>
#include <memory>
#include <vector>

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

namespace {

bool solve_local_matrix(const arma::mat& system,
                        const arma::mat& rhs,
                        arma::mat& solution,
                        int& status) {
  status = 0;
  bool ok = arma::solve(
    solution, system, rhs,
    arma::solve_opts::likely_sympd + arma::solve_opts::no_approx
  );
  if (ok && solution.is_finite()) return true;

  status = 1;
  ok = arma::solve(solution, system, rhs, arma::solve_opts::no_approx);
  if (ok && solution.is_finite()) return true;

  status = 2;
  const double scale = std::max(1.0, arma::norm(system, "inf"));
  solution = arma::pinv(system, scale * 1e-10) * rhs;
  if (solution.is_finite()) return true;

  status = 3;
  solution.zeros(system.n_cols, rhs.n_cols);
  return false;
}

bool solve_local_vector(const arma::mat& system,
                        const arma::vec& rhs,
                        arma::vec& solution,
                        int& status) {
  arma::mat matrix_solution;
  const bool ok = solve_local_matrix(
    system, arma::mat(rhs), matrix_solution, status
  );
  solution = matrix_solution.col(0);
  return ok;
}

double local_objective(const arma::vec& residual,
                       const arma::vec& normalized_weights,
                       const arma::vec& beta,
                       const double lambda,
                       const double gamma,
                       const int penalty) {
  double value = 0.5 * arma::dot(
    normalized_weights % residual, residual
  );
  for (arma::uword column = 0; column < beta.n_elem; ++column) {
    value += gwrs::nonconvex_penalty_value(
      beta[column], lambda, gamma, penalty
    );
  }
  return value;
}

// Compact per-target/path counters; no trajectory storage or global settings.
struct BlockCounts {
  arma::imat checks, accepted;
  arma::mat gain;
  BlockCounts(arma::uword n, arma::uword l)
    : checks(n, l, arma::fill::zeros), accepted(n, l, arma::fill::zeros),
      gain(n, l, arma::fill::zeros) {}
};

struct NonconvexPathWorker : public RcppParallel::Worker {
  const arma::mat& x_train;
  const arma::vec& y_train;
  const arma::mat& x_target;
  const arma::umat& neighbor_index;
  const arma::mat& neighbor_distance;
  const int kernel;
  const bool adaptive;
  const double fixed_bandwidth;
  const arma::vec& lambda;
  const int penalty;
  const double gamma;
  const double tolerance;
  const int max_iterations;
  const bool screening;
  const bool keep_coefficients;
  const bool return_diagnostics;
  const bool active_block_enabled;

  arma::mat& predictions;
  arma::imat& nonzero;
  arma::imat& iterations;
  arma::imat& converged;
  arma::mat& stationarity;
  arma::mat& coordinate_gap;
  arma::mat& objective;
  arma::imat& screened_count;
  arma::imat& active_count;
  arma::imat& full_scans;
  arma::imat& pair_attempts;
  arma::imat& pair_accepted;
  arma::mat& pair_gain;
  arma::imat& path_state;
  arma::ivec& status;
  arma::vec& weight_sum;
  arma::vec& bandwidth_used;
  arma::cube& coefficients;
  arma::mat& hat_diag;
  arma::mat& s2_diag;
  arma::mat& post_predictions;
  arma::mat& post_hat_diag;
  arma::mat& post_s2_diag;
  arma::imat& local_convex;
  arma::mat& minimum_hessian;
  arma::imat& penalty_knot_count;
  std::atomic<unsigned int>* surface_active_count;
  BlockCounts& blocks;

  NonconvexPathWorker(
      const arma::mat& x_train,
      const arma::vec& y_train,
      const arma::mat& x_target,
      const arma::umat& neighbor_index,
      const arma::mat& neighbor_distance,
      const int kernel,
      const bool adaptive,
      const double fixed_bandwidth,
      const arma::vec& lambda,
      const int penalty,
      const double gamma,
      const double tolerance,
      const int max_iterations,
      const bool screening,
      const bool keep_coefficients,
      const bool return_diagnostics,
      const bool active_block_enabled,
      arma::mat& predictions,
      arma::imat& nonzero,
      arma::imat& iterations,
      arma::imat& converged,
      arma::mat& stationarity,
      arma::mat& coordinate_gap,
      arma::mat& objective,
      arma::imat& screened_count,
      arma::imat& active_count,
      arma::imat& full_scans,
      arma::imat& pair_attempts,
      arma::imat& pair_accepted,
      arma::mat& pair_gain,
      arma::imat& path_state,
      arma::ivec& status,
      arma::vec& weight_sum,
      arma::vec& bandwidth_used,
      arma::cube& coefficients,
      arma::mat& hat_diag,
      arma::mat& s2_diag,
      arma::mat& post_predictions,
      arma::mat& post_hat_diag,
      arma::mat& post_s2_diag,
      arma::imat& local_convex,
      arma::mat& minimum_hessian,
      arma::imat& penalty_knot_count,
      std::atomic<unsigned int>* surface_active_count,
      BlockCounts& blocks)
    : x_train(x_train), y_train(y_train), x_target(x_target),
      neighbor_index(neighbor_index), neighbor_distance(neighbor_distance),
      kernel(kernel), adaptive(adaptive), fixed_bandwidth(fixed_bandwidth),
      lambda(lambda), penalty(penalty), gamma(gamma), tolerance(tolerance),
      max_iterations(max_iterations), screening(screening),
      keep_coefficients(keep_coefficients),
      return_diagnostics(return_diagnostics), active_block_enabled(active_block_enabled), predictions(predictions),
      nonzero(nonzero), iterations(iterations), converged(converged),
      stationarity(stationarity), coordinate_gap(coordinate_gap),
      objective(objective), screened_count(screened_count),
      active_count(active_count), full_scans(full_scans),
      pair_attempts(pair_attempts), pair_accepted(pair_accepted),
      pair_gain(pair_gain), path_state(path_state), status(status),
      weight_sum(weight_sum), bandwidth_used(bandwidth_used),
      coefficients(coefficients), hat_diag(hat_diag), s2_diag(s2_diag),
      post_predictions(post_predictions), post_hat_diag(post_hat_diag),
      post_s2_diag(post_s2_diag), local_convex(local_convex),
      minimum_hessian(minimum_hessian),
      penalty_knot_count(penalty_knot_count),
      surface_active_count(surface_active_count), blocks(blocks) {}

  void operator()(std::size_t begin, std::size_t end) {
    const arma::uword k = neighbor_index.n_rows;
    const arma::uword p = x_train.n_cols;
    const arma::uword n_lambda = lambda.n_elem;

    arma::mat x_local(k, p, arma::fill::zeros);
    arma::mat x_centered(k, p, arma::fill::zeros);
    arma::vec y_local(k, arma::fill::zeros);
    arma::vec y_centered(k, arma::fill::zeros);
    arma::vec weights(k, arma::fill::zeros);
    arma::vec normalized_weights(k, arma::fill::zeros);
    arma::vec sqrt_weights(k, arma::fill::zeros);
    arma::vec curvature(p, arma::fill::zeros);
    arma::vec beta(p, arma::fill::zeros);
    arma::vec residual(k, arma::fill::zeros);
    arma::vec previous_score(p, arma::fill::zeros);
    arma::vec scalar_deltas(active_block_enabled ? p : 0, arma::fill::zeros);
    std::vector<gwrs::NonconvexCoordinateCache> coordinate_cache;
    coordinate_cache.reserve(p);
    arma::rowvec x_mean(p, arma::fill::zeros);
    std::vector<unsigned char> active(p, 0);
    std::vector<unsigned char> ever_active(p, 0);

    for (std::size_t target = begin; target < end; ++target) {
      bool bad_index = false;
      for (arma::uword row = 0; row < k; ++row) {
        const arma::uword encoded = neighbor_index(row, target);
        if (encoded == 0 || encoded > x_train.n_rows) {
          bad_index = true;
          break;
        }
        const arma::uword source = encoded - 1;
        y_local[row] = y_train[source];
        x_local.row(row) = x_train.row(source);
      }
      if (bad_index) {
        status[target] = 4;
        continue;
      }

      double local_bandwidth = 0.0;
      const double local_weight_sum = gwrs::compute_kernel_weights(
        neighbor_distance, target, kernel, adaptive, fixed_bandwidth,
        weights, local_bandwidth
      );
      weight_sum[target] = local_weight_sum;
      bandwidth_used[target] = local_bandwidth;
      if (!std::isfinite(local_weight_sum) ||
          local_weight_sum <= gwrs::kNumericalEpsilon) {
        status[target] = 3;
        continue;
      }

      normalized_weights = weights / local_weight_sum;
      sqrt_weights = arma::sqrt(normalized_weights);
      const double y_mean = arma::dot(normalized_weights, y_local);
      x_mean = normalized_weights.t() * x_local;
      y_centered = y_local - y_mean;
      x_centered = x_local;
      x_centered.each_row() -= x_mean;
      for (arma::uword column = 0; column < p; ++column) {
        curvature[column] = arma::dot(
          normalized_weights, arma::square(x_centered.col(column))
        );
      }

      beta.zeros();
      residual = y_centered;
      previous_score = x_centered.t() * (normalized_weights % residual);
      std::fill(active.begin(), active.end(), static_cast<unsigned char>(0));
      std::fill(ever_active.begin(), ever_active.end(),
                static_cast<unsigned char>(0));
      int target_status = 0;

      // New cache for every target; never shared across locations or workers.
      gwrs::ActiveBlockQRCache qr_cache;

      for (arma::uword path_index = 0; path_index < n_lambda; ++path_index) {
        const double current_lambda = lambda[path_index];
        int total_iterations = 0;
        int scan_count = 0;
        bool coefficient_convergence = false;
        bool needs_refit = false;
        coordinate_cache.clear();
        if (current_lambda > 0.0) {
          for (arma::uword column = 0; column < p; ++column) {
            coordinate_cache.emplace_back(curvature[column], current_lambda, gamma, penalty);
          }
        }

        if (current_lambda <= 0.0) {
          arma::mat weighted_x = x_centered;
          weighted_x.each_col() %= sqrt_weights;
          const arma::mat gram = weighted_x.t() * weighted_x;
          const arma::vec score = x_centered.t() *
            (normalized_weights % y_centered);
          int solve_status = 0;
          gwrs::solve_local_system(gram, score, beta, solve_status);
          target_status = std::max(target_status, solve_status);
          residual = y_centered - x_centered * beta;
          std::fill(active.begin(), active.end(),
                    static_cast<unsigned char>(1));
          std::fill(ever_active.begin(), ever_active.end(),
                    static_cast<unsigned char>(1));
          coefficient_convergence = true;
        } else {
          if (screening) {
            const double threshold = path_index == 0 ? current_lambda :
              std::max(0.0, 2.0 * current_lambda - lambda[path_index - 1]);
            for (arma::uword column = 0; column < p; ++column) {
              active[column] = (ever_active[column] ||
                                std::abs(beta[column]) > 1e-12 ||
                                std::abs(previous_score[column]) >= threshold)
                ? 1 : 0;
            }
          } else {
            std::fill(active.begin(), active.end(),
                      static_cast<unsigned char>(1));
          }

          int initially_active = 0;
          for (arma::uword column = 0; column < p; ++column) {
            initially_active += active[column] ? 1 : 0;
          }
          screened_count(target, path_index) =
            static_cast<int>(p) - initially_active;
          needs_refit = true;

          while (needs_refit && total_iterations < max_iterations) {
            needs_refit = false;
            coefficient_convergence = false;
            while (total_iterations < max_iterations) {
              double maximum_delta = 0.0;
              double maximum_beta = 0.0;
              scalar_deltas.zeros();
              for (arma::uword column = 0; column < p; ++column) {
                if (!active[column]) continue;
                const double local_curvature = curvature[column];
                if (!std::isfinite(local_curvature) ||
                    local_curvature <= gwrs::kNumericalEpsilon) {
                  if (beta[column] != 0.0) {
                    residual += x_centered.col(column) * beta[column];
                    beta[column] = 0.0;
                  }
                  continue;
                }
                const double partial_score = arma::dot(
                  x_centered.col(column), normalized_weights % residual
                ) + local_curvature * beta[column];
                const double updated = gwrs::nonconvex_coordinate_minimum(
                  partial_score, local_curvature, current_lambda, gamma,
                  penalty, &coordinate_cache[column]
                );
                const double delta = updated - beta[column];
                if (active_block_enabled) scalar_deltas[column] = std::abs(delta);
                if (delta != 0.0) {
                  beta[column] = updated;
                  residual -= x_centered.col(column) * delta;
                  maximum_delta = std::max(maximum_delta, std::abs(delta));
                }
                if (std::abs(beta[column]) > 1e-12) {
                  ever_active[column] = 1;
                }
                maximum_beta = std::max(maximum_beta,
                                        std::abs(beta[column]));
              }
              ++total_iterations;
              const bool scalar_stopping =
                maximum_delta <= tolerance * (1.0 + maximum_beta);
              if (active_block_enabled &&
                  (total_iterations % gwrs::kPairEverySweeps == 0 || scalar_stopping)) {
                arma::uword first = 0, second = 0;
                double cross_curvature = 0.0;
                if (gwrs::select_nonconvex_pair(x_centered, normalized_weights,
                      curvature, beta, scalar_deltas, first, second, cross_curvature)) {
                  ++pair_attempts(target, path_index);
                  const double z1 = arma::dot(x_centered.col(first), normalized_weights % residual) +
                    curvature[first] * beta[first] + cross_curvature * beta[second];
                  const double z2 = arma::dot(x_centered.col(second), normalized_weights % residual) +
                    curvature[second] * beta[second] + cross_curvature * beta[first];
                  const auto proposal = gwrs::nonconvex_pair_candidate(
                    curvature[first], cross_curvature, curvature[second], z1, z2,
                    current_lambda, gamma, penalty, beta[first], beta[second]);
                  double gain = 0.0, block_delta = 0.0;
                  if (gwrs::apply_nonconvex_pair(x_centered, normalized_weights, beta,
                        residual, first, second, proposal, current_lambda, gamma,
                        penalty, gain, block_delta)) {
                    ++pair_accepted(target, path_index);
                    pair_gain(target, path_index) += gain;
                    maximum_delta = std::max(maximum_delta, block_delta);
                    maximum_beta = arma::abs(beta).max();
                    if (std::abs(beta[first]) > 1e-12) ever_active[first] = 1;
                    if (std::abs(beta[second]) > 1e-12) ever_active[second] = 1;
                  }
                }
                if (active_block_enabled) {
                  const auto block = gwrs::apply_nonconvex_active_block(
                    x_centered, normalized_weights, beta, residual, current_lambda, gamma, penalty,
                    &qr_cache);
                  ++blocks.checks(target, path_index);
                  if (block.status == 7) {
                    ++blocks.accepted(target, path_index);
                    blocks.gain(target, path_index) += block.gain;
                    maximum_delta = std::max(maximum_delta, block.delta);
                    maximum_beta = arma::abs(beta).max();
                  }
                }
              }
              if (maximum_delta <= tolerance * (1.0 + maximum_beta)) {
                coefficient_convergence = true;
                break;
              }
            }

            if (screening && coefficient_convergence) {
              ++scan_count;
              for (arma::uword column = 0; column < p; ++column) {
                if (active[column]) continue;
                const double local_curvature = curvature[column];
                if (!std::isfinite(local_curvature) ||
                    local_curvature <= gwrs::kNumericalEpsilon) continue;
                const double partial_score = arma::dot(
                  x_centered.col(column), normalized_weights % residual
                );
                const double gap = gwrs::nonconvex_coordinate_improvement(
                  0.0, partial_score, local_curvature, current_lambda, gamma,
                  penalty, &coordinate_cache[column]
                );
                const double scale = 1.0 + std::abs(
                  gwrs::coordinate_objective(
                    0.0, partial_score, local_curvature, current_lambda,
                    gamma, penalty
                  )
                );
                if (gap > tolerance * scale) {
                  active[column] = 1;
                  ever_active[column] = 1;
                  needs_refit = true;
                }
              }
            }
          }
        }

        double maximum_stationarity = 0.0;
        double maximum_coordinate_gap = 0.0;
        int local_nonzero = 0;
        int local_active_count = 0;
        for (arma::uword column = 0; column < p; ++column) {
          const double gradient_score = arma::dot(
            x_centered.col(column), normalized_weights % residual
          );
          previous_score[column] = gradient_score;
          maximum_stationarity = std::max(
            maximum_stationarity,
            gwrs::nonconvex_stationarity_violation(
              beta[column], gradient_score, current_lambda, gamma, penalty
            )
          );
          const double partial_score = gradient_score +
            curvature[column] * beta[column];
          maximum_coordinate_gap = std::max(
            maximum_coordinate_gap,
            gwrs::nonconvex_coordinate_improvement(
              beta[column], partial_score, curvature[column], current_lambda,
              gamma, penalty, current_lambda > 0.0 ? &coordinate_cache[column] : nullptr
            )
          );
          if (std::abs(beta[column]) > 1e-10) ++local_nonzero;
          if (active[column]) ++local_active_count;
        }

        const double intercept = y_mean - arma::dot(x_mean, beta.t());
        predictions(target, path_index) = intercept +
          arma::dot(x_target.row(target), beta.t());
        nonzero(target, path_index) = local_nonzero;
        iterations(target, path_index) = total_iterations;
        const double local_objective_value = local_objective(
          residual, normalized_weights, beta, current_lambda, gamma, penalty
        );
        const double convergence_scale = 1.0 +
          std::abs(local_objective_value);
        converged(target, path_index) =
          (coefficient_convergence && !needs_refit &&
           maximum_coordinate_gap <= 10.0 * tolerance * convergence_scale)
          ? 1 : 0;
        stationarity(target, path_index) = maximum_stationarity;
        coordinate_gap(target, path_index) = maximum_coordinate_gap;
        objective(target, path_index) = local_objective_value;
        active_count(target, path_index) = local_active_count;
        full_scans(target, path_index) = scan_count;

        if (keep_coefficients) {
          coefficients(0, target, path_index) = intercept;
          coefficients.slice(path_index).submat(1, target, p, target) = beta;
        }

        // Explicit target-local failure: preserve this solution and do not
        // warm-start later lambdas from it. Later allocated entries stay NaN/0.
        // 0 = unattempted, 1 = converged under the original rule, 2 = failed.
        path_state(target, path_index) = converged(target, path_index) ? 1 : 2;
        if (active_block_enabled && !converged(target, path_index)) break;

        if (!return_diagnostics) continue;

        std::vector<arma::uword> selected;
        selected.reserve(p);
        int knot_count = 0;
        for (arma::uword column = 0; column < p; ++column) {
          if (current_lambda <= 0.0 || std::abs(beta[column]) > 1e-10) {
            selected.push_back(column);
            surface_active_count[column + p * path_index].fetch_add(
              1u, std::memory_order_relaxed
            );
            if (gwrs::nonconvex_penalty_knot(
                  std::abs(beta[column]), current_lambda, gamma, penalty)) {
              ++knot_count;
            }
          }
        }
        penalty_knot_count(target, path_index) = knot_count;

        arma::uvec selected_index(selected.size());
        for (arma::uword index = 0; index < selected_index.n_elem; ++index) {
          selected_index[index] = selected[index];
        }
        arma::mat derivative(p, k, arma::fill::zeros);
        arma::mat post_derivative(p, k, arma::fill::zeros);
        arma::vec post_beta(p, arma::fill::zeros);
        if (selected_index.n_elem > 0) {
          const arma::mat active_x = x_centered.cols(selected_index);
          arma::mat weighted_active_x = active_x;
          weighted_active_x.each_col() %= sqrt_weights;
          const arma::mat active_gram =
            weighted_active_x.t() * weighted_active_x;
          arma::mat hessian = active_gram;
          for (arma::uword index = 0; index < selected_index.n_elem; ++index) {
            const arma::uword column = selected_index[index];
            hessian(index, index) +=
              gwrs::nonconvex_penalty_second_derivative(
                std::abs(beta[column]), current_lambda, gamma, penalty
              );
          }
          arma::vec eigenvalues;
          if (arma::eig_sym(eigenvalues, hessian) &&
              eigenvalues.is_finite()) {
            const double minimum = eigenvalues.min();
            const double scale = std::max(1.0, eigenvalues.max());
            minimum_hessian(target, path_index) = minimum;
            local_convex(target, path_index) =
              minimum > scale * 1e-10 ? 1 : 0;
          }

          arma::mat score_map = active_x.t();
          score_map.each_row() %= normalized_weights.t();
          arma::mat active_derivative;
          int derivative_status = 0;
          solve_local_matrix(
            hessian, score_map, active_derivative, derivative_status
          );
          target_status = std::max(target_status, derivative_status);
          derivative.rows(selected_index) = active_derivative;

          arma::mat active_post_derivative;
          int post_status = 0;
          solve_local_matrix(
            active_gram, score_map, active_post_derivative, post_status
          );
          post_derivative.rows(selected_index) = active_post_derivative;
          const arma::vec active_score = score_map * y_centered;
          arma::vec active_beta;
          solve_local_vector(
            active_gram, active_score, active_beta, post_status
          );
          post_beta.elem(selected_index) = active_beta;
          target_status = std::max(target_status, post_status);
        } else {
          minimum_hessian(target, path_index) = arma::datum::inf;
          local_convex(target, path_index) = 1;
        }

        const arma::rowvec focus_centered =
          x_target.row(target) - x_mean;
        const arma::vec smoother = normalized_weights +
          derivative.t() * focus_centered.t();
        const arma::vec post_smoother = normalized_weights +
          post_derivative.t() * focus_centered.t();
        double leverage = 0.0;
        double post_leverage = 0.0;
        for (arma::uword row = 0; row < k; ++row) {
          if (neighbor_index(row, target) == target + 1) {
            leverage += smoother[row];
            post_leverage += post_smoother[row];
          }
        }
        hat_diag(target, path_index) = leverage;
        s2_diag(target, path_index) = arma::dot(smoother, smoother);
        post_hat_diag(target, path_index) = post_leverage;
        post_s2_diag(target, path_index) =
          arma::dot(post_smoother, post_smoother);
        const double post_intercept = y_mean -
          arma::dot(x_mean, post_beta.t());
        post_predictions(target, path_index) = post_intercept +
          arma::dot(x_target.row(target), post_beta.t());
      }
      status[target] = target_status;
    }
  }
};

double coordinate_zero_anchor(const double score,
                              const double curvature,
                              const double gamma,
                              const int penalty) {
  const double absolute_score = std::abs(score);
  if (!std::isfinite(absolute_score) ||
      !std::isfinite(curvature) ||
      absolute_score <= gwrs::kNumericalEpsilon ||
      curvature <= gwrs::kNumericalEpsilon) return 0.0;

  auto has_zero_minimum = [&](const double lambda_value) {
    const double solution = gwrs::nonconvex_coordinate_minimum(
      score, curvature, lambda_value, gamma, penalty
    );
    return std::abs(solution) <= 1e-12 *
      (1.0 + absolute_score / curvature);
  };

  double lower = 0.0;
  double upper = absolute_score;
  for (int iteration = 0; iteration < 64 && !has_zero_minimum(upper);
       ++iteration) {
    lower = upper;
    upper *= 2.0;
    if (!std::isfinite(upper)) {
      return arma::datum::inf;
    }
  }
  if (!has_zero_minimum(upper)) {
    return arma::datum::inf;
  }

  for (int iteration = 0; iteration < 80; ++iteration) {
    const double middle = lower + 0.5 * (upper - lower);
    if (has_zero_minimum(middle)) {
      upper = middle;
    } else {
      lower = middle;
    }
  }
  return upper * (1.0 + 1e-10);
}

struct NonconvexLambdaMaxWorker : public RcppParallel::Worker {
  const arma::mat& x;
  const arma::vec& y;
  const arma::umat& neighbor_index;
  const arma::mat& neighbor_distance;
  const int kernel;
  const bool adaptive;
  const double fixed_bandwidth;
  const int penalty;
  const double gamma;
  arma::vec& local_maximum;

  NonconvexLambdaMaxWorker(const arma::mat& x,
                           const arma::vec& y,
                           const arma::umat& neighbor_index,
                           const arma::mat& neighbor_distance,
                           const int kernel,
                           const bool adaptive,
                           const double fixed_bandwidth,
                           const int penalty,
                           const double gamma,
                           arma::vec& local_maximum)
    : x(x), y(y), neighbor_index(neighbor_index),
      neighbor_distance(neighbor_distance), kernel(kernel),
      adaptive(adaptive), fixed_bandwidth(fixed_bandwidth),
      penalty(penalty), gamma(gamma), local_maximum(local_maximum) {}

  void operator()(std::size_t begin, std::size_t end) {
    const arma::uword k = neighbor_index.n_rows;
    const arma::uword p = x.n_cols;
    arma::mat x_local(k, p, arma::fill::zeros);
    arma::mat x_centered(k, p, arma::fill::zeros);
    arma::vec y_local(k, arma::fill::zeros);
    arma::vec weights(k, arma::fill::zeros);
    arma::vec normalized_weights(k, arma::fill::zeros);
    arma::rowvec x_mean(p, arma::fill::zeros);

    for (std::size_t target = begin; target < end; ++target) {
      bool valid = true;
      for (arma::uword row = 0; row < k; ++row) {
        const arma::uword encoded = neighbor_index(row, target);
        if (encoded == 0 || encoded > x.n_rows) {
          valid = false;
          break;
        }
        const arma::uword source = encoded - 1;
        y_local[row] = y[source];
        x_local.row(row) = x.row(source);
      }
      if (!valid) continue;

      double bandwidth_used = 0.0;
      const double weight_sum = gwrs::compute_kernel_weights(
        neighbor_distance, target, kernel, adaptive, fixed_bandwidth,
        weights, bandwidth_used
      );
      if (!std::isfinite(weight_sum) ||
          weight_sum <= gwrs::kNumericalEpsilon) continue;

      normalized_weights = weights / weight_sum;
      const double y_mean = arma::dot(normalized_weights, y_local);
      x_mean = normalized_weights.t() * x_local;
      x_centered = x_local;
      x_centered.each_row() -= x_mean;
      const arma::vec y_centered = y_local - y_mean;

      double maximum = 0.0;
      for (arma::uword column = 0; column < p; ++column) {
        const double curvature = arma::dot(
          normalized_weights, arma::square(x_centered.col(column))
        );
        const double score = arma::dot(
          x_centered.col(column), normalized_weights % y_centered
        );
        maximum = std::max(
          maximum,
          coordinate_zero_anchor(score, curvature, gamma, penalty)
        );
      }
      local_maximum[target] = maximum;
    }
  }
};

void validate_nonconvex_path_inputs(
    const arma::mat& x_train,
    const arma::vec& y_train,
    const arma::mat& x_target,
    const arma::umat& neighbor_index,
    const arma::mat& neighbor_distance,
    const int kernel,
    const bool adaptive,
    const double fixed_bandwidth,
    const arma::vec& lambda,
    const int penalty,
    const double gamma,
    const double tolerance,
    const int max_iterations) {
  if (y_train.n_elem != x_train.n_rows) {
    Rcpp::stop("x_train and y_train have incompatible dimensions");
  }
  if (x_target.n_cols != x_train.n_cols) {
    Rcpp::stop("training and target design matrices must have equal columns");
  }
  if (neighbor_index.n_cols != x_target.n_rows ||
      neighbor_distance.n_cols != x_target.n_rows ||
      neighbor_index.n_rows != neighbor_distance.n_rows) {
    Rcpp::stop("neighbor matrices must have k rows and one column per target");
  }
  if (!x_train.is_finite() || !y_train.is_finite() ||
      !x_target.is_finite() || !neighbor_distance.is_finite()) {
    Rcpp::stop("model inputs and neighbor distances must be finite");
  }
  if (kernel < 0 || kernel > 4) Rcpp::stop("unknown kernel code");
  if (!adaptive && (!std::isfinite(fixed_bandwidth) ||
                    fixed_bandwidth <= 0.0)) {
    Rcpp::stop("fixed bandwidth must be positive");
  }
  if (lambda.n_elem < 1 || !lambda.is_finite() || arma::any(lambda < 0.0)) {
    Rcpp::stop("lambda must be a finite nonnegative vector");
  }
  for (arma::uword index = 1; index < lambda.n_elem; ++index) {
    if (lambda[index] > lambda[index - 1]) {
      Rcpp::stop("lambda must be ordered from largest to smallest");
    }
  }
  if (!gwrs::valid_nonconvex_penalty(penalty)) {
    Rcpp::stop("unknown nonconvex penalty code");
  }
  if (!std::isfinite(gamma) ||
      (penalty == gwrs::kMcpPenalty && gamma <= 1.0) ||
      (penalty == gwrs::kScadPenalty && gamma <= 2.0)) {
    Rcpp::stop("gamma is outside the admissible range for SCAD/MCP");
  }
  if (!std::isfinite(tolerance) || tolerance <= 0.0 ||
      max_iterations < 1) {
    Rcpp::stop("invalid solver tolerance or iteration limit");
  }
}

} // anonymous namespace

// Internal scalar oracle entry point used by package tests. It deliberately
// shares the same exact candidate evaluator as the path solver.
// [[Rcpp::export]]
double cpp_nonconvex_coordinate_minimum(const double score,
                                        const double curvature,
                                        const double lambda,
                                        const int penalty,
                                        const double gamma) {
  if (!gwrs::valid_nonconvex_penalty(penalty)) {
    Rcpp::stop("unknown nonconvex penalty code");
  }
  if (!std::isfinite(gamma) ||
      (penalty == gwrs::kMcpPenalty && gamma <= 1.0) ||
      (penalty == gwrs::kScadPenalty && gamma <= 2.0)) {
    Rcpp::stop("gamma is outside the admissible range for SCAD/MCP");
  }
  return gwrs::nonconvex_coordinate_minimum(
    score, curvature, lambda, gamma, penalty
  );
}

// [[Rcpp::export]]
double cpp_gwr_nonconvex_lambda_max(
    const arma::mat& x,
    const arma::vec& y,
    const arma::umat& neighbor_index,
    const arma::mat& neighbor_distance,
    const int kernel,
    const bool adaptive,
    const double fixed_bandwidth,
    const int penalty,
    const double gamma,
    const int n_threads,
    const int grain_size) {
  const arma::vec validation_lambda(1, arma::fill::zeros);
  validate_nonconvex_path_inputs(
    x, y, x, neighbor_index, neighbor_distance, kernel, adaptive,
    fixed_bandwidth, validation_lambda, penalty, gamma, 1e-7, 1
  );
  arma::vec local_maximum(x.n_rows, arma::fill::zeros);
  NonconvexLambdaMaxWorker worker(
    x, y, neighbor_index, neighbor_distance, kernel, adaptive,
    fixed_bandwidth, penalty, gamma, local_maximum
  );
  if (n_threads == 1 || x.n_rows < 2) {
    worker(0, x.n_rows);
  } else {
    RcppParallel::parallelFor(
      0, x.n_rows, worker,
      static_cast<std::size_t>(std::max(1, grain_size)), n_threads
    );
  }
  const double result = local_maximum.max();
  if (!std::isfinite(result)) {
    Rcpp::stop("failed to compute a finite nonconvex lambda anchor");
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::List cpp_gwr_nonconvex_path_predict(
    const arma::mat& x_train,
    const arma::vec& y_train,
    const arma::mat& x_target,
    const arma::umat& neighbor_index,
    const arma::mat& neighbor_distance,
    const int kernel,
    const bool adaptive,
    const double fixed_bandwidth,
    const arma::vec& lambda,
    const int penalty,
    const double gamma,
    const double tolerance,
    const int max_iterations,
    const bool screening,
    const bool keep_coefficients,
    const bool return_diagnostics,
    const int n_threads,
    const int grain_size,
    const bool guarded_block = false) {
  validate_nonconvex_path_inputs(
    x_train, y_train, x_target, neighbor_index, neighbor_distance, kernel,
    adaptive, fixed_bandwidth, lambda, penalty, gamma, tolerance,
    max_iterations
  );

  const arma::uword n_target = x_target.n_rows;
  const arma::uword n_lambda = lambda.n_elem;
  const arma::uword p = x_train.n_cols;
  arma::mat predictions(n_target, n_lambda,
                        arma::fill::value(arma::datum::nan));
  arma::imat nonzero(n_target, n_lambda, arma::fill::zeros);
  arma::imat iterations(n_target, n_lambda, arma::fill::zeros);
  arma::imat converged(n_target, n_lambda, arma::fill::zeros);
  arma::mat stationarity(n_target, n_lambda,
                         arma::fill::value(arma::datum::nan));
  arma::mat coordinate_gap(n_target, n_lambda,
                           arma::fill::value(arma::datum::nan));
  arma::mat objective(n_target, n_lambda,
                      arma::fill::value(arma::datum::nan));
  arma::imat screened_count(n_target, n_lambda, arma::fill::zeros);
  arma::imat active_count(n_target, n_lambda, arma::fill::zeros);
  arma::imat full_scans(n_target, n_lambda, arma::fill::zeros);
  arma::imat pair_attempts(n_target, n_lambda, arma::fill::zeros);
  arma::imat pair_accepted(n_target, n_lambda, arma::fill::zeros);
  arma::mat pair_gain(n_target, n_lambda, arma::fill::zeros);
  arma::imat path_state(n_target, n_lambda, arma::fill::zeros);
  arma::ivec status(n_target, arma::fill::zeros);
  arma::vec weight_sum(n_target, arma::fill::value(arma::datum::nan));
  arma::vec bandwidth_used(n_target, arma::fill::value(arma::datum::nan));

  arma::cube coefficients;
  if (keep_coefficients) {
    coefficients.set_size(p + 1, n_target, n_lambda);
    coefficients.fill(arma::datum::nan);
  }

  arma::mat hat_diag;
  arma::mat s2_diag;
  arma::mat post_predictions;
  arma::mat post_hat_diag;
  arma::mat post_s2_diag;
  arma::imat local_convex;
  arma::mat minimum_hessian;
  arma::imat penalty_knot_count;
  std::unique_ptr<std::atomic<unsigned int>[]> selected_count;
  if (return_diagnostics) {
    hat_diag.set_size(n_target, n_lambda);
    s2_diag.set_size(n_target, n_lambda);
    post_predictions.set_size(n_target, n_lambda);
    post_hat_diag.set_size(n_target, n_lambda);
    post_s2_diag.set_size(n_target, n_lambda);
    local_convex.set_size(n_target, n_lambda);
    minimum_hessian.set_size(n_target, n_lambda);
    penalty_knot_count.set_size(n_target, n_lambda);
    hat_diag.fill(arma::datum::nan);
    s2_diag.fill(arma::datum::nan);
    post_predictions.fill(arma::datum::nan);
    post_hat_diag.fill(arma::datum::nan);
    post_s2_diag.fill(arma::datum::nan);
    local_convex.zeros();
    minimum_hessian.fill(arma::datum::nan);
    penalty_knot_count.zeros();
    selected_count.reset(new std::atomic<unsigned int>[p * n_lambda]);
    for (arma::uword index = 0; index < p * n_lambda; ++index) {
      selected_count[index].store(0u, std::memory_order_relaxed);
    }
  }

  BlockCounts blocks(n_target, n_lambda);
  NonconvexPathWorker worker(
    x_train, y_train, x_target, neighbor_index, neighbor_distance, kernel,
    adaptive, fixed_bandwidth, lambda, penalty, gamma, tolerance,
    max_iterations, screening, keep_coefficients, return_diagnostics, guarded_block,
    predictions, nonzero, iterations, converged, stationarity,
    coordinate_gap, objective, screened_count, active_count, full_scans,
    pair_attempts, pair_accepted, pair_gain, path_state,
    status, weight_sum, bandwidth_used, coefficients, hat_diag, s2_diag,
    post_predictions, post_hat_diag, post_s2_diag, local_convex,
    minimum_hessian, penalty_knot_count, selected_count.get(), blocks
  );
  if (n_threads == 1 || n_target < 2) {
    worker(0, n_target);
  } else {
    RcppParallel::parallelFor(
      0, n_target, worker,
      static_cast<std::size_t>(std::max(1, grain_size)), n_threads
    );
  }

  arma::umat surface_active_count;
  if (return_diagnostics) {
    surface_active_count.set_size(p, n_lambda);
    for (arma::uword path_index = 0; path_index < n_lambda; ++path_index) {
      for (arma::uword column = 0; column < p; ++column) {
        surface_active_count(column, path_index) =
          selected_count[column + p * path_index].load(
            std::memory_order_relaxed
          );
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::_ ["predictions"] = predictions,
    Rcpp::_ ["nonzero"] = nonzero,
    Rcpp::_ ["iterations"] = iterations,
    Rcpp::_ ["converged"] = converged,
    Rcpp::_ ["stationarity"] = stationarity,
    Rcpp::_ ["coordinate_gap"] = coordinate_gap,
    Rcpp::_ ["objective"] = objective,
    Rcpp::_ ["screened_count"] = screened_count,
    Rcpp::_ ["active_count"] = active_count,
    Rcpp::_ ["full_scans"] = full_scans,
    Rcpp::_ ["pair_attempts"] = pair_attempts,
    Rcpp::_ ["pair_accepted"] = pair_accepted,
    Rcpp::_ ["pair_gain"] = pair_gain,
    Rcpp::_ ["path_state"] = path_state,
    Rcpp::_ ["status"] = status,
    Rcpp::_ ["weight_sum"] = weight_sum,
    Rcpp::_ ["bandwidth"] = bandwidth_used,
    Rcpp::_ ["coefficients"] = coefficients,
    Rcpp::_ ["hat_diag_path"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(hat_diag)) : R_NilValue,
    Rcpp::_ ["s2_diag_path"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(s2_diag)) : R_NilValue,
    Rcpp::_ ["post_fitted_path"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(post_predictions)) : R_NilValue,
    Rcpp::_ ["post_hat_diag_path"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(post_hat_diag)) : R_NilValue,
    Rcpp::_ ["post_s2_diag_path"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(post_s2_diag)) : R_NilValue,
    Rcpp::_ ["local_convex"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(local_convex)) : R_NilValue,
    Rcpp::_ ["minimum_hessian"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(minimum_hessian)) : R_NilValue,
    Rcpp::_ ["penalty_knot_count"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(penalty_knot_count)) : R_NilValue,
    Rcpp::_ ["surface_active_count"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(surface_active_count)) : R_NilValue,
    Rcpp::_ ["block_checks"] = blocks.checks,
    Rcpp::_ ["block_accepted"] = blocks.accepted,
    Rcpp::_ ["block_gain"] = blocks.gain,
    Rcpp::_ ["active_block_policy"] = Rcpp::List::create(
      Rcpp::_ ["enabled"] = guarded_block, Rcpp::_ ["min_size"] = gwrs::kActiveBlockMin,
      Rcpp::_ ["max_size"] = gwrs::kActiveBlockMax, Rcpp::_ ["every_sweeps"] = gwrs::kPairEverySweeps,
      Rcpp::_ ["schema"] = "gwrs-flat-tail-qr-v14",
      Rcpp::_ ["candidate_region_guard"] = false)
  );
}
