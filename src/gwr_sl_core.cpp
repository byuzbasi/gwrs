#include "gwrs_common.h"
#include "gwrs_ridge.h"

#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

namespace {

struct GwrSlWorker : public RcppParallel::Worker {
  const arma::mat& x;
  const arma::vec& y;
  const arma::umat& neighbor_index;
  const arma::mat& neighbor_distance;
  const int kernel;
  const bool adaptive;
  const double fixed_bandwidth;
  const double lambda;
  const double alpha;
  const double d;
  const double tolerance;
  const int max_iterations;

  arma::mat& coefficients;
  arma::mat& gwr_center;
  arma::vec& fitted;
  arma::vec& residuals;
  arma::vec& objective;
  arma::vec& weight_sum;
  arma::vec& bandwidth_used;
  arma::vec& kkt;
  arma::ivec& iterations;
  arma::ivec& converged;
  arma::ivec& status;

  GwrSlWorker(const arma::mat& x,
              const arma::vec& y,
              const arma::umat& neighbor_index,
              const arma::mat& neighbor_distance,
              const int kernel,
              const bool adaptive,
              const double fixed_bandwidth,
              const double lambda,
              const double alpha,
              const double d,
              const double tolerance,
              const int max_iterations,
              arma::mat& coefficients,
              arma::mat& gwr_center,
              arma::vec& fitted,
              arma::vec& residuals,
              arma::vec& objective,
              arma::vec& weight_sum,
              arma::vec& bandwidth_used,
              arma::vec& kkt,
              arma::ivec& iterations,
              arma::ivec& converged,
              arma::ivec& status)
    : x(x), y(y), neighbor_index(neighbor_index),
      neighbor_distance(neighbor_distance), kernel(kernel),
      adaptive(adaptive), fixed_bandwidth(fixed_bandwidth), lambda(lambda),
      alpha(alpha), d(d), tolerance(tolerance),
      max_iterations(max_iterations), coefficients(coefficients),
      gwr_center(gwr_center), fitted(fitted), residuals(residuals),
      objective(objective), weight_sum(weight_sum),
      bandwidth_used(bandwidth_used), kkt(kkt), iterations(iterations),
      converged(converged), status(status) {}

  void operator()(std::size_t begin, std::size_t end) {
    const arma::uword k = neighbor_index.n_rows;
    const arma::uword p = x.n_cols;
    const double lambda1 = lambda * alpha;
    const double lambda2 = lambda * (1.0 - alpha);
    const bool pure_positive_ridge = lambda > 0.0 && alpha == 0.0 && d == 0.0;
    gwrs::RidgeSystem ridge_system;

    arma::mat x_local(k, p, arma::fill::zeros);
    arma::mat x_centered(k, p, arma::fill::zeros);
    arma::mat weighted_x(k, p, arma::fill::zeros);
    arma::mat gram(p, p, arma::fill::zeros);
    arma::vec y_local(k, arma::fill::zeros);
    arma::vec y_centered(k, arma::fill::zeros);
    arma::vec weights(k, arma::fill::zeros);
    arma::vec normalized_weights(k, arma::fill::zeros);
    arma::vec sqrt_weights(k, arma::fill::zeros);
    arma::vec score(p, arma::fill::zeros);
    arma::vec beta_gwr(p, arma::fill::zeros);
    arma::vec beta(p, arma::fill::zeros);
    arma::vec local_residual(k, arma::fill::zeros);
    arma::rowvec x_mean(p, arma::fill::zeros);

    for (std::size_t target = begin; target < end; ++target) {
      int local_status = 0;
      bool bad_index = false;
      for (arma::uword row = 0; row < k; ++row) {
        const arma::uword encoded = neighbor_index(row, target);
        if (encoded == 0 || encoded > x.n_rows) {
          bad_index = true;
          break;
        }
        const arma::uword source = encoded - 1;
        y_local[row] = y[source];
        for (arma::uword column = 0; column < p; ++column) {
          x_local(row, column) = x(source, column);
        }
      }

      if (bad_index) {
        status[target] = 4;
        converged[target] = 0;
        continue;
      }

      double local_bandwidth = 0.0;
      const double local_weight_sum = gwrs::compute_kernel_weights(
        neighbor_distance, target, kernel, adaptive, fixed_bandwidth,
        weights, local_bandwidth);
      weight_sum[target] = local_weight_sum;
      bandwidth_used[target] = local_bandwidth;

      if (!std::isfinite(local_weight_sum) ||
          local_weight_sum <= gwrs::kNumericalEpsilon) {
        status[target] = 3;
        converged[target] = 0;
        continue;
      }

      normalized_weights = weights / local_weight_sum;
      const double y_mean = arma::dot(normalized_weights, y_local);
      x_mean = normalized_weights.t() * x_local;
      y_centered = y_local - y_mean;
      x_centered = x_local;
      x_centered.each_row() -= x_mean;

      sqrt_weights = arma::sqrt(normalized_weights);
      weighted_x = x_centered;
      weighted_x.each_col() %= sqrt_weights;
      gram = weighted_x.t() * weighted_x;
      score = x_centered.t() * (normalized_weights % y_centered);

      gwrs::solve_local_system(gram, score, beta_gwr, local_status);
      gwr_center.col(target) = beta_gwr;

      if (lambda <= 0.0) {
        beta = beta_gwr;
        iterations[target] = 0;
        converged[target] = 1;
      } else if (pure_positive_ridge) {
        // The legacy GWR-center output above is retained for API compatibility;
        // it does not enter the pure Ridge solve, which has no fallback.
        ridge_system.prepare(weighted_x, sqrt_weights % y_centered);
        if (!ridge_system.solve(weighted_x, lambda, beta)) {
          status[target] = 2;
          converged[target] = 0;
          continue;
        }
        local_status = 0;
        iterations[target] = 1;
        converged[target] = 1;
      } else {
        beta = beta_gwr;
        local_residual = y_centered - x_centered * beta;
        int iteration = 0;
        bool did_converge = false;

        for (; iteration < max_iterations; ++iteration) {
          double maximum_delta = 0.0;
          double maximum_beta = 0.0;
          for (arma::uword column = 0; column < p; ++column) {
            const double curvature = gram(column, column);
            if (!std::isfinite(curvature) ||
                curvature <= gwrs::kNumericalEpsilon) {
              if (beta[column] != 0.0) {
                local_residual += x_centered.col(column) * beta[column];
                beta[column] = 0.0;
              }
              continue;
            }

            const double partial_score =
              arma::dot(x_centered.col(column),
                        normalized_weights % local_residual) +
              curvature * beta[column];
            const double shifted_score = partial_score +
              lambda2 * d * beta_gwr[column];
            const double updated = gwrs::soft_threshold(shifted_score,
                                                        lambda1) /
              (curvature + lambda2);
            const double delta = updated - beta[column];
            if (delta != 0.0) {
              beta[column] = updated;
              local_residual -= x_centered.col(column) * delta;
              maximum_delta = std::max(maximum_delta, std::abs(delta));
            }
            maximum_beta = std::max(maximum_beta, std::abs(beta[column]));
          }

          if (maximum_delta <= tolerance * (1.0 + maximum_beta)) {
            did_converge = true;
            ++iteration;
            break;
          }
        }
        iterations[target] = iteration;
        converged[target] = did_converge ? 1 : 0;
      }

      local_residual = y_centered - x_centered * beta;
      const double intercept = y_mean - arma::dot(x_mean, beta.t());
      coefficients(0, target) = intercept;
      coefficients.submat(1, target, p, target) = beta;

      const double target_fitted = intercept +
        arma::dot(x.row(target), beta.t());
      fitted[target] = target_fitted;
      residuals[target] = y[target] - target_fitted;

      const arma::vec difference = beta - d * beta_gwr;
      objective[target] = 0.5 * arma::dot(normalized_weights % local_residual,
                                          local_residual) +
        lambda1 * arma::accu(arma::abs(beta)) +
        0.5 * lambda2 * arma::dot(difference, difference);
      kkt[target] = gwrs::kkt_residual(x_centered, local_residual,
                                       normalized_weights, beta,
                                       d * beta_gwr, lambda1, lambda2);
      status[target] = local_status;
    }
  }
};

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List cpp_gwr_sl_fit(const arma::mat& x,
                          const arma::vec& y,
                          const arma::umat& neighbor_index,
                          const arma::mat& neighbor_distance,
                          const int kernel,
                          const bool adaptive,
                          const double fixed_bandwidth,
                          const double lambda,
                          const double alpha,
                          const double d,
                          const double tolerance,
                          const int max_iterations,
                          const int n_threads,
                          const int grain_size) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (y.n_elem != n) Rcpp::stop("x and y have incompatible dimensions");
  if (neighbor_index.n_cols != n || neighbor_distance.n_cols != n ||
      neighbor_index.n_rows != neighbor_distance.n_rows) {
    Rcpp::stop("neighbor matrices must have k rows and n columns");
  }
  if (!x.is_finite() || !y.is_finite() || !neighbor_distance.is_finite()) {
    Rcpp::stop("x, y, and neighbor distances must be finite");
  }
  if (kernel < 0 || kernel > 4) Rcpp::stop("unknown kernel code");
  if (!adaptive && (!std::isfinite(fixed_bandwidth) || fixed_bandwidth <= 0.0)) {
    Rcpp::stop("fixed bandwidth must be positive");
  }
  if (lambda < 0.0 || alpha < 0.0 || alpha > 1.0 || d < 0.0 || d > 1.0) {
    Rcpp::stop("lambda, alpha, or d is outside its valid range");
  }

  arma::mat coefficients(p + 1, n, arma::fill::value(arma::datum::nan));
  arma::mat gwr_center(p, n, arma::fill::value(arma::datum::nan));
  arma::vec fitted(n, arma::fill::value(arma::datum::nan));
  arma::vec residuals(n, arma::fill::value(arma::datum::nan));
  arma::vec objective(n, arma::fill::value(arma::datum::nan));
  arma::vec weight_sum(n, arma::fill::value(arma::datum::nan));
  arma::vec bandwidth_used(n, arma::fill::value(arma::datum::nan));
  arma::vec kkt(n, arma::fill::value(arma::datum::nan));
  arma::ivec iterations(n, arma::fill::zeros);
  arma::ivec converged(n, arma::fill::zeros);
  arma::ivec status(n, arma::fill::zeros);

  GwrSlWorker worker(x, y, neighbor_index, neighbor_distance, kernel,
                     adaptive, fixed_bandwidth, lambda, alpha, d, tolerance,
                     max_iterations, coefficients, gwr_center, fitted,
                     residuals, objective, weight_sum, bandwidth_used, kkt,
                     iterations, converged, status);

  if (n_threads == 1 || n < 2) {
    worker(0, n);
  } else {
    RcppParallel::parallelFor(0, n, worker,
                              static_cast<std::size_t>(std::max(1, grain_size)),
                              n_threads);
  }

  return Rcpp::List::create(
    Rcpp::_["coefficients"] = coefficients,
    Rcpp::_["gwr_center"] = gwr_center,
    Rcpp::_["fitted"] = fitted,
    Rcpp::_["residuals"] = residuals,
    Rcpp::_["objective"] = objective,
    Rcpp::_["weight_sum"] = weight_sum,
    Rcpp::_["bandwidth"] = bandwidth_used,
    Rcpp::_["kkt"] = kkt,
    Rcpp::_["iterations"] = iterations,
    Rcpp::_["converged"] = converged,
    Rcpp::_["status"] = status
  );
}
