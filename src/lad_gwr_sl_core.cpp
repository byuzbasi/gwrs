#include "gwrs_common.h"

#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

namespace {

struct LadGwrSlWorker : public RcppParallel::Worker {
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
  const double rho_residual;
  const double rho_penalty;
  const double tolerance;
  const int max_iterations;

  arma::mat& coefficients;
  arma::mat& gwr_center;
  arma::vec& fitted;
  arma::vec& residuals;
  arma::vec& objective;
  arma::vec& primal_residual;
  arma::vec& dual_residual;
  arma::vec& weight_sum;
  arma::vec& bandwidth_used;
  arma::ivec& iterations;
  arma::ivec& converged;
  arma::ivec& status;

  LadGwrSlWorker(const arma::mat& x,
                 const arma::vec& y,
                 const arma::umat& neighbor_index,
                 const arma::mat& neighbor_distance,
                 const int kernel,
                 const bool adaptive,
                 const double fixed_bandwidth,
                 const double lambda,
                 const double alpha,
                 const double d,
                 const double rho_residual,
                 const double rho_penalty,
                 const double tolerance,
                 const int max_iterations,
                 arma::mat& coefficients,
                 arma::mat& gwr_center,
                 arma::vec& fitted,
                 arma::vec& residuals,
                 arma::vec& objective,
                 arma::vec& primal_residual,
                 arma::vec& dual_residual,
                 arma::vec& weight_sum,
                 arma::vec& bandwidth_used,
                 arma::ivec& iterations,
                 arma::ivec& converged,
                 arma::ivec& status)
    : x(x), y(y), neighbor_index(neighbor_index),
      neighbor_distance(neighbor_distance), kernel(kernel),
      adaptive(adaptive), fixed_bandwidth(fixed_bandwidth), lambda(lambda),
      alpha(alpha), d(d), rho_residual(rho_residual),
      rho_penalty(rho_penalty), tolerance(tolerance),
      max_iterations(max_iterations), coefficients(coefficients),
      gwr_center(gwr_center), fitted(fitted), residuals(residuals),
      objective(objective), primal_residual(primal_residual),
      dual_residual(dual_residual), weight_sum(weight_sum),
      bandwidth_used(bandwidth_used), iterations(iterations),
      converged(converged), status(status) {}

  void operator()(std::size_t begin, std::size_t end) {
    const arma::uword k = neighbor_index.n_rows;
    const arma::uword p = x.n_cols;
    const arma::uword p1 = p + 1;
    const double lambda1 = lambda * alpha;
    const double lambda2 = lambda * (1.0 - alpha);

    arma::mat x_local(k, p, arma::fill::zeros);
    arma::mat x_centered(k, p, arma::fill::zeros);
    arma::mat weighted_x(k, p, arma::fill::zeros);
    arma::mat gram(p, p, arma::fill::zeros);
    arma::mat design(k, p1, arma::fill::ones);
    arma::mat system(p1, p1, arma::fill::zeros);
    arma::mat cholesky;
    arma::mat inverse_system;
    arma::vec y_local(k, arma::fill::zeros);
    arma::vec weights(k, arma::fill::zeros);
    arma::vec normalized_weights(k, arma::fill::zeros);
    arma::vec sqrt_weights(k, arma::fill::zeros);
    arma::vec score(p, arma::fill::zeros);
    arma::vec beta_gwr(p, arma::fill::zeros);
    arma::vec theta(p1, arma::fill::zeros);
    arma::vec center(p1, arma::fill::zeros);
    arma::vec residual_split(k, arma::fill::zeros);
    arma::vec residual_previous(k, arma::fill::zeros);
    arma::vec coefficient_split(p1, arma::fill::zeros);
    arma::vec coefficient_previous(p1, arma::fill::zeros);
    arma::vec dual_residual_split(k, arma::fill::zeros);
    arma::vec dual_coefficient_split(p1, arma::fill::zeros);
    arma::vec rhs(p1, arma::fill::zeros);
    arma::vec temporary(p1, arma::fill::zeros);
    arma::vec constraint_residual(k, arma::fill::zeros);
    arma::vec constraint_coefficient(p1, arma::fill::zeros);
    arma::rowvec x_mean(p, arma::fill::zeros);

    for (std::size_t target = begin; target < end; ++target) {
      bool bad_index = false;
      for (arma::uword row = 0; row < k; ++row) {
        const arma::uword encoded = neighbor_index(row, target);
        if (encoded == 0 || encoded > x.n_rows) {
          bad_index = true;
          break;
        }
        const arma::uword source = encoded - 1;
        y_local[row] = y[source];
        x_local.row(row) = x.row(source);
      }
      if (bad_index) {
        status[target] = 4;
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
        continue;
      }

      normalized_weights = weights / local_weight_sum;
      const double y_mean = arma::dot(normalized_weights, y_local);
      x_mean = normalized_weights.t() * x_local;
      x_centered = x_local;
      x_centered.each_row() -= x_mean;
      sqrt_weights = arma::sqrt(normalized_weights);
      weighted_x = x_centered;
      weighted_x.each_col() %= sqrt_weights;
      gram = weighted_x.t() * weighted_x;
      score = x_centered.t() *
        (normalized_weights % (y_local - y_mean));

      int local_status = 0;
      gwrs::solve_local_system(gram, score, beta_gwr, local_status);
      status[target] = local_status;
      gwr_center.col(target) = beta_gwr;

      design.col(0).ones();
      design.cols(1, p) = x_local;
      center.zeros();
      center.subvec(1, p) = d * beta_gwr;
      theta.zeros();
      theta[0] = y_mean - arma::dot(x_mean, beta_gwr.t());
      theta.subvec(1, p) = beta_gwr;
      residual_split = y_local - design * theta;
      coefficient_split = theta;
      dual_residual_split.zeros();
      dual_coefficient_split.zeros();

      system = rho_residual * (design.t() * design);
      system.diag() += rho_penalty;
      for (arma::uword column = 1; column < p1; ++column) {
        system(column, column) += lambda2;
      }
      const bool chol_ok = arma::chol(cholesky, system, "lower");
      if (!chol_ok) inverse_system = arma::pinv(system);

      bool did_converge = false;
      int iteration = 0;
      double local_primal = arma::datum::inf;
      double local_dual = arma::datum::inf;
      for (; iteration < max_iterations; ++iteration) {
        residual_previous = residual_split;
        coefficient_previous = coefficient_split;

        rhs = rho_residual * design.t() *
          (y_local - residual_split - dual_residual_split) +
          rho_penalty *
          (coefficient_split - dual_coefficient_split) +
          lambda2 * center;
        if (chol_ok) {
          temporary = arma::solve(arma::trimatl(cholesky), rhs,
                                  arma::solve_opts::fast);
          theta = arma::solve(arma::trimatu(cholesky.t()), temporary,
                              arma::solve_opts::fast);
        } else {
          theta = inverse_system * rhs;
        }

        const arma::vec residual_argument =
          y_local - design * theta - dual_residual_split;
        for (arma::uword row = 0; row < k; ++row) {
          residual_split[row] = gwrs::soft_threshold(
            residual_argument[row], normalized_weights[row] / rho_residual);
        }

        coefficient_split[0] = theta[0] + dual_coefficient_split[0];
        for (arma::uword column = 1; column < p1; ++column) {
          coefficient_split[column] = gwrs::soft_threshold(
            theta[column] + dual_coefficient_split[column],
            lambda1 / rho_penalty);
        }

        constraint_residual = residual_split - y_local + design * theta;
        constraint_coefficient = theta - coefficient_split;
        dual_residual_split += constraint_residual;
        dual_coefficient_split += constraint_coefficient;

        local_primal = std::sqrt(
          arma::dot(constraint_residual, constraint_residual) +
          arma::dot(constraint_coefficient, constraint_coefficient));
        const arma::vec residual_dual = rho_residual * design.t() *
          (residual_split - residual_previous);
        const arma::vec coefficient_dual = rho_penalty *
          (coefficient_split - coefficient_previous);
        local_dual = std::sqrt(
          arma::dot(residual_dual, residual_dual) +
          arma::dot(coefficient_dual, coefficient_dual));

        const double primal_scale = 1.0 +
          arma::norm(residual_split, 2) + arma::norm(y_local, 2) +
          arma::norm(design * theta, 2) + arma::norm(theta, 2) +
          arma::norm(coefficient_split, 2);
        const double dual_scale = 1.0 +
          arma::norm(rho_residual * design.t() * dual_residual_split, 2) +
          arma::norm(rho_penalty * dual_coefficient_split, 2);
        if (local_primal <= tolerance * primal_scale &&
            local_dual <= tolerance * dual_scale) {
          did_converge = true;
          ++iteration;
          break;
        }
      }

      iterations[target] = iteration;
      converged[target] = did_converge ? 1 : 0;
      primal_residual[target] = local_primal;
      dual_residual[target] = local_dual;
      coefficients.col(target) = theta;
      const double target_fitted = theta[0] +
        arma::dot(x.row(target), theta.subvec(1, p).t());
      fitted[target] = target_fitted;
      residuals[target] = y[target] - target_fitted;

      const arma::vec local_fit_residual = y_local - design * theta;
      const arma::vec center_difference = theta.subvec(1, p) -
        d * beta_gwr;
      objective[target] =
        arma::dot(normalized_weights, arma::abs(local_fit_residual)) +
        lambda1 * arma::accu(arma::abs(theta.subvec(1, p))) +
        0.5 * lambda2 * arma::dot(center_difference, center_difference);
    }
  }
};

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List cpp_lad_gwr_sl_fit(const arma::mat& x,
                              const arma::vec& y,
                              const arma::umat& neighbor_index,
                              const arma::mat& neighbor_distance,
                              const int kernel,
                              const bool adaptive,
                              const double fixed_bandwidth,
                              const double lambda,
                              const double alpha,
                              const double d,
                              const double rho_residual,
                              const double rho_penalty,
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
  if (!adaptive && (!std::isfinite(fixed_bandwidth) ||
                    fixed_bandwidth <= 0.0)) {
    Rcpp::stop("fixed bandwidth must be positive");
  }
  if (!std::isfinite(lambda) || lambda < 0.0 ||
      !std::isfinite(alpha) || alpha < 0.0 || alpha > 1.0 ||
      !std::isfinite(d) || d < 0.0 || d > 1.0 ||
      !std::isfinite(rho_residual) || rho_residual <= 0.0 ||
      !std::isfinite(rho_penalty) || rho_penalty <= 0.0) {
    Rcpp::stop("penalty or ADMM parameters are outside their valid ranges");
  }

  arma::mat coefficients(p + 1, n, arma::fill::value(arma::datum::nan));
  arma::mat gwr_center(p, n, arma::fill::value(arma::datum::nan));
  arma::vec fitted(n, arma::fill::value(arma::datum::nan));
  arma::vec residuals(n, arma::fill::value(arma::datum::nan));
  arma::vec objective(n, arma::fill::value(arma::datum::nan));
  arma::vec primal_residual(n, arma::fill::value(arma::datum::nan));
  arma::vec dual_residual(n, arma::fill::value(arma::datum::nan));
  arma::vec weight_sum(n, arma::fill::value(arma::datum::nan));
  arma::vec bandwidth_used(n, arma::fill::value(arma::datum::nan));
  arma::ivec iterations(n, arma::fill::zeros);
  arma::ivec converged(n, arma::fill::zeros);
  arma::ivec status(n, arma::fill::zeros);

  LadGwrSlWorker worker(
    x, y, neighbor_index, neighbor_distance, kernel, adaptive,
    fixed_bandwidth, lambda, alpha, d, rho_residual, rho_penalty, tolerance,
    max_iterations, coefficients, gwr_center, fitted, residuals, objective,
    primal_residual, dual_residual, weight_sum, bandwidth_used, iterations,
    converged, status);
  if (n_threads == 1 || n < 2) {
    worker(0, n);
  } else {
    RcppParallel::parallelFor(
      0, n, worker,
      static_cast<std::size_t>(std::max(1, grain_size)), n_threads);
  }

  return Rcpp::List::create(
    Rcpp::_["coefficients"] = coefficients,
    Rcpp::_["gwr_center"] = gwr_center,
    Rcpp::_["fitted"] = fitted,
    Rcpp::_["residuals"] = residuals,
    Rcpp::_["objective"] = objective,
    Rcpp::_["primal_residual"] = primal_residual,
    Rcpp::_["dual_residual"] = dual_residual,
    Rcpp::_["weight_sum"] = weight_sum,
    Rcpp::_["bandwidth"] = bandwidth_used,
    Rcpp::_["iterations"] = iterations,
    Rcpp::_["converged"] = converged,
    Rcpp::_["status"] = status
  );
}
