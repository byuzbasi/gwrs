#include "gwrs_nonconvex.h"

#include <Rcpp.h>

#include <vector>

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

namespace {

bool solve_symmetric_matrix(const arma::mat& system,
                            const arma::mat& rhs,
                            arma::mat& solution,
                            int& status) {
  status = 0;
  bool ok = arma::solve(solution, system, rhs,
                        arma::solve_opts::likely_sympd +
                          arma::solve_opts::no_approx);
  if (ok && solution.is_finite()) return true;
  status = 1;
  const double scale = std::max(1.0, arma::norm(system, "inf"));
  solution = arma::pinv(system, scale * 1e-10) * rhs;
  if (solution.is_finite()) return true;
  status = 2;
  solution.zeros(system.n_cols, rhs.n_cols);
  return false;
}

bool solve_symmetric_vector(const arma::mat& system,
                            const arma::vec& rhs,
                            arma::vec& solution,
                            int& status) {
  arma::mat matrix_solution;
  const bool ok = solve_symmetric_matrix(
    system, arma::mat(rhs), matrix_solution, status
  );
  solution = matrix_solution.col(0);
  return ok;
}

arma::uvec conditional_active_set(const arma::mat& coefficients,
                                  const arma::uword target,
                                  const arma::uword p,
                                  const double lambda,
                                  const double alpha,
                                  const int penalty_mode) {
  if (lambda <= 0.0 || (penalty_mode == 0 && alpha <= 0.0)) {
    return arma::regspace<arma::uvec>(0, p - 1);
  }
  std::vector<arma::uword> active;
  active.reserve(p);
  for (arma::uword column = 0; column < p; ++column) {
    if (std::abs(coefficients(column + 1, target)) > 1e-10) {
      active.push_back(column);
    }
  }
  arma::uvec result(active.size());
  for (arma::uword index = 0; index < result.n_elem; ++index) {
    result[index] = active[index];
  }
  return result;
}

arma::mat original_scale_map(const arma::mat& standard_map,
                             const arma::vec& x_center,
                             const arma::vec& x_scale) {
  const arma::uword p = x_center.n_elem;
  arma::mat result = standard_map;
  arma::rowvec intercept = standard_map.row(0);
  for (arma::uword column = 0; column < p; ++column) {
    const arma::rowvec slope = standard_map.row(column + 1) / x_scale[column];
    result.row(column + 1) = slope;
    intercept -= x_center[column] * slope;
  }
  result.row(0) = intercept;
  return result;
}

struct GwrDiagnosticWorker : public RcppParallel::Worker {
  const arma::mat& x;
  const arma::vec& y;
  const arma::vec& fitted;
  const arma::mat& coefficients;
  const arma::umat& neighbor_index;
  const arma::mat& neighbor_distance;
  const arma::vec& x_center;
  const arma::vec& x_scale;
  const int kernel;
  const bool adaptive;
  const double fixed_bandwidth;
  const double lambda;
  const double alpha;
  const double d;
  const int penalty_mode;
  const double gamma;
  const bool return_inference;

  arma::vec& hat_diag;
  arma::vec& s2_diag;
  arma::vec& local_r2;
  arma::vec& local_neff;
  arma::ivec& local_n_positive;
  arma::vec& local_kappa;
  arma::ivec& local_rank;
  arma::ivec& operator_status;
  arma::ivec& local_convex;
  arma::vec& local_min_hessian;
  arma::ivec& local_penalty_knot_count;
  arma::mat& coefficient_map_ss;
  arma::mat& post_coefficient_map_ss;
  arma::mat& post_coefficients;
  arma::vec& post_fitted;
  arma::vec& post_hat_diag;
  arma::vec& post_s2_diag;

  GwrDiagnosticWorker(
      const arma::mat& x,
      const arma::vec& y,
      const arma::vec& fitted,
      const arma::mat& coefficients,
      const arma::umat& neighbor_index,
      const arma::mat& neighbor_distance,
      const arma::vec& x_center,
      const arma::vec& x_scale,
      const int kernel,
      const bool adaptive,
      const double fixed_bandwidth,
      const double lambda,
      const double alpha,
      const double d,
      const int penalty_mode,
      const double gamma,
      const bool return_inference,
      arma::vec& hat_diag,
      arma::vec& s2_diag,
      arma::vec& local_r2,
      arma::vec& local_neff,
      arma::ivec& local_n_positive,
      arma::vec& local_kappa,
      arma::ivec& local_rank,
      arma::ivec& operator_status,
      arma::ivec& local_convex,
      arma::vec& local_min_hessian,
      arma::ivec& local_penalty_knot_count,
      arma::mat& coefficient_map_ss,
      arma::mat& post_coefficient_map_ss,
      arma::mat& post_coefficients,
      arma::vec& post_fitted,
      arma::vec& post_hat_diag,
      arma::vec& post_s2_diag)
    : x(x), y(y), fitted(fitted), coefficients(coefficients),
      neighbor_index(neighbor_index), neighbor_distance(neighbor_distance),
      x_center(x_center), x_scale(x_scale), kernel(kernel),
      adaptive(adaptive), fixed_bandwidth(fixed_bandwidth), lambda(lambda),
      alpha(alpha), d(d), penalty_mode(penalty_mode), gamma(gamma),
      return_inference(return_inference),
      hat_diag(hat_diag), s2_diag(s2_diag), local_r2(local_r2),
      local_neff(local_neff), local_n_positive(local_n_positive),
      local_kappa(local_kappa), local_rank(local_rank),
      operator_status(operator_status), local_convex(local_convex),
      local_min_hessian(local_min_hessian),
      local_penalty_knot_count(local_penalty_knot_count),
      coefficient_map_ss(coefficient_map_ss),
      post_coefficient_map_ss(post_coefficient_map_ss),
      post_coefficients(post_coefficients), post_fitted(post_fitted),
      post_hat_diag(post_hat_diag), post_s2_diag(post_s2_diag) {}

  void operator()(std::size_t begin, std::size_t end) {
    const arma::uword k = neighbor_index.n_rows;
    const arma::uword p = x.n_cols;
    const double lambda2 = lambda * (1.0 - alpha);

    arma::mat x_local(k, p, arma::fill::zeros);
    arma::mat x_centered(k, p, arma::fill::zeros);
    arma::mat weighted_x(k, p, arma::fill::zeros);
    arma::mat gram(p, p, arma::fill::zeros);
    arma::mat score_map(p, k, arma::fill::zeros);
    arma::mat gwr_map(p, k, arma::fill::zeros);
    arma::vec y_local(k, arma::fill::zeros);
    arma::vec fitted_local(k, arma::fill::zeros);
    arma::vec weights(k, arma::fill::zeros);
    arma::vec normalized_weights(k, arma::fill::zeros);
    arma::vec sqrt_weights(k, arma::fill::zeros);
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
        x_local.row(row) = x.row(source);
        y_local[row] = y[source];
        fitted_local[row] = fitted[source];
      }
      if (!valid) {
        operator_status[target] = 4;
        continue;
      }

      double bandwidth_used = 0.0;
      const double weight_sum = gwrs::compute_kernel_weights(
        neighbor_distance, target, kernel, adaptive, fixed_bandwidth,
        weights, bandwidth_used
      );
      if (!std::isfinite(weight_sum) ||
          weight_sum <= gwrs::kNumericalEpsilon) {
        operator_status[target] = 3;
        continue;
      }
      normalized_weights = weights / weight_sum;
      local_neff[target] = 1.0 /
        std::max(arma::dot(normalized_weights, normalized_weights),
                 gwrs::kNumericalEpsilon);
      local_n_positive[target] = static_cast<int>(
        arma::accu(weights > gwrs::kNumericalEpsilon)
      );

      const double y_mean = arma::dot(normalized_weights, y_local);
      x_mean = normalized_weights.t() * x_local;
      x_centered = x_local;
      x_centered.each_row() -= x_mean;
      sqrt_weights = arma::sqrt(normalized_weights);
      weighted_x = x_centered;
      weighted_x.each_col() %= sqrt_weights;
      gram = weighted_x.t() * weighted_x;
      score_map = x_centered.t();
      score_map.each_row() %= normalized_weights.t();

      arma::vec eigenvalues;
      if (arma::eig_sym(eigenvalues, gram) && eigenvalues.is_finite()) {
        const double maximum = eigenvalues.max();
        const double tolerance = std::max(1.0, maximum) * 1e-10;
        const arma::uword rank = arma::accu(eigenvalues > tolerance);
        local_rank[target] = static_cast<int>(rank);
        if (rank < p || maximum <= tolerance) {
          local_kappa[target] = arma::datum::inf;
        } else {
          local_kappa[target] = std::sqrt(maximum / eigenvalues.min());
        }
      }

      int local_status = 0;
      if (penalty_mode == 0) {
        solve_symmetric_matrix(gram, score_map, gwr_map, local_status);
      }

      const arma::uvec active = conditional_active_set(
        coefficients, target, p, lambda, alpha, penalty_mode
      );
      arma::mat derivative(p, k, arma::fill::zeros);
      if (active.n_elem > 0) {
        arma::mat active_system = gram.submat(active, active);
        arma::mat active_rhs = score_map.rows(active);
        if (penalty_mode == 0) {
          active_system.diag() += lambda2;
          active_rhs += lambda2 * d * gwr_map.rows(active);
        } else {
          const int nonconvex_penalty = penalty_mode - 1;
          for (arma::uword index = 0; index < active.n_elem; ++index) {
            const double coefficient =
              coefficients(active[index] + 1, target);
            active_system(index, index) +=
              gwrs::nonconvex_penalty_second_derivative(
                std::abs(coefficient), lambda, gamma, nonconvex_penalty
              );
            if (gwrs::nonconvex_penalty_knot(
                  std::abs(coefficient), lambda, gamma,
                  nonconvex_penalty)) {
              ++local_penalty_knot_count[target];
            }
          }
        }
        arma::vec active_eigenvalues;
        if (arma::eig_sym(active_eigenvalues, active_system) &&
            active_eigenvalues.is_finite()) {
          const double minimum = active_eigenvalues.min();
          const double scale = std::max(1.0, active_eigenvalues.max());
          local_min_hessian[target] = minimum;
          local_convex[target] = minimum > scale * 1e-10 ? 1 : 0;
        }
        arma::mat active_derivative;
        int active_status = 0;
        solve_symmetric_matrix(active_system, active_rhs,
                               active_derivative, active_status);
        derivative.rows(active) = active_derivative;
        local_status = std::max(local_status, active_status);
      } else {
        local_convex[target] = 1;
        local_min_hessian[target] = arma::datum::inf;
      }

      const arma::rowvec focus_centered = x.row(target) - x_mean;
      const arma::vec smoother_row = normalized_weights +
        derivative.t() * focus_centered.t();
      s2_diag[target] = arma::dot(smoother_row, smoother_row);
      double leverage = 0.0;
      for (arma::uword row = 0; row < k; ++row) {
        if (neighbor_index(row, target) == target + 1) {
          leverage += smoother_row[row];
        }
      }
      hat_diag[target] = leverage;
      operator_status[target] = local_status;

      const double local_tss = arma::dot(
        normalized_weights, arma::square(y_local - y_mean)
      );
      const double local_rss = arma::dot(
        normalized_weights, arma::square(y_local - fitted_local)
      );
      if (local_tss > gwrs::kNumericalEpsilon) {
        local_r2[target] = 1.0 - local_rss / local_tss;
      }

      if (!return_inference) continue;

      arma::mat coefficient_map(p + 1, k, arma::fill::zeros);
      coefficient_map.rows(1, p) = derivative;
      coefficient_map.row(0) = normalized_weights.t() -
        x_mean * derivative;
      const arma::mat original_map = original_scale_map(
        coefficient_map, x_center, x_scale
      );
      coefficient_map_ss.row(target) =
        arma::sum(arma::square(original_map), 1).t();

      arma::mat post_derivative(p, k, arma::fill::zeros);
      arma::vec post_beta(p, arma::fill::zeros);
      if (active.n_elem > 0) {
        const arma::mat active_gram = gram.submat(active, active);
        arma::mat active_post_derivative;
        int post_status = 0;
        solve_symmetric_matrix(active_gram, score_map.rows(active),
                               active_post_derivative, post_status);
        post_derivative.rows(active) = active_post_derivative;
        const arma::vec active_score = score_map.rows(active) * y_local;
        arma::vec active_beta;
        solve_symmetric_vector(active_gram, active_score, active_beta,
                               post_status);
        post_beta.elem(active) = active_beta;
        local_status = std::max(local_status, post_status);
      }
      operator_status[target] = local_status;
      const double post_intercept = y_mean -
        arma::dot(x_mean, post_beta.t());
      post_coefficients(0, target) = post_intercept;
      post_coefficients.submat(1, target, p, target) = post_beta;
      post_fitted[target] = post_intercept +
        arma::dot(x.row(target), post_beta.t());

      const arma::vec post_smoother = normalized_weights +
        post_derivative.t() * focus_centered.t();
      post_s2_diag[target] = arma::dot(post_smoother, post_smoother);
      double post_leverage = 0.0;
      for (arma::uword row = 0; row < k; ++row) {
        if (neighbor_index(row, target) == target + 1) {
          post_leverage += post_smoother[row];
        }
      }
      post_hat_diag[target] = post_leverage;

      arma::mat post_map(p + 1, k, arma::fill::zeros);
      post_map.rows(1, p) = post_derivative;
      post_map.row(0) = normalized_weights.t() -
        x_mean * post_derivative;
      const arma::mat original_post_map = original_scale_map(
        post_map, x_center, x_scale
      );
      post_coefficient_map_ss.row(target) =
        arma::sum(arma::square(original_post_map), 1).t();
    }
  }
};

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List cpp_gwr_model_diagnostics(
    const arma::mat& x,
    const arma::vec& y,
    const arma::vec& fitted,
    const arma::mat& coefficients,
    const arma::umat& neighbor_index,
    const arma::mat& neighbor_distance,
    const arma::vec& x_center,
    const arma::vec& x_scale,
    const int kernel,
    const bool adaptive,
    const double fixed_bandwidth,
    const double lambda,
    const double alpha,
    const double d,
    const int penalty_mode,
    const double gamma,
    const bool return_inference,
    const int n_threads,
    const int grain_size) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (y.n_elem != n || fitted.n_elem != n ||
      coefficients.n_rows != p + 1 || coefficients.n_cols != n ||
      neighbor_index.n_cols != n || neighbor_distance.n_cols != n ||
      neighbor_index.n_rows != neighbor_distance.n_rows ||
      x_center.n_elem != p || x_scale.n_elem != p) {
    Rcpp::stop("incompatible dimensions in diagnostic inputs");
  }
  if (!x.is_finite() || !y.is_finite() || !fitted.is_finite() ||
      !coefficients.is_finite() || !neighbor_distance.is_finite() ||
      !x_center.is_finite() || !x_scale.is_finite() ||
      arma::any(x_scale <= 0.0)) {
    Rcpp::stop("diagnostic inputs must be finite with positive scales");
  }
  if (penalty_mode < 0 || penalty_mode > 2 ||
      (penalty_mode > 0 &&
       (!std::isfinite(gamma) ||
        (penalty_mode == 1 && gamma <= 1.0) ||
        (penalty_mode == 2 && gamma <= 2.0)))) {
    Rcpp::stop("invalid diagnostic penalty mode or gamma");
  }

  arma::vec hat_diag(n, arma::fill::value(arma::datum::nan));
  arma::vec s2_diag(n, arma::fill::value(arma::datum::nan));
  arma::vec local_r2(n, arma::fill::value(arma::datum::nan));
  arma::vec local_neff(n, arma::fill::value(arma::datum::nan));
  arma::ivec local_n_positive(n, arma::fill::zeros);
  arma::vec local_kappa(n, arma::fill::value(arma::datum::nan));
  arma::ivec local_rank(n, arma::fill::zeros);
  arma::ivec operator_status(n, arma::fill::zeros);
  arma::ivec local_convex(n, arma::fill::zeros);
  arma::vec local_min_hessian(n, arma::fill::value(arma::datum::nan));
  arma::ivec local_penalty_knot_count(n, arma::fill::zeros);
  arma::mat coefficient_map_ss;
  arma::mat post_coefficient_map_ss;
  arma::mat post_coefficients;
  arma::vec post_fitted;
  arma::vec post_hat_diag;
  arma::vec post_s2_diag;
  if (return_inference) {
    coefficient_map_ss.set_size(n, p + 1);
    coefficient_map_ss.fill(arma::datum::nan);
    post_coefficient_map_ss.set_size(n, p + 1);
    post_coefficient_map_ss.fill(arma::datum::nan);
    post_coefficients.set_size(p + 1, n);
    post_coefficients.fill(arma::datum::nan);
    post_fitted.set_size(n);
    post_fitted.fill(arma::datum::nan);
    post_hat_diag.set_size(n);
    post_hat_diag.fill(arma::datum::nan);
    post_s2_diag.set_size(n);
    post_s2_diag.fill(arma::datum::nan);
  }

  GwrDiagnosticWorker worker(
    x, y, fitted, coefficients, neighbor_index, neighbor_distance,
    x_center, x_scale, kernel, adaptive, fixed_bandwidth, lambda, alpha, d,
    penalty_mode, gamma, return_inference, hat_diag, s2_diag, local_r2, local_neff,
    local_n_positive, local_kappa, local_rank, operator_status,
    local_convex, local_min_hessian, local_penalty_knot_count,
    coefficient_map_ss, post_coefficient_map_ss, post_coefficients,
    post_fitted, post_hat_diag, post_s2_diag
  );
  if (n_threads == 1 || n < 2) {
    worker(0, n);
  } else {
    RcppParallel::parallelFor(
      0, n, worker,
      static_cast<std::size_t>(std::max(1, grain_size)), n_threads
    );
  }

  return Rcpp::List::create(
    Rcpp::_ ["hat_diag"] = hat_diag,
    Rcpp::_ ["s2_diag"] = s2_diag,
    Rcpp::_ ["local_r2"] = local_r2,
    Rcpp::_ ["local_neff"] = local_neff,
    Rcpp::_ ["local_n_positive"] = local_n_positive,
    Rcpp::_ ["local_kappa"] = local_kappa,
    Rcpp::_ ["local_rank"] = local_rank,
    Rcpp::_ ["operator_status"] = operator_status,
    Rcpp::_ ["local_convex"] = local_convex,
    Rcpp::_ ["local_min_hessian"] = local_min_hessian,
    Rcpp::_ ["local_penalty_knot_count"] = local_penalty_knot_count,
    Rcpp::_ ["coefficient_map_ss"] = return_inference ?
      static_cast<SEXP>(Rcpp::wrap(coefficient_map_ss)) : R_NilValue,
    Rcpp::_ ["post_coefficient_map_ss"] = return_inference ?
      static_cast<SEXP>(Rcpp::wrap(post_coefficient_map_ss)) : R_NilValue,
    Rcpp::_ ["post_coefficients"] = return_inference ?
      static_cast<SEXP>(Rcpp::wrap(post_coefficients)) : R_NilValue,
    Rcpp::_ ["post_fitted"] = return_inference ?
      static_cast<SEXP>(Rcpp::wrap(post_fitted)) : R_NilValue,
    Rcpp::_ ["post_hat_diag"] = return_inference ?
      static_cast<SEXP>(Rcpp::wrap(post_hat_diag)) : R_NilValue,
    Rcpp::_ ["post_s2_diag"] = return_inference ?
      static_cast<SEXP>(Rcpp::wrap(post_s2_diag)) : R_NilValue
  );
}
