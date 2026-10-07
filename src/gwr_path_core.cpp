#include "gwrs_common.h"
#include "gwrs_ridge.h"

#include <Rcpp.h>

#include <atomic>
#include <memory>
#include <vector>
#include <chrono>

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

namespace {

arma::mat solve_path_matrix(const arma::mat& system,
                            const arma::mat& rhs) {
  arma::mat solution;
  const bool ok = arma::solve(
    solution, system, rhs,
    arma::solve_opts::likely_sympd + arma::solve_opts::no_approx
  );
  if (ok && solution.is_finite()) return solution;
  const double scale = std::max(1.0, arma::norm(system, "inf"));
  return arma::pinv(system, scale * 1e-10) * rhs;
}

struct LambdaMaxWorker : public RcppParallel::Worker {
  const arma::mat& x;
  const arma::vec& y;
  const arma::umat& neighbor_index;
  const arma::mat& neighbor_distance;
  const int kernel;
  const bool adaptive;
  const double fixed_bandwidth;
  const double alpha;
  arma::vec& local_maximum;

  LambdaMaxWorker(const arma::mat& x,
                  const arma::vec& y,
                  const arma::umat& neighbor_index,
                  const arma::mat& neighbor_distance,
                  const int kernel,
                  const bool adaptive,
                  const double fixed_bandwidth,
                  const double alpha,
                  arma::vec& local_maximum)
    : x(x), y(y), neighbor_index(neighbor_index),
      neighbor_distance(neighbor_distance), kernel(kernel),
      adaptive(adaptive), fixed_bandwidth(fixed_bandwidth), alpha(alpha),
      local_maximum(local_maximum) {}

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
        weights, bandwidth_used);
      if (!std::isfinite(weight_sum) ||
          weight_sum <= gwrs::kNumericalEpsilon) continue;

      normalized_weights = weights / weight_sum;
      const double y_mean = arma::dot(normalized_weights, y_local);
      x_mean = normalized_weights.t() * x_local;
      x_centered = x_local;
      x_centered.each_row() -= x_mean;
      const arma::vec score = x_centered.t() *
        (normalized_weights % (y_local - y_mean));
      local_maximum[target] = arma::max(arma::abs(score)) / alpha;
    }
  }
};

struct GwrSlPathWorker : public RcppParallel::Worker {
  const arma::mat& x_train;
  const arma::vec& y_train;
  const arma::mat& x_target;
  const arma::umat& neighbor_index;
  const arma::mat& neighbor_distance;
  const int kernel;
  const bool adaptive;
  const double fixed_bandwidth;
  const arma::vec& lambda;
  const double alpha;
  const double d;
  const double tolerance;
  const int max_iterations;
  const bool screening;
  const bool keep_coefficients;
  const bool return_diagnostics;

  arma::mat& predictions;
  arma::imat& nonzero;
  arma::imat& iterations;
  arma::imat& converged;
  arma::mat& kkt;
  arma::ivec& status;
  arma::vec& weight_sum;
  arma::vec& bandwidth_used;
  arma::cube& coefficients;
  arma::mat& hat_diag;
  arma::mat& s2_diag;
  arma::mat& post_predictions;
  arma::mat& post_hat_diag;
  arma::mat& post_s2_diag;
  std::atomic<unsigned int>* surface_active_count;

  GwrSlPathWorker(const arma::mat& x_train,
                  const arma::vec& y_train,
                  const arma::mat& x_target,
                  const arma::umat& neighbor_index,
                  const arma::mat& neighbor_distance,
                  const int kernel,
                  const bool adaptive,
                  const double fixed_bandwidth,
                  const arma::vec& lambda,
                  const double alpha,
                  const double d,
                  const double tolerance,
                  const int max_iterations,
                  const bool screening,
                  const bool keep_coefficients,
                  const bool return_diagnostics,
                  arma::mat& predictions,
                  arma::imat& nonzero,
                  arma::imat& iterations,
                  arma::imat& converged,
                  arma::mat& kkt,
                  arma::ivec& status,
                  arma::vec& weight_sum,
                  arma::vec& bandwidth_used,
                  arma::cube& coefficients,
                  arma::mat& hat_diag,
                  arma::mat& s2_diag,
                  arma::mat& post_predictions,
                  arma::mat& post_hat_diag,
                  arma::mat& post_s2_diag,
                  std::atomic<unsigned int>* surface_active_count)
    : x_train(x_train), y_train(y_train), x_target(x_target),
      neighbor_index(neighbor_index), neighbor_distance(neighbor_distance),
      kernel(kernel), adaptive(adaptive), fixed_bandwidth(fixed_bandwidth),
      lambda(lambda), alpha(alpha), d(d), tolerance(tolerance),
      max_iterations(max_iterations), screening(screening),
      keep_coefficients(keep_coefficients),
      return_diagnostics(return_diagnostics), predictions(predictions),
      nonzero(nonzero), iterations(iterations), converged(converged),
      kkt(kkt), status(status), weight_sum(weight_sum),
      bandwidth_used(bandwidth_used), coefficients(coefficients),
      hat_diag(hat_diag), s2_diag(s2_diag),
      post_predictions(post_predictions), post_hat_diag(post_hat_diag),
      post_s2_diag(post_s2_diag),
      surface_active_count(surface_active_count) {}

  void operator()(std::size_t begin, std::size_t end) {
    const arma::uword k = neighbor_index.n_rows;
    const arma::uword p = x_train.n_cols;
    const arma::uword n_lambda = lambda.n_elem;
    const bool pure_ridge = alpha == 0.0 && d == 0.0;
    const bool ridge_only_positive = pure_ridge && lambda.min() > 0.0;
    // Preserve other estimators and optional post-OLS diagnostics exactly. Pure
    // positive Ridge without diagnostics needs neither a p-by-p Gram nor GWR.
    const bool need_legacy_gram = !ridge_only_positive || return_diagnostics;

    arma::mat x_local(k, p, arma::fill::zeros);
    arma::mat x_centered(k, p, arma::fill::zeros);
    arma::mat weighted_x(k, p, arma::fill::zeros);
    arma::mat gram;
    if (need_legacy_gram) gram.zeros(p, p);
    gwrs::RidgeSystem ridge_system;
    arma::mat score_map(p, k, arma::fill::zeros);
    arma::mat gwr_map(p, k, arma::fill::zeros);
    arma::vec y_local(k, arma::fill::zeros);
    arma::vec y_centered(k, arma::fill::zeros);
    arma::vec weights(k, arma::fill::zeros);
    arma::vec normalized_weights(k, arma::fill::zeros);
    arma::vec sqrt_weights(k, arma::fill::zeros);
    arma::vec score(p, arma::fill::zeros);
    arma::vec beta_gwr(p, arma::fill::zeros);
    arma::vec beta(p, arma::fill::zeros);
    arma::vec residual(k, arma::fill::zeros);
    arma::vec previous_kkt_score(p, arma::fill::zeros);
    arma::rowvec x_mean(p, arma::fill::zeros);
    std::vector<unsigned char> active(p, 1);

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
      y_centered = y_local - y_mean;
      x_centered = x_local;
      x_centered.each_row() -= x_mean;
      sqrt_weights = arma::sqrt(normalized_weights);
      weighted_x = x_centered;
      weighted_x.each_col() %= sqrt_weights;
      if (need_legacy_gram) {
        gram = weighted_x.t() * weighted_x;
        score = x_centered.t() * (normalized_weights % y_centered);
      }
      if (pure_ridge) ridge_system.prepare(weighted_x, sqrt_weights % y_centered);
      if (return_diagnostics) {
        score_map = x_centered.t();
        score_map.each_row() %= normalized_weights.t();
        gwr_map = solve_path_matrix(gram, score_map);
      }

      int local_status = 0;
      if (!ridge_only_positive) {
        gwrs::solve_local_system(gram, score, beta_gwr, local_status);
      } else {
        beta_gwr.zeros(); // d=0: the unpenalized center is not part of Ridge.
      }
      status[target] = local_status;
      beta.zeros();
      residual = y_centered;
      std::fill(active.begin(), active.end(), static_cast<unsigned char>(1));

      for (arma::uword path_index = 0; path_index < n_lambda; ++path_index) {
        const double current_lambda = lambda[path_index];
        const double lambda1 = current_lambda * alpha;
        const double lambda2 = current_lambda * (1.0 - alpha);

        if (screening && path_index > 0) {
          const double previous_lambda1 = lambda[path_index - 1] * alpha;
          const double threshold = std::max(0.0,
                                             2.0 * lambda1 - previous_lambda1);
          for (arma::uword column = 0; column < p; ++column) {
            active[column] = (std::abs(beta[column]) > 1e-12 ||
                              std::abs(previous_kkt_score[column]) >= threshold)
              ? 1 : 0;
          }
        } else {
          std::fill(active.begin(), active.end(),
                    static_cast<unsigned char>(1));
        }

        int total_iterations = 0;
        bool coefficient_convergence = false;
        bool needs_refit = true;

        if (current_lambda <= 0.0) {
          beta = beta_gwr;
          residual = y_centered - x_centered * beta;
          coefficient_convergence = true;
          needs_refit = false;
        } else if (pure_ridge) {
          if (!ridge_system.solve(weighted_x, current_lambda, beta)) {
            status[target] = 2; // Failed SPD solve; retain NA outputs and flag.
            converged(target, path_index) = 0;
            continue;
          }
          residual = y_centered - x_centered * beta;
          coefficient_convergence = true;
          needs_refit = false;
          total_iterations = 1; // One direct factorization/solve, not CD sweeps.
        }

        while (needs_refit && total_iterations < max_iterations) {
          needs_refit = false;
          coefficient_convergence = false;
          while (total_iterations < max_iterations) {
            double maximum_delta = 0.0;
            double maximum_beta = 0.0;
            for (arma::uword column = 0; column < p; ++column) {
              if (!active[column]) continue;
              const double curvature = gram(column, column);
              if (!std::isfinite(curvature) ||
                  curvature <= gwrs::kNumericalEpsilon) {
                if (beta[column] != 0.0) {
                  residual += x_centered.col(column) * beta[column];
                  beta[column] = 0.0;
                }
                continue;
              }

              const double partial_score =
                arma::dot(x_centered.col(column),
                          normalized_weights % residual) +
                curvature * beta[column];
              const double shifted_score = partial_score +
                lambda2 * d * beta_gwr[column];
              const double updated = gwrs::soft_threshold(shifted_score,
                                                           lambda1) /
                (curvature + lambda2);
              const double delta = updated - beta[column];
              if (delta != 0.0) {
                beta[column] = updated;
                residual -= x_centered.col(column) * delta;
                maximum_delta = std::max(maximum_delta, std::abs(delta));
              }
              maximum_beta = std::max(maximum_beta, std::abs(beta[column]));
            }
            ++total_iterations;
            if (maximum_delta <= tolerance * (1.0 + maximum_beta)) {
              coefficient_convergence = true;
              break;
            }
          }

          for (arma::uword column = 0; column < p; ++column) {
            const double current_score =
              arma::dot(x_centered.col(column),
                        normalized_weights % residual) +
              lambda2 * (d * beta_gwr[column] - beta[column]);
            double violation = 0.0;
            if (std::abs(beta[column]) <= 1e-12) {
              violation = std::max(0.0,
                                   std::abs(current_score) - lambda1);
            } else {
              violation = std::abs(current_score -
                lambda1 * (beta[column] > 0.0 ? 1.0 : -1.0));
            }
            if (!active[column] &&
                violation > tolerance * (1.0 + lambda1)) {
              active[column] = 1;
              needs_refit = true;
            }
          }
        }

        double maximum_kkt = 0.0;
        int local_nonzero = 0;
        for (arma::uword column = 0; column < p; ++column) {
          previous_kkt_score[column] =
            arma::dot(x_centered.col(column),
                      normalized_weights % residual) +
            lambda2 * (d * beta_gwr[column] - beta[column]);
          double violation = 0.0;
          if (std::abs(beta[column]) <= 1e-12) {
            violation = std::max(0.0,
                                 std::abs(previous_kkt_score[column]) -
                                   lambda1);
          } else {
            ++local_nonzero;
            violation = std::abs(previous_kkt_score[column] -
              lambda1 * (beta[column] > 0.0 ? 1.0 : -1.0));
          }
          maximum_kkt = std::max(maximum_kkt, violation);
        }

        const double intercept = y_mean - arma::dot(x_mean, beta.t());
        predictions(target, path_index) = intercept +
          arma::dot(x_target.row(target), beta.t());
        nonzero(target, path_index) = local_nonzero;
        iterations(target, path_index) = total_iterations;
        converged(target, path_index) =
          (coefficient_convergence && !needs_refit) ? 1 : 0;
        kkt(target, path_index) = maximum_kkt;

        if (return_diagnostics) {
          std::vector<arma::uword> selected;
          selected.reserve(p);
          for (arma::uword column = 0; column < p; ++column) {
            if (current_lambda <= 0.0 || alpha <= 0.0 ||
                std::abs(beta[column]) > 1e-10) {
              selected.push_back(column);
              surface_active_count[column + p * path_index].fetch_add(
                1u, std::memory_order_relaxed
              );
            }
          }
          arma::uvec selected_index(selected.size());
          for (arma::uword index = 0; index < selected_index.n_elem; ++index) {
            selected_index[index] = selected[index];
          }
          arma::mat derivative(p, k, arma::fill::zeros);
          arma::mat post_derivative(p, k, arma::fill::zeros);
          arma::vec post_beta(p, arma::fill::zeros);
          if (selected_index.n_elem > 0) {
            arma::mat system = gram.submat(selected_index, selected_index);
            system.diag() += lambda2;
            const arma::mat rhs = score_map.rows(selected_index) +
              lambda2 * d * gwr_map.rows(selected_index);
            derivative.rows(selected_index) = solve_path_matrix(system, rhs);
            const arma::mat active_gram =
              gram.submat(selected_index, selected_index);
            post_derivative.rows(selected_index) = solve_path_matrix(
              active_gram, score_map.rows(selected_index)
            );
            post_beta.elem(selected_index) = solve_path_matrix(
              active_gram, arma::mat(score.elem(selected_index))
            ).col(0);
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

        if (keep_coefficients) {
          coefficients(0, target, path_index) = intercept;
          coefficients.subcube(1, target, path_index,
                               p, target, path_index) = beta;
        }
      }
    }
  }
};

void validate_path_inputs(const arma::mat& x_train,
                          const arma::vec& y_train,
                          const arma::mat& x_target,
                          const arma::umat& neighbor_index,
                          const arma::mat& neighbor_distance,
                          const int kernel,
                          const bool adaptive,
                          const double fixed_bandwidth) {
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
}

} // anonymous namespace

// Serial, bounded-memory CV grid. Parallelism belongs to the replication
// scheduler. Each target's sufficient statistics are shared across all paths;
// paths themselves retain independent zero starts and decreasing-lambda warms.
// [[Rcpp::export]]
Rcpp::List cpp_gwr_sl_grid_predict(
    const arma::mat& x_train, const arma::vec& y_train,
    const arma::mat& x_target, const arma::umat& neighbor_index,
    const arma::mat& neighbor_distance, const int kernel,
    const bool adaptive, const double fixed_bandwidth,
    const arma::vec& lambda, const arma::vec& alpha, const arma::vec& d,
    const double tolerance, const int max_iterations,
    const bool screening, const bool keep_coefficients,
    const bool active_solve = false) {
  validate_path_inputs(x_train, y_train, x_target, neighbor_index,
    neighbor_distance, kernel, adaptive, fixed_bandwidth);
  if (x_train.n_cols == 0 || neighbor_index.n_rows == 0 ||
      lambda.n_elem == 0 || !lambda.is_finite() || arma::any(lambda < 0) ||
      alpha.n_elem == 0 || alpha.n_elem != d.n_elem ||
      !alpha.is_finite() || !d.is_finite() || arma::any(alpha <= 0) ||
      arma::any(alpha > 1) || arma::any(d < 0) || arma::any(d > 1) ||
      !std::isfinite(tolerance) || tolerance <= 0 || max_iterations < 1)
    Rcpp::stop("Invalid SL grid dimensions, penalties or solver controls");
  for (arma::uword l = 1; l < lambda.n_elem; ++l)
    if (lambda[l] > lambda[l-1]) Rcpp::stop("lambda must be decreasing");

  using Clock = std::chrono::steady_clock;
  const arma::uword p = x_train.n_cols, k = neighbor_index.n_rows;
  const arma::uword nt = x_target.n_rows, nl = lambda.n_elem, nb = alpha.n_elem;
  arma::mat predictions(nt, nl*nb, arma::fill::value(arma::datum::nan));
  arma::mat kkt(nt, nl*nb, arma::fill::value(arma::datum::nan));
  arma::imat converged(nt, nl*nb, arma::fill::zeros);
  arma::imat iterations(nt, nl*nb, arma::fill::zeros);
  arma::imat solve_attempts(nt, nl*nb, arma::fill::zeros);
  arma::imat solve_accepted(nt, nl*nb, arma::fill::zeros);
  arma::ivec status(nt, arma::fill::zeros);
  arma::vec path_seconds(nb, arma::fill::zeros);
  double preparation_seconds = 0;
  arma::cube coefficients;
  if (keep_coefficients) {
    coefficients.set_size(p+1, nt, nl*nb);
    coefficients.fill(arma::datum::nan);
  }
  arma::mat xc(k,p), weighted(k,p), gram(p,p);
  arma::vec yl(k), yc(k), weights(k), score(p), pilot(p), beta(p), gradient(p);
  arma::vec previous(p), residual(k);
  std::vector<unsigned char> active(p,1);
  for (arma::uword target = 0; target < nt; ++target) {
    Rcpp::checkUserInterrupt();
    auto start = Clock::now();
    bool valid = true;
    for (arma::uword row = 0; row < k; ++row) {
      const arma::uword encoded = neighbor_index(row,target);
      if (encoded == 0 || encoded > x_train.n_rows) { valid=false; break; }
      xc.row(row) = x_train.row(encoded-1); yl[row] = y_train[encoded-1];
    }
    if (!valid) { status[target]=4; continue; }
    double bandwidth = 0;
    const double sumw = gwrs::compute_kernel_weights(neighbor_distance,target,
      kernel,adaptive,fixed_bandwidth,weights,bandwidth);
    if (!std::isfinite(sumw) || sumw <= gwrs::kNumericalEpsilon) {
      status[target]=3; continue;
    }
    weights /= sumw;
    const arma::rowvec xm = weights.t()*xc;
    const double ym = arma::dot(weights,yl);
    xc.each_row() -= xm; yc = yl-ym;
    weighted = xc; weighted.each_col() %= arma::sqrt(weights);
    gram = weighted.t()*weighted;
    score = xc.t()*(weights % yc);
    int local_status = 0;
    gwrs::solve_local_system(gram,score,pilot,local_status);
    status[target] = local_status;
    preparation_seconds += std::chrono::duration<double>(Clock::now()-start).count();
    for (arma::uword block = 0; block < nb; ++block) {
      start = Clock::now();
      beta.zeros(); previous.zeros();
      for (arma::uword l = 0; l < nl; ++l) {
        const arma::uword col = block*nl+l;
        const double l1 = lambda[l]*alpha[block], l2 = lambda[l]*(1-alpha[block]);
        const arma::vec shift = l2*d[block]*pilot;
        gradient = score-gram*beta;
        for (arma::uword j = 0; j < p; ++j) {
          const double threshold = l > 0 ? std::max(0.0,2*l1-lambda[l-1]*alpha[block]) : 0;
          active[j] = !screening || l==0 || std::abs(beta[j])>1e-12 ||
            std::abs(previous[j])>=threshold;
        }
        int sweeps = 0;
        bool small_step = false, needs_refit = true;
        if (lambda[l] <= 0) {
          beta=pilot; small_step=true; needs_refit=false;
        }
        while (needs_refit && sweeps < max_iterations) {
          needs_refit=false; small_step=false;
          while (sweeps < max_iterations) {
            double maximum_delta=0, maximum_beta=0;
            for (arma::uword j=0; j<p; ++j) {
              if (!active[j]) continue;
              const double curvature=gram(j,j);
              if (!std::isfinite(curvature) || curvature<=gwrs::kNumericalEpsilon) {
                gradient += gram.col(j)*beta[j]; beta[j]=0; continue;
              }
              const double updated=gwrs::soft_threshold(
                gradient[j]+curvature*beta[j]+shift[j],l1)/(curvature+l2);
              const double delta=updated-beta[j];
              if (delta!=0) {
                beta[j]=updated; gradient -= gram.col(j)*delta;
                maximum_delta=std::max(maximum_delta,std::abs(delta));
              }
              maximum_beta=std::max(maximum_beta,std::abs(beta[j]));
            }
            ++sweeps;
            if (maximum_delta<=tolerance*(1+maximum_beta)) { small_step=true; break; }
            // Optional active-face solve: no additional penalty or approximate
            // inverse. Quick paths keep ordinary CD; only after 32 sweeps try
            // at geometric intervals, then retain ordinary
            // CD whenever signs, the original-data objective or full KKT fail.
            // Positive l2 makes this system strictly convex, including when
            // the local unpenalized design is rank deficient.
            if (active_solve && l2>0 && sweeps>=32 && sweeps<max_iterations &&
                (sweeps & (sweeps-1))==0) {
              const arma::uvec selected=arma::find(arma::abs(beta)>1e-12);
              if (selected.n_elem>0) {
                ++solve_attempts(target,col);
                arma::mat system=gram.submat(selected,selected);
                system.diag() += l2;
                const arma::vec signs=arma::sign(beta.elem(selected));
                const arma::vec rhs=score.elem(selected)+shift.elem(selected)-l1*signs;
                arma::vec solution;
                const bool solved=arma::solve(solution,system,rhs,
                  arma::solve_opts::likely_sympd+arma::solve_opts::no_approx);
                if (solved && solution.is_finite() && arma::all(solution % signs>1e-12)) {
                  arma::vec candidate(p,arma::fill::zeros);
                  candidate.elem(selected)=solution;
                  const arma::vec trial_residual=yc-xc*candidate;
                  const arma::vec trial_score=xc.t()*(weights % trial_residual)+shift-l2*candidate;
                  double trial_kkt=0;
                  for (arma::uword j=0;j<p;++j) {
                    const double violation=std::abs(candidate[j])<=1e-12 ?
                      std::max(0.0,std::abs(trial_score[j])-l1) :
                      std::abs(trial_score[j]-l1*(candidate[j]>0 ? 1 : -1));
                    trial_kkt=std::max(trial_kkt,violation);
                  }
                  if (trial_kkt<=tolerance*(1+l1)) {
                    const arma::vec old_residual=yc-xc*beta;
                    const arma::vec old_shift=beta-d[block]*pilot;
                    const arma::vec new_shift=candidate-d[block]*pilot;
                    const double old_objective=.5*arma::dot(weights,arma::square(old_residual))+
                      l1*arma::accu(arma::abs(beta))+.5*l2*arma::dot(old_shift,old_shift);
                    const double new_objective=.5*arma::dot(weights,arma::square(trial_residual))+
                      l1*arma::accu(arma::abs(candidate))+.5*l2*arma::dot(new_shift,new_shift);
                    const double slack=64*std::numeric_limits<double>::epsilon()*
                      std::max(1.0,std::abs(old_objective));
                    if (std::isfinite(new_objective) && new_objective<=old_objective+slack) {
                      beta=candidate;
                      gradient=xc.t()*(weights % trial_residual);
                      ++solve_accepted(target,col);
                      // The next CD sweep must still satisfy the original
                      // coefficient-change stopping rule before convergence.
                    }
                  }
                }
              }
            }
          }
          // Recompute from original weighted data, avoiding cancellation in the
          // Gram score for the screening correction and reported KKT residual.
          residual=yc-xc*beta;
          gradient=xc.t()*(weights % residual);
          for (arma::uword j=0; j<p; ++j) {
            const double s=gradient[j]+shift[j]-l2*beta[j];
            const double violation=std::abs(beta[j])<=1e-12 ?
              std::max(0.0,std::abs(s)-l1) : std::abs(s-l1*(beta[j]>0 ? 1 : -1));
            if (!active[j] && violation>tolerance*(1+l1)) {
              active[j]=1; needs_refit=true;
            }
          }
        }
        residual=yc-xc*beta;
        previous=xc.t()*(weights % residual)+shift-l2*beta;
        double worst=0;
        for (arma::uword j=0; j<p; ++j) {
          const double violation=std::abs(beta[j])<=1e-12 ?
            std::max(0.0,std::abs(previous[j])-l1) :
            std::abs(previous[j]-l1*(beta[j]>0 ? 1 : -1));
          worst=std::max(worst,violation);
        }
        const double intercept=ym-arma::dot(xm,beta.t());
        predictions(target,col)=intercept+arma::dot(x_target.row(target),beta.t());
        kkt(target,col)=worst; iterations(target,col)=sweeps;
        // Preserve the legacy convergence policy; KKT is independently exposed
        // for audits, not silently substituted for the approved stopping rule.
        converged(target,col)=(small_step && !needs_refit) ? 1 : 0;
        if (keep_coefficients) {
          coefficients(0,target,col)=intercept;
          coefficients.subcube(1,target,col,p,target,col)=beta;
        }
      }
      path_seconds[block] += std::chrono::duration<double>(Clock::now()-start).count();
    }
  }
  return Rcpp::List::create(Rcpp::_["predictions"]=predictions,
    Rcpp::_["converged"]=converged,Rcpp::_["iterations"]=iterations,
    Rcpp::_["kkt"]=kkt,Rcpp::_["status"]=status,
    Rcpp::_["active_solve_attempts"]=solve_attempts,
    Rcpp::_["active_solve_accepted"]=solve_accepted,
    Rcpp::_["coefficients"]=coefficients,Rcpp::_["path_seconds"]=path_seconds,
    Rcpp::_["preparation_seconds"]=preparation_seconds);
}

// [[Rcpp::export]]
double cpp_gwr_sl_lambda_max(const arma::mat& x,
                             const arma::vec& y,
                             const arma::umat& neighbor_index,
                             const arma::mat& neighbor_distance,
                             const int kernel,
                             const bool adaptive,
                             const double fixed_bandwidth,
                             const double alpha,
                             const int n_threads,
                             const int grain_size) {
  validate_path_inputs(x, y, x, neighbor_index, neighbor_distance, kernel,
                       adaptive, fixed_bandwidth);
  if (!std::isfinite(alpha) || alpha <= 0.0 || alpha > 1.0) {
    Rcpp::stop("alpha must be in (0, 1] for lambda_max");
  }

  arma::vec local_maximum(x.n_rows, arma::fill::zeros);
  LambdaMaxWorker worker(x, y, neighbor_index, neighbor_distance, kernel,
                         adaptive, fixed_bandwidth, alpha, local_maximum);
  if (n_threads == 1 || x.n_rows < 2) {
    worker(0, x.n_rows);
  } else {
    RcppParallel::parallelFor(
      0, x.n_rows, worker,
      static_cast<std::size_t>(std::max(1, grain_size)), n_threads);
  }
  return local_maximum.max();
}

// [[Rcpp::export]]
Rcpp::List cpp_gwr_sl_path_predict(
    const arma::mat& x_train,
    const arma::vec& y_train,
    const arma::mat& x_target,
    const arma::umat& neighbor_index,
    const arma::mat& neighbor_distance,
    const int kernel,
    const bool adaptive,
    const double fixed_bandwidth,
    const arma::vec& lambda,
    const double alpha,
    const double d,
    const double tolerance,
    const int max_iterations,
    const bool screening,
    const bool keep_coefficients,
    const bool return_diagnostics,
    const int n_threads,
    const int grain_size) {
  validate_path_inputs(x_train, y_train, x_target, neighbor_index,
                       neighbor_distance, kernel, adaptive, fixed_bandwidth);
  if (lambda.n_elem < 1 || !lambda.is_finite() || arma::any(lambda < 0.0)) {
    Rcpp::stop("lambda must be a finite nonnegative vector");
  }
  for (arma::uword index = 1; index < lambda.n_elem; ++index) {
    if (lambda[index] > lambda[index - 1]) {
      Rcpp::stop("lambda must be ordered from largest to smallest");
    }
  }
  if (!std::isfinite(alpha) || alpha < 0.0 || alpha > 1.0 ||
      !std::isfinite(d) || d < 0.0 || d > 1.0) {
    Rcpp::stop("alpha or d is outside its valid range");
  }

  const arma::uword n_target = x_target.n_rows;
  const arma::uword n_lambda = lambda.n_elem;
  const arma::uword p = x_train.n_cols;
  arma::mat predictions(n_target, n_lambda,
                        arma::fill::value(arma::datum::nan));
  arma::imat nonzero(n_target, n_lambda, arma::fill::zeros);
  arma::imat iterations(n_target, n_lambda, arma::fill::zeros);
  arma::imat converged(n_target, n_lambda, arma::fill::zeros);
  arma::mat kkt(n_target, n_lambda,
                arma::fill::value(arma::datum::nan));
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
  std::unique_ptr<std::atomic<unsigned int>[]> active_count;
  if (return_diagnostics) {
    hat_diag.set_size(n_target, n_lambda);
    s2_diag.set_size(n_target, n_lambda);
    post_predictions.set_size(n_target, n_lambda);
    post_hat_diag.set_size(n_target, n_lambda);
    post_s2_diag.set_size(n_target, n_lambda);
    hat_diag.fill(arma::datum::nan);
    s2_diag.fill(arma::datum::nan);
    post_predictions.fill(arma::datum::nan);
    post_hat_diag.fill(arma::datum::nan);
    post_s2_diag.fill(arma::datum::nan);
    active_count.reset(new std::atomic<unsigned int>[p * n_lambda]);
    for (arma::uword index = 0; index < p * n_lambda; ++index) {
      active_count[index].store(0u, std::memory_order_relaxed);
    }
  }

  GwrSlPathWorker worker(
    x_train, y_train, x_target, neighbor_index, neighbor_distance, kernel,
    adaptive, fixed_bandwidth, lambda, alpha, d, tolerance, max_iterations,
    screening, keep_coefficients, return_diagnostics, predictions, nonzero,
    iterations, converged, kkt, status, weight_sum, bandwidth_used,
    coefficients, hat_diag, s2_diag, post_predictions, post_hat_diag,
    post_s2_diag, active_count.get());
  if (n_threads == 1 || n_target < 2) {
    worker(0, n_target);
  } else {
    RcppParallel::parallelFor(
      0, n_target, worker,
      static_cast<std::size_t>(std::max(1, grain_size)), n_threads);
  }

  arma::umat surface_active_count;
  if (return_diagnostics) {
    surface_active_count.set_size(p, n_lambda);
    for (arma::uword path_index = 0; path_index < n_lambda; ++path_index) {
      for (arma::uword column = 0; column < p; ++column) {
        surface_active_count(column, path_index) =
          active_count[column + p * path_index].load(
            std::memory_order_relaxed
          );
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::_["predictions"] = predictions,
    Rcpp::_["nonzero"] = nonzero,
    Rcpp::_["iterations"] = iterations,
    Rcpp::_["converged"] = converged,
    Rcpp::_["kkt"] = kkt,
    Rcpp::_["status"] = status,
    Rcpp::_["weight_sum"] = weight_sum,
    Rcpp::_["bandwidth"] = bandwidth_used,
    Rcpp::_["coefficients"] = coefficients,
    Rcpp::_["hat_diag_path"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(hat_diag)) : R_NilValue,
    Rcpp::_["s2_diag_path"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(s2_diag)) : R_NilValue,
    Rcpp::_["post_fitted_path"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(post_predictions)) : R_NilValue,
    Rcpp::_["post_hat_diag_path"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(post_hat_diag)) : R_NilValue,
    Rcpp::_["post_s2_diag_path"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(post_s2_diag)) : R_NilValue,
    Rcpp::_["surface_active_count"] = return_diagnostics ?
      static_cast<SEXP>(Rcpp::wrap(surface_active_count)) : R_NilValue
  );
}
