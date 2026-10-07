#include "gwrs_nonconvex.h"

#include <Rcpp.h>
#include <limits>
#include <vector>

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

namespace {

// Values are stored only on neighbor edges, never in a dense n-by-n operator.
struct LocalMap {
  arma::uvec active;
  arma::mat coefficient; // k by (intercept + active slopes), original scale
  arma::vec smoother;
  int status = 0;
};

bool spd_solve(const arma::mat& system, const arma::mat& rhs, arma::mat& answer) {
  arma::mat upper;
  if (!arma::chol(upper, system)) return false;
  const bool ok = arma::solve(answer, system, rhs,
    arma::solve_opts::likely_sympd + arma::solve_opts::no_approx);
  return ok && answer.is_finite();
}

arma::uvec active_set(const arma::mat& beta, const arma::vec& scale,
                     arma::uword target, double lambda, double alpha, int mode) {
  const arma::uword p = scale.n_elem;
  if (lambda <= 0 || (mode == 0 && alpha <= 0)) {
    return arma::regspace<arma::uvec>(0, p - 1);
  }
  std::vector<arma::uword> selected;
  for (arma::uword j = 0; j < p; ++j) {
    // The fitted active set and penalty are defined on the working scale.
    if (std::abs(beta(j + 1, target) * scale[j]) > 1e-10) selected.push_back(j);
  }
  return arma::uvec(selected);
}

struct MapWorker : public RcppParallel::Worker {
  const arma::mat& x;
  const arma::mat& beta;
  const arma::umat& index;
  const arma::mat& distance;
  const arma::vec& center;
  const arma::vec& scale;
  int kernel, mode;
  bool adaptive;
  double bandwidth, lambda, alpha, d, gamma;
  std::vector<LocalMap>& maps;

  MapWorker(const arma::mat& x, const arma::mat& beta, const arma::umat& index,
            const arma::mat& distance, const arma::vec& center, const arma::vec& scale,
            int kernel, bool adaptive, double bandwidth, double lambda,
            double alpha, double d, int mode, double gamma, std::vector<LocalMap>& maps)
    : x(x), beta(beta), index(index), distance(distance), center(center), scale(scale),
      kernel(kernel), mode(mode), adaptive(adaptive), bandwidth(bandwidth),
      lambda(lambda), alpha(alpha), d(d), gamma(gamma), maps(maps) {}

  void operator()(std::size_t begin, std::size_t end) {
    const arma::uword k = index.n_rows;
    for (std::size_t i = begin; i < end; ++i) {
      LocalMap& out = maps[i];
      out.active = active_set(beta, scale, i, lambda, alpha, mode);
      const arma::uword a = out.active.n_elem;
      arma::vec weight(k);
      double used = 0;
      const double mass = gwrs::compute_kernel_weights(
        distance, i, kernel, adaptive, bandwidth, weight, used);
      if (!std::isfinite(mass) || mass <= gwrs::kNumericalEpsilon) {
        out.status = 1;
        continue;
      }
      weight /= mass;
      arma::mat local(k, a);
      for (arma::uword j = 0; j < a; ++j)
        for (arma::uword row = 0; row < k; ++row)
          local(row, j) = x(index(row, i) - 1, out.active[j]);
      arma::rowvec mean = weight.t() * local;
      local.each_row() -= mean;
      arma::mat rhs = local.t();
      rhs.each_row() %= weight.t();
      arma::mat derivative(a, k, arma::fill::zeros);
      if (a > 0) {
        arma::mat system = rhs * local;
        if (mode == 0) {
          const double lambda2 = lambda * (1 - alpha);
          system.diag() += lambda2;
          if (lambda2 * d != 0) {
            // Centered shrinkage uses the existing full-GWR pilot.
            arma::mat full(k, x.n_cols);
            for (arma::uword row = 0; row < k; ++row)
              full.row(row) = x.row(index(row, i) - 1);
            const arma::rowvec full_mean = weight.t() * full;
            full.each_row() -= full_mean;
            arma::mat pilot_rhs = full.t();
            pilot_rhs.each_row() %= weight.t();
            arma::mat pilot;
            if (!spd_solve(pilot_rhs * full, pilot_rhs, pilot)) {
              out.status = 2;
              continue;
            }
            rhs += lambda2 * d * pilot.rows(out.active);
          }
        } else {
          for (arma::uword j = 0; j < a; ++j) {
            const double working = std::abs(beta(out.active[j] + 1, i) * scale[out.active[j]]);
            if (gwrs::nonconvex_penalty_knot(working, lambda, gamma, mode - 1)) {
              out.status = 3;
              break;
            }
            system(j, j) += gwrs::nonconvex_penalty_second_derivative(
              working, lambda, gamma, mode - 1);
          }
          if (out.status) continue;
        }
        if (!spd_solve(system, rhs, derivative)) {
          out.status = 2;
          continue;
        }
      }
      out.coefficient.set_size(k, a + 1);
      out.coefficient.col(0) = weight - derivative.t() * mean.t();
      arma::vec focus(a);
      for (arma::uword j = 0; j < a; ++j) {
        const arma::uword column = out.active[j];
        out.coefficient.col(j + 1) = derivative.row(j).t() / scale[column];
        out.coefficient.col(0) -= center[column] * out.coefficient.col(j + 1);
        focus[j] = x(i, column) - mean[j];
      }
      out.smoother = weight + derivative.t() * focus;
      if (!out.coefficient.is_finite() || !out.smoother.is_finite()) out.status = 4;
    }
  }
};

// Right multiplication by an implicit unit-vector block.
arma::mat extract_columns(const arma::mat& values, const arma::umat& index,
                          arma::uword first, arma::uword width) {
  arma::mat answer(index.n_cols, width, arma::fill::zeros);
  for (arma::uword i = 0; i < index.n_cols; ++i) {
    for (arma::uword row = 0; row < index.n_rows; ++row) {
      const arma::uword source = index(row, i) - 1;
      if (source >= first && source < first + width)
        answer(i, source - first) += values(row, i);
    }
  }
  return answer;
}

arma::mat transpose_apply(const arma::mat& values, const arma::umat& index,
                          const arma::mat& rhs) {
  arma::mat answer(index.n_cols, rhs.n_cols, arma::fill::zeros);
  for (arma::uword column = 0; column < rhs.n_cols; ++column) {
    double* dest = answer.colptr(column);
    const double* input = rhs.colptr(column);
    for (arma::uword i = 0; i < index.n_cols; ++i)
      for (arma::uword row = 0; row < index.n_rows; ++row)
        dest[index(row, i) - 1] += values(row, i) * input[i];
  }
  return answer;
}

struct MomentWorker : public RcppParallel::Worker {
  const arma::mat& values;
  const arma::umat& index;
  arma::uword block;
  bool residual;
  arma::vec& diagonal;
  arma::vec& squares;

  MomentWorker(const arma::mat& values, const arma::umat& index, arma::uword block,
               bool residual, arma::vec& diagonal, arma::vec& squares)
    : values(values), index(index), block(block), residual(residual),
      diagonal(diagonal), squares(squares) {}

  void operator()(std::size_t begin, std::size_t end) {
    const arma::uword n = index.n_cols;
    for (std::size_t task = begin; task < end; ++task) {
      const arma::uword first = task * block, width = std::min(block, n - first);
      arma::mat right = extract_columns(values, index, first, width);
      arma::mat result;
      if (residual) {
        right *= -1;
        for (arma::uword c = 0; c < width; ++c) right(first + c, c) += 1;
        // Q E = (I-S)' (I-S) E, including all off-diagonal contributions.
        result = right - transpose_apply(values, index, right);
      } else {
        right.each_row() -= arma::mean(right, 0);
        // B_j E = C_j' J C_j E / n.
        result = transpose_apply(values, index, right) / static_cast<double>(n);
      }
      for (arma::uword c = 0; c < width; ++c) {
        diagonal[first + c] = result(first + c, c);
        squares[first + c] = arma::dot(result.col(c), result.col(c));
      }
    }
  }
};

void moments(const arma::mat& values, const arma::umat& index, arma::uword block,
             bool residual, int threads, arma::vec& diagonal, double& second) {
  const arma::uword n = index.n_cols, tasks = (n + block - 1) / block;
  arma::vec squares(n, arma::fill::zeros);
  diagonal.zeros(n);
  MomentWorker worker(values, index, block, residual, diagonal, squares);
  // Bounded batches allow interruption between parallel calls.
  const arma::uword batch = std::max(1, threads);
  for (arma::uword start = 0; start < tasks; start += batch) {
    const arma::uword stop = std::min(tasks, start + batch);
    Rcpp::checkUserInterrupt();
    if (threads == 1) worker(start, stop);
    else RcppParallel::parallelFor(start, stop, worker, 1, threads);
  }
  second = arma::accu(squares);
}

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List cpp_gwr_f_test_components(
    const arma::mat& x, const arma::mat& coefficients,
    const arma::umat& neighbor_index, const arma::mat& neighbor_distance,
    const arma::vec& x_center, const arma::vec& x_scale,
    const int kernel, const bool adaptive, const double fixed_bandwidth,
    const double lambda, const double alpha, const double d,
    const int penalty_mode, const double gamma,
    const int n_threads, const int grain_size,
    const arma::mat& global_basis, const int block_size, const double memory_limit_mb) {
  const arma::uword n = x.n_rows, p = x.n_cols, k = neighbor_index.n_rows;
  if (n < 2 || p < 1 || k < 1 || k > n ||
      coefficients.n_rows != p + 1 || coefficients.n_cols != n ||
      neighbor_index.n_cols != n || neighbor_distance.n_cols != n ||
      neighbor_distance.n_rows != k || x_center.n_elem != p || x_scale.n_elem != p ||
      global_basis.n_rows != n || global_basis.n_cols >= n)
    Rcpp::stop("incompatible F-test input dimensions");
  if (!x.is_finite() || !coefficients.is_finite() || !x_center.is_finite() ||
      !x_scale.is_finite() || arma::any(x_scale <= 0) || !global_basis.is_finite() ||
      !neighbor_distance.is_finite() || arma::any(arma::vectorise(neighbor_distance) < 0) ||
      !std::isfinite(lambda) || lambda < 0 || !std::isfinite(alpha) || alpha < 0 || alpha > 1 ||
      !std::isfinite(d) || kernel < 0 || kernel > 4 ||
      (!adaptive && (!std::isfinite(fixed_bandwidth) || fixed_bandwidth <= 0)) ||
      n_threads < 1 || grain_size < 1 || block_size < 1 ||
      !std::isfinite(memory_limit_mb) || memory_limit_mb <= 0)
    Rcpp::stop("invalid F-test inputs or resource limits");
  if (penalty_mode < 0 || penalty_mode > 2 ||
      (penalty_mode > 0 && (!std::isfinite(gamma) ||
       (penalty_mode == 1 && gamma <= 1) || (penalty_mode == 2 && gamma <= 2))))
    Rcpp::stop("invalid F-test penalty mode or gamma");
  const arma::uword block = std::min<arma::uword>(block_size, n - 1);
  const int threads = std::min<int>(n_threads, static_cast<int>(std::min<arma::uword>(
    n, std::numeric_limits<int>::max())));
  // Long double arithmetic prevents overflow before any large allocation.
  long double active_total = 0, max_active = 0;
  for (arma::uword i = 0; i < n; ++i) {
    arma::uword count = 0;
    for (arma::uword j = 0; j < p; ++j)
      if (lambda <= 0 || (penalty_mode == 0 && alpha <= 0) ||
          std::abs(coefficients(j + 1, i) * x_scale[j]) > 1e-10) ++count;
    active_total += count;
    max_active = std::max(max_active, static_cast<long double>(count));
  }
  const long double nn = n, pp = p, kk = k;
  const long double a = (penalty_mode == 0 && lambda * (1 - alpha) * d != 0) ? pp : max_active;
  const long double cache = 8 * kk * (active_total + 2 * nn) + 8 * active_total +
    nn * (sizeof(LocalMap) + 64);
  const long double shared = 8 * (2 * nn * kk + 6 * nn +
    8 * nn * (pp + 1) + 2 * nn * global_basis.n_cols + 8 * pp);
  const long double work = threads * 8 * (4 * nn * block + 10 * kk * (a + 1) +
    6 * (a + 1) * (a + 1));
  const double estimated_mb = static_cast<double>((cache + shared + work) / (1024 * 1024));
  if (!std::isfinite(estimated_mb) || estimated_mb > memory_limit_mb)
    Rcpp::stop("Exact F-test workspace estimate %.2f MiB exceeds memory_limit_mb %.2f; "
               "reduce block_size/threads or explicitly raise the limit. No approximation was used.",
               estimated_mb, memory_limit_mb);
  // Duplicate edges would obscure the per-edge cache contract.
  arma::uvec seen(n, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i)
    for (arma::uword row = 0; row < k; ++row) {
      const arma::uword encoded = neighbor_index(row, i);
      if (encoded == 0 || encoded > n || seen[encoded - 1] == i + 1)
        Rcpp::stop("invalid or duplicate F-test neighbor index");
      seen[encoded - 1] = i + 1;
    }
  std::vector<LocalMap> maps(n);
  MapWorker worker(x, coefficients, neighbor_index, neighbor_distance, x_center, x_scale,
    kernel, adaptive, fixed_bandwidth, lambda, alpha, d, penalty_mode, gamma, maps);
  if (threads == 1) worker(0, n);
  else RcppParallel::parallelFor(0, n, worker, grain_size, threads);
  for (arma::uword i = 0; i < n; ++i) {
    if (maps[i].status)
      Rcpp::stop("F-test local operator unavailable at location %d (status %d: "
                 "1 invalid weights; 2 singular/non-positive Hessian or pilot; "
                 "3 penalty knot; 4 non-finite map). No pseudoinverse was substituted.",
                 static_cast<int>(i + 1), maps[i].status);
  }
  arma::mat values(k, n);
  double tr_s = 0, tr_sts = 0, projection = 0;
  for (arma::uword i = 0; i < n; ++i) {
    values.col(i) = maps[i].smoother;
    tr_sts += arma::dot(values.col(i), values.col(i));
    for (arma::uword row = 0; row < k; ++row)
      if (neighbor_index(row, i) - 1 == i) tr_s += values(row, i);
    // tr(QH) = ||(I-S) U||_F^2 for the global QR basis U.
    for (arma::uword c = 0; c < global_basis.n_cols; ++c) {
      double residual = global_basis(i, c);
      for (arma::uword row = 0; row < k; ++row)
        residual -= values(row, i) * global_basis(neighbor_index(row, i) - 1, c);
      projection += residual * residual;
    }
  }
  arma::vec q_diag;
  double delta2 = 0;
  moments(values, neighbor_index, block, true, threads, q_diag, delta2);
  const double delta1 = arma::accu(q_diag);
  arma::vec gamma1(p + 1, arma::fill::zeros), gamma2(p + 1, arma::fill::zeros);
  arma::vec diagonal;
  for (arma::uword j = 0; j <= p; ++j) {
    Rcpp::checkUserInterrupt();
    values.zeros();
    for (arma::uword i = 0; i < n; ++i) {
      if (j == 0) values.col(i) = maps[i].coefficient.col(0);
      else {
        const auto& active = maps[i].active;
        const auto found = std::lower_bound(active.begin(), active.end(), j - 1);
        if (found != active.end() && *found == j - 1)
          values.col(i) = maps[i].coefficient.col(1 + (found - active.begin()));
      }
    }
    const double energy = arma::accu(arma::square(values)) / n;
    if (energy == 0) continue;
    moments(values, neighbor_index, block, false, threads, diagonal, gamma2[j]);
    gamma1[j] = arma::accu(diagonal);
    if (gamma1[j] <= 64 * std::numeric_limits<double>::epsilon() * energy) {
      gamma1[j] = 0;
      gamma2[j] = 0;
    }
  }
  return Rcpp::List::create(
    Rcpp::_["q_diag"] = q_diag, Rcpp::_["delta1"] = delta1,
    Rcpp::_["delta2"] = delta2, Rcpp::_["gamma1"] = gamma1,
    Rcpp::_["gamma2"] = gamma2, Rcpp::_["trS"] = tr_s, Rcpp::_["trStS"] = tr_sts,
    Rcpp::_["trQH"] = projection, Rcpp::_["block_size"] = block,
    Rcpp::_["threads"] = threads, Rcpp::_["workspace_estimate_mb"] = estimated_mb,
    Rcpp::_["moment_method"] = "exact_neighbor_block");
}
