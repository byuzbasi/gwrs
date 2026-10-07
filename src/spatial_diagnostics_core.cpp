#include "gwrs_common.h"

#include <Rcpp.h>

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <random>

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

namespace {

arma::mat symmetric_inverse(const arma::mat& matrix, bool& full_rank) {
  arma::vec eigenvalues;
  full_rank = false;
  if (arma::eig_sym(eigenvalues, matrix) && eigenvalues.is_finite()) {
    const double maximum = eigenvalues.max();
    const double tolerance = std::max(1.0, maximum) * 1e-10;
    full_rank = arma::all(eigenvalues > tolerance);
  }
  arma::mat inverse;
  if (full_rank && arma::inv_sympd(inverse, matrix) && inverse.is_finite()) {
    return inverse;
  }
  const double scale = std::max(1.0, arma::norm(matrix, "inf"));
  return arma::pinv(matrix, scale * 1e-10);
}

struct CollinearityWorker : public RcppParallel::Worker {
  const arma::mat& x;
  const arma::umat& neighbor_index;
  const arma::mat& neighbor_distance;
  const int kernel;
  const bool adaptive;
  const double fixed_bandwidth;
  const bool center;
  const bool scale;
  const bool return_pairwise;

  arma::vec& condition_number;
  arma::mat& vif;
  arma::vec& max_abs_correlation;
  arma::vec& effective_sample_size;
  arma::ivec& positive_weights;
  arma::cube& correlations;

  CollinearityWorker(
      const arma::mat& x,
      const arma::umat& neighbor_index,
      const arma::mat& neighbor_distance,
      const int kernel,
      const bool adaptive,
      const double fixed_bandwidth,
      const bool center,
      const bool scale,
      const bool return_pairwise,
      arma::vec& condition_number,
      arma::mat& vif,
      arma::vec& max_abs_correlation,
      arma::vec& effective_sample_size,
      arma::ivec& positive_weights,
      arma::cube& correlations)
    : x(x), neighbor_index(neighbor_index),
      neighbor_distance(neighbor_distance), kernel(kernel),
      adaptive(adaptive), fixed_bandwidth(fixed_bandwidth), center(center),
      scale(scale), return_pairwise(return_pairwise),
      condition_number(condition_number), vif(vif),
      max_abs_correlation(max_abs_correlation),
      effective_sample_size(effective_sample_size),
      positive_weights(positive_weights), correlations(correlations) {}

  void operator()(std::size_t begin, std::size_t end) {
    const arma::uword k = neighbor_index.n_rows;
    const arma::uword p = x.n_cols;
    arma::mat local_x(k, p, arma::fill::zeros);
    arma::mat z(k, p, arma::fill::zeros);
    arma::vec weights(k, arma::fill::zeros);
    arma::vec normalized_weights(k, arma::fill::zeros);
    arma::rowvec means(p, arma::fill::zeros);
    arma::rowvec scales(p, arma::fill::ones);

    for (std::size_t target = begin; target < end; ++target) {
      bool valid = true;
      for (arma::uword row = 0; row < k; ++row) {
        const arma::uword encoded = neighbor_index(row, target);
        if (encoded == 0 || encoded > x.n_rows) {
          valid = false;
          break;
        }
        local_x.row(row) = x.row(encoded - 1);
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
      effective_sample_size[target] = 1.0 /
        std::max(arma::dot(normalized_weights, normalized_weights),
                 gwrs::kNumericalEpsilon);
      positive_weights[target] = static_cast<int>(
        arma::accu(weights > gwrs::kNumericalEpsilon)
      );

      means.zeros();
      if (center) means = normalized_weights.t() * local_x;
      z = local_x;
      z.each_row() -= means;
      scales.ones();
      if (scale) {
        for (arma::uword column = 0; column < p; ++column) {
          const double variance = arma::dot(
            normalized_weights, arma::square(z.col(column))
          );
          scales[column] = std::sqrt(std::max(variance, 0.0));
          if (!std::isfinite(scales[column]) ||
              scales[column] <= gwrs::kNumericalEpsilon) {
            scales[column] = arma::datum::nan;
            z.col(column).fill(arma::datum::nan);
          } else {
            z.col(column) /= scales[column];
          }
        }
      }
      if (!z.is_finite()) {
        condition_number[target] = arma::datum::inf;
        vif.row(target).fill(arma::datum::inf);
        continue;
      }

      arma::mat weighted_z = z;
      weighted_z.each_col() %= normalized_weights;
      arma::mat correlation = z.t() * weighted_z;
      correlation = 0.5 * (correlation + correlation.t());
      if (!scale) {
        if (arma::any(correlation.diag() <=
                      gwrs::kNumericalEpsilon * gwrs::kNumericalEpsilon)) {
          condition_number[target] = arma::datum::inf;
          vif.row(target).fill(arma::datum::inf);
          continue;
        }
        const arma::vec diagonal = arma::sqrt(
          arma::clamp(correlation.diag(), gwrs::kNumericalEpsilon,
                      arma::datum::inf)
        );
        correlation.each_col() /= diagonal;
        correlation.each_row() /= diagonal.t();
      }
      correlation.diag().ones();

      arma::vec eigenvalues;
      bool full_rank = false;
      if (arma::eig_sym(eigenvalues, correlation) &&
          eigenvalues.is_finite()) {
        const double maximum = eigenvalues.max();
        const double tolerance = std::max(1.0, maximum) * 1e-10;
        full_rank = arma::all(eigenvalues > tolerance);
        condition_number[target] = full_rank ?
          std::sqrt(maximum / eigenvalues.min()) : arma::datum::inf;
      }
      bool inverse_full_rank = false;
      const arma::mat inverse = symmetric_inverse(
        correlation, inverse_full_rank
      );
      if (full_rank && inverse_full_rank && inverse.is_finite()) {
        vif.row(target) = inverse.diag().t();
      } else {
        vif.row(target).fill(arma::datum::inf);
      }
      double maximum_correlation = 0.0;
      for (arma::uword row = 0; row < p; ++row) {
        for (arma::uword column = row + 1; column < p; ++column) {
          maximum_correlation = std::max(
            maximum_correlation, std::abs(correlation(row, column))
          );
        }
      }
      max_abs_correlation[target] = maximum_correlation;
      if (return_pairwise) correlations.slice(target) = correlation;
    }
  }
};

struct LocalResidualWorker : public RcppParallel::Worker {
  const arma::vec& y;
  const arma::vec& fitted;
  const arma::umat& neighbor_index;
  const arma::mat& neighbor_distance;
  const int kernel;
  const bool adaptive;
  const double fixed_bandwidth;
  arma::vec& local_r2;
  arma::vec& local_rmse;
  arma::vec& local_mae;
  arma::vec& effective_sample_size;
  arma::ivec& positive_weights;

  LocalResidualWorker(
      const arma::vec& y,
      const arma::vec& fitted,
      const arma::umat& neighbor_index,
      const arma::mat& neighbor_distance,
      const int kernel,
      const bool adaptive,
      const double fixed_bandwidth,
      arma::vec& local_r2,
      arma::vec& local_rmse,
      arma::vec& local_mae,
      arma::vec& effective_sample_size,
      arma::ivec& positive_weights)
    : y(y), fitted(fitted), neighbor_index(neighbor_index),
      neighbor_distance(neighbor_distance), kernel(kernel),
      adaptive(adaptive), fixed_bandwidth(fixed_bandwidth),
      local_r2(local_r2), local_rmse(local_rmse), local_mae(local_mae),
      effective_sample_size(effective_sample_size),
      positive_weights(positive_weights) {}

  void operator()(std::size_t begin, std::size_t end) {
    const arma::uword k = neighbor_index.n_rows;
    arma::vec weights(k, arma::fill::zeros);
    arma::vec normalized(k, arma::fill::zeros);
    arma::vec local_y(k, arma::fill::zeros);
    arma::vec local_fitted(k, arma::fill::zeros);
    for (std::size_t target = begin; target < end; ++target) {
      bool valid = true;
      for (arma::uword row = 0; row < k; ++row) {
        const arma::uword encoded = neighbor_index(row, target);
        if (encoded == 0 || encoded > y.n_elem) {
          valid = false;
          break;
        }
        local_y[row] = y[encoded - 1];
        local_fitted[row] = fitted[encoded - 1];
      }
      if (!valid) continue;
      double bandwidth_used = 0.0;
      const double weight_sum = gwrs::compute_kernel_weights(
        neighbor_distance, target, kernel, adaptive, fixed_bandwidth,
        weights, bandwidth_used
      );
      if (!std::isfinite(weight_sum) ||
          weight_sum <= gwrs::kNumericalEpsilon) continue;
      normalized = weights / weight_sum;
      const double mean_y = arma::dot(normalized, local_y);
      const arma::vec residual = local_y - local_fitted;
      const double tss = arma::dot(
        normalized, arma::square(local_y - mean_y)
      );
      const double rss = arma::dot(normalized, arma::square(residual));
      if (tss > gwrs::kNumericalEpsilon) local_r2[target] = 1.0 - rss / tss;
      local_rmse[target] = std::sqrt(std::max(0.0, rss));
      local_mae[target] = arma::dot(normalized, arma::abs(residual));
      effective_sample_size[target] = 1.0 /
        std::max(arma::dot(normalized, normalized),
                 gwrs::kNumericalEpsilon);
      positive_weights[target] = static_cast<int>(
        arma::accu(weights > gwrs::kNumericalEpsilon)
      );
    }
  }
};

std::uint64_t mix_seed(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

void permuted_moran_values(const arma::vec& z,
                           const arma::umat& neighbor_index,
                           const arma::mat& weights,
                           const double m2,
                           const std::uint64_t seed,
                           const arma::uword permutation,
                           double& global_value,
                           arma::vec& local_values) {
  arma::vec permuted = z;
  std::mt19937_64 generator(mix_seed(seed + permutation));
  std::shuffle(permuted.begin(), permuted.end(), generator);
  const arma::uword n = z.n_elem;
  const arma::uword k = neighbor_index.n_rows;
  double numerator = 0.0;
  for (arma::uword target = 0; target < n; ++target) {
    double lag = 0.0;
    for (arma::uword row = 0; row < k; ++row) {
      const arma::uword source = neighbor_index(row, target) - 1;
      lag += weights(row, target) * permuted[source];
    }
    local_values[target] = permuted[target] * lag / m2;
    numerator += permuted[target] * lag;
  }
  global_value = numerator / arma::dot(permuted, permuted);
}

constexpr std::size_t kDeterministicMoranMomentBlocks = 64;

struct MoranMomentBlockWorker : public RcppParallel::Worker {
  const arma::vec& z;
  const arma::umat& neighbor_index;
  const arma::mat& weights;
  const double m2;
  const std::uint64_t seed;
  const std::size_t permutations;
  const std::size_t blocks;
  arma::vec& global_permuted;
  arma::mat& block_sum;
  arma::mat& block_sumsq;

  MoranMomentBlockWorker(const arma::vec& z,
                         const arma::umat& neighbor_index,
                         const arma::mat& weights,
                         const double m2,
                         const std::uint64_t seed,
                         const std::size_t permutations,
                         const std::size_t blocks,
                         arma::vec& global_permuted,
                         arma::mat& block_sum,
                         arma::mat& block_sumsq)
    : z(z), neighbor_index(neighbor_index), weights(weights), m2(m2),
      seed(seed), permutations(permutations), blocks(blocks),
      global_permuted(global_permuted), block_sum(block_sum),
      block_sumsq(block_sumsq) {}

  void operator()(std::size_t begin, std::size_t end) {
    arma::vec local_values(z.n_elem, arma::fill::zeros);
    arma::vec local_sum(z.n_elem, arma::fill::zeros);
    arma::vec local_sumsq(z.n_elem, arma::fill::zeros);
    for (std::size_t block = begin; block < end; ++block) {
      local_sum.zeros();
      local_sumsq.zeros();
      const std::size_t permutation_begin = block * permutations / blocks;
      const std::size_t permutation_end =
        (block + 1) * permutations / blocks;
      for (std::size_t permutation = permutation_begin;
           permutation < permutation_end; ++permutation) {
        double global_value = 0.0;
        permuted_moran_values(
          z, neighbor_index, weights, m2, seed,
          static_cast<arma::uword>(permutation), global_value, local_values
        );
        global_permuted[permutation] = global_value;
        local_sum += local_values;
        local_sumsq += arma::square(local_values);
      }
      block_sum.col(block) = local_sum;
      block_sumsq.col(block) = local_sumsq;
    }
  }
};

struct MoranCountWorker : public RcppParallel::Worker {
  const arma::vec& z;
  const arma::umat& neighbor_index;
  const arma::mat& weights;
  const double m2;
  const std::uint64_t seed;
  const arma::vec& observed;
  const arma::vec& center;
  arma::uvec exceedance;

  MoranCountWorker(const arma::vec& z,
                   const arma::umat& neighbor_index,
                   const arma::mat& weights,
                   const double m2,
                   const std::uint64_t seed,
                   const arma::vec& observed,
                   const arma::vec& center)
    : z(z), neighbor_index(neighbor_index), weights(weights), m2(m2),
      seed(seed), observed(observed), center(center),
      exceedance(z.n_elem, arma::fill::zeros) {}

  MoranCountWorker(MoranCountWorker& other, RcppParallel::Split)
    : z(other.z), neighbor_index(other.neighbor_index), weights(other.weights),
      m2(other.m2), seed(other.seed), observed(other.observed),
      center(other.center), exceedance(z.n_elem, arma::fill::zeros) {}

  void operator()(std::size_t begin, std::size_t end) {
    arma::vec local_values(z.n_elem, arma::fill::zeros);
    for (std::size_t permutation = begin; permutation < end; ++permutation) {
      double global_value = 0.0;
      permuted_moran_values(
        z, neighbor_index, weights, m2, seed,
        static_cast<arma::uword>(permutation), global_value, local_values
      );
      for (arma::uword target = 0; target < z.n_elem; ++target) {
        if (std::abs(local_values[target] - center[target]) >=
            std::abs(observed[target] - center[target])) {
          ++exceedance[target];
        }
      }
    }
  }

  void join(const MoranCountWorker& rhs) {
    exceedance += rhs.exceedance;
  }
};

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List cpp_gwr_global_collinearity(const arma::mat& x) {
  if (x.n_rows < 2 || x.n_cols < 1 || !x.is_finite()) {
    Rcpp::stop("global collinearity requires at least two finite rows");
  }
  const arma::uword p = x.n_cols;
  arma::mat z = x;
  z.each_row() -= arma::mean(x, 0);
  const arma::rowvec scales = arma::sqrt(arma::mean(arma::square(z), 0));
  arma::uvec invalid = arma::find(scales.t() <= gwrs::kNumericalEpsilon);
  arma::mat correlation(p, p, arma::fill::value(arma::datum::nan));
  arma::vec eigenvalues(p, arma::fill::value(arma::datum::nan));
  arma::vec vif(p, arma::fill::value(arma::datum::inf));
  double cn = arma::datum::inf;
  double max_correlation = arma::datum::nan;
  int rank = NA_INTEGER;
  bool full_rank = false;
  if (scales.is_finite() && invalid.is_empty()) {
    z.each_row() /= scales;
    correlation = z.t() * z / static_cast<double>(x.n_rows);
    correlation = 0.5 * (correlation + correlation.t());
    correlation.diag().ones();
    if (arma::eig_sym(eigenvalues, correlation) && eigenvalues.is_finite()) {
      const double maximum = eigenvalues.max();
      const double tolerance = std::max(1.0, maximum) * 1e-10;
      rank = static_cast<int>(arma::accu(eigenvalues > tolerance));
      full_rank = rank == static_cast<int>(p);
      if (full_rank) {
        cn = std::sqrt(maximum / eigenvalues.min());
        arma::mat inverse;
        if (arma::inv_sympd(inverse, correlation) && inverse.is_finite()) {
          vif = inverse.diag();
        }
      }
    }
    max_correlation = 0;
    for (arma::uword a = 0; a < p; ++a) {
      for (arma::uword b = a + 1; b < p; ++b) {
        max_correlation = std::max(max_correlation, std::abs(correlation(a, b)));
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::_["condition_number"] = cn, Rcpp::_["vif"] = vif,
    Rcpp::_["correlation"] = correlation, Rcpp::_["eigenvalues"] = eigenvalues,
    Rcpp::_["rank"] = rank, Rcpp::_["full_rank"] = full_rank,
    Rcpp::_["constant_columns"] = invalid + 1,
    Rcpp::_["max_abs_correlation"] = max_correlation
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_gwr_local_collinearity(
    const arma::mat& x,
    const arma::umat& neighbor_index,
    const arma::mat& neighbor_distance,
    const int kernel,
    const bool adaptive,
    const double fixed_bandwidth,
    const bool center,
    const bool scale,
    const bool return_pairwise,
    const int n_threads,
    const int grain_size) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (neighbor_index.n_cols != n || neighbor_distance.n_cols != n ||
      neighbor_index.n_rows != neighbor_distance.n_rows) {
    Rcpp::stop("neighbor dimensions are incompatible with x");
  }
  arma::vec condition_number(n, arma::fill::value(arma::datum::nan));
  arma::mat vif(n, p, arma::fill::value(arma::datum::nan));
  arma::vec max_abs_correlation(n, arma::fill::value(arma::datum::nan));
  arma::vec effective_sample_size(n, arma::fill::value(arma::datum::nan));
  arma::ivec positive_weights(n, arma::fill::zeros);
  arma::cube correlations;
  if (return_pairwise) {
    correlations.set_size(p, p, n);
    correlations.fill(arma::datum::nan);
  }
  CollinearityWorker worker(
    x, neighbor_index, neighbor_distance, kernel, adaptive, fixed_bandwidth,
    center, scale, return_pairwise, condition_number, vif,
    max_abs_correlation, effective_sample_size, positive_weights, correlations
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
    Rcpp::_ ["condition_number"] = condition_number,
    Rcpp::_ ["vif"] = vif,
    Rcpp::_ ["max_abs_correlation"] = max_abs_correlation,
    Rcpp::_ ["effective_sample_size"] = effective_sample_size,
    Rcpp::_ ["positive_weights"] = positive_weights,
    Rcpp::_ ["correlations"] = return_pairwise ?
      static_cast<SEXP>(Rcpp::wrap(correlations)) : R_NilValue
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_gwr_local_residual_diagnostics(
    const arma::vec& y,
    const arma::vec& fitted,
    const arma::umat& neighbor_index,
    const arma::mat& neighbor_distance,
    const int kernel,
    const bool adaptive,
    const double fixed_bandwidth,
    const int n_threads,
    const int grain_size) {
  const arma::uword n = y.n_elem;
  if (fitted.n_elem != n || neighbor_index.n_cols != n ||
      neighbor_distance.n_cols != n ||
      neighbor_index.n_rows != neighbor_distance.n_rows) {
    Rcpp::stop("incompatible local residual diagnostic dimensions");
  }
  arma::vec local_r2(n, arma::fill::value(arma::datum::nan));
  arma::vec local_rmse(n, arma::fill::value(arma::datum::nan));
  arma::vec local_mae(n, arma::fill::value(arma::datum::nan));
  arma::vec effective_sample_size(n, arma::fill::value(arma::datum::nan));
  arma::ivec positive_weights(n, arma::fill::zeros);
  LocalResidualWorker worker(
    y, fitted, neighbor_index, neighbor_distance, kernel, adaptive,
    fixed_bandwidth, local_r2, local_rmse, local_mae,
    effective_sample_size, positive_weights
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
    Rcpp::_ ["local_r2"] = local_r2,
    Rcpp::_ ["local_rmse"] = local_rmse,
    Rcpp::_ ["local_mae"] = local_mae,
    Rcpp::_ ["effective_sample_size"] = effective_sample_size,
    Rcpp::_ ["positive_weights"] = positive_weights
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_gwr_moran(
    const arma::vec& residuals,
    const arma::umat& neighbor_index,
    const arma::mat& neighbor_distance,
    const int kernel,
    const bool adaptive,
    const double fixed_bandwidth,
    const bool row_standardize,
    const int permutations,
    const double seed,
    const int n_threads,
    const int grain_size) {
  const arma::uword n = residuals.n_elem;
  const arma::uword k = neighbor_index.n_rows;
  if (n < 2 || neighbor_index.n_cols != n ||
      neighbor_distance.n_cols != n || neighbor_distance.n_rows != k) {
    Rcpp::stop("invalid residual or neighbor dimensions");
  }
  if (permutations < 0 || !std::isfinite(seed) || seed < 0.0) {
    Rcpp::stop("permutations and seed must be nonnegative");
  }
  arma::mat weights(k, n, arma::fill::zeros);
  arma::vec raw_weights(k, arma::fill::zeros);
  double total_weight = 0.0;
  for (arma::uword target = 0; target < n; ++target) {
    double bandwidth_used = 0.0;
    gwrs::compute_kernel_weights(
      neighbor_distance, target, kernel, adaptive, fixed_bandwidth,
      raw_weights, bandwidth_used
    );
    double row_sum = 0.0;
    for (arma::uword row = 0; row < k; ++row) {
      const arma::uword encoded = neighbor_index(row, target);
      if (encoded == 0 || encoded > n) {
        Rcpp::stop("invalid neighbor index");
      }
      const double value = encoded == target + 1 ? 0.0 : raw_weights[row];
      weights(row, target) = value;
      row_sum += value;
    }
    if (row_standardize && row_sum > gwrs::kNumericalEpsilon) {
      weights.col(target) /= row_sum;
      row_sum = 1.0;
    }
    total_weight += row_sum;
  }
  if (total_weight <= gwrs::kNumericalEpsilon) {
    Rcpp::stop("spatial weights contain no non-self connections");
  }

  const arma::vec z = residuals - arma::mean(residuals);
  const double denominator = arma::dot(z, z);
  if (denominator <= gwrs::kNumericalEpsilon) {
    Rcpp::stop("Moran's I is undefined for constant residuals");
  }
  const double m2 = denominator / static_cast<double>(n);
  arma::vec lag(n, arma::fill::zeros);
  for (arma::uword target = 0; target < n; ++target) {
    for (arma::uword row = 0; row < k; ++row) {
      lag[target] += weights(row, target) *
        z[neighbor_index(row, target) - 1];
    }
  }
  const double observed_global = static_cast<double>(n) / total_weight *
    arma::dot(z, lag) / denominator;
  const arma::vec observed_local = z % lag / m2;
  arma::vec global_permuted;
  arma::vec local_p(n, arma::fill::value(arma::datum::nan));
  arma::vec local_z(n, arma::fill::value(arma::datum::nan));
  double global_p = arma::datum::nan;

  if (permutations > 0) {
    global_permuted.set_size(permutations);
    const std::size_t moment_blocks = std::min(
      static_cast<std::size_t>(permutations),
      kDeterministicMoranMomentBlocks
    );
    arma::mat block_sum(n, moment_blocks, arma::fill::zeros);
    arma::mat block_sumsq(n, moment_blocks, arma::fill::zeros);
    MoranMomentBlockWorker moments(
      z, neighbor_index, weights, m2,
      static_cast<std::uint64_t>(seed),
      static_cast<std::size_t>(permutations), moment_blocks,
      global_permuted, block_sum, block_sumsq
    );
    if (n_threads == 1 || moment_blocks < 2) {
      moments(0, moment_blocks);
    } else {
      RcppParallel::parallelFor(
        0, moment_blocks, moments, 1, n_threads
      );
    }
    global_permuted *= static_cast<double>(n) / total_weight;
    arma::vec local_sum(n, arma::fill::zeros);
    arma::vec local_sumsq(n, arma::fill::zeros);
    for (std::size_t block = 0; block < moment_blocks; ++block) {
      local_sum += block_sum.col(block);
      local_sumsq += block_sumsq.col(block);
    }
    const arma::vec local_mean = local_sum /
      static_cast<double>(permutations);
    arma::vec local_variance(n, arma::fill::zeros);
    if (permutations > 1) {
      local_variance = (local_sumsq -
        static_cast<double>(permutations) * arma::square(local_mean)) /
        static_cast<double>(permutations - 1);
      local_variance = arma::clamp(local_variance, 0.0, arma::datum::inf);
      const arma::vec local_sd = arma::sqrt(local_variance);
      const arma::uvec valid_sd = arma::find(local_sd >
        gwrs::kNumericalEpsilon);
      local_z.elem(valid_sd) =
        (observed_local.elem(valid_sd) - local_mean.elem(valid_sd)) /
        local_sd.elem(valid_sd);
    }
    MoranCountWorker counts(
      z, neighbor_index, weights, m2,
      static_cast<std::uint64_t>(seed), observed_local, local_mean
    );
    if (n_threads == 1 || permutations < 2) {
      counts(0, permutations);
    } else {
      RcppParallel::parallelReduce(
        0, permutations, counts,
        static_cast<std::size_t>(std::max(1, grain_size)), n_threads
      );
    }
    local_p = (arma::conv_to<arma::vec>::from(counts.exceedance) + 1.0) /
      static_cast<double>(permutations + 1);
    const double expected = -1.0 / static_cast<double>(n - 1);
    const arma::uword global_exceed = arma::accu(
      arma::abs(global_permuted - expected) >=
        std::abs(observed_global - expected)
    );
    global_p = static_cast<double>(global_exceed + 1) /
      static_cast<double>(permutations + 1);
  }

  return Rcpp::List::create(
    Rcpp::_ ["I"] = observed_global,
    Rcpp::_ ["expected_I"] = -1.0 / static_cast<double>(n - 1),
    Rcpp::_ ["p_value"] = global_p,
    Rcpp::_ ["local_I"] = observed_local,
    Rcpp::_ ["local_z"] = local_z,
    Rcpp::_ ["local_p"] = local_p,
    Rcpp::_ ["centered_residual"] = z,
    Rcpp::_ ["spatial_lag"] = lag,
    Rcpp::_ ["permuted_I"] = global_permuted,
    Rcpp::_ ["total_weight"] = total_weight
  );
}
