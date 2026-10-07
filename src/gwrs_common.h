#ifndef GWRS_COMMON_H
#define GWRS_COMMON_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace gwrs {

constexpr double kNumericalEpsilon = 1e-12;

inline double soft_threshold(const double value, const double threshold) {
  if (value > threshold) return value - threshold;
  if (value < -threshold) return value + threshold;
  return 0.0;
}

inline double kernel_weight(const double distance,
                            const double bandwidth,
                            const int kernel) {
  if (!std::isfinite(distance) || distance < 0.0 || bandwidth <= 0.0) {
    return 0.0;
  }

  const double ratio = distance / bandwidth;
  switch (kernel) {
  case 0: // Gaussian
    return std::exp(-0.5 * ratio * ratio);
  case 1: // Bisquare
    if (ratio >= 1.0) return 0.0;
    return std::pow(1.0 - ratio * ratio, 2.0);
  case 2: // Exponential
    return std::exp(-ratio);
  case 3: // Tricube
    if (ratio >= 1.0) return 0.0;
    return std::pow(1.0 - std::pow(ratio, 3.0), 3.0);
  case 4: // Boxcar
    return ratio <= 1.0 ? 1.0 : 0.0;
  default:
    return 0.0;
  }
}

inline double adaptive_bandwidth(const arma::mat& distances,
                                 const arma::uword target) {
  double bandwidth = 0.0;
  for (arma::uword j = 0; j < distances.n_rows; ++j) {
    const double value = distances(j, target);
    if (std::isfinite(value) && value > bandwidth) bandwidth = value;
  }
  if (bandwidth <= kNumericalEpsilon) bandwidth = 1.0;
  // Keep the most distant retained neighbor inside compact-support kernels.
  return bandwidth * (1.0 + 1e-10);
}

inline double compute_kernel_weights(const arma::mat& distances,
                                     const arma::uword target,
                                     const int kernel,
                                     const bool adaptive,
                                     const double fixed_bandwidth,
                                     arma::vec& weights,
                                     double& bandwidth_used) {
  bandwidth_used = adaptive ? adaptive_bandwidth(distances, target)
                            : fixed_bandwidth;
  double weight_sum = 0.0;
  for (arma::uword j = 0; j < distances.n_rows; ++j) {
    const double weight = kernel_weight(distances(j, target), bandwidth_used,
                                        kernel);
    weights[j] = weight;
    weight_sum += weight;
  }
  return weight_sum;
}

inline bool finite_matrix(const arma::mat& x) {
  return x.is_finite();
}

inline bool finite_vector(const arma::vec& x) {
  return x.is_finite();
}

inline bool solve_local_system(const arma::mat& gram,
                               const arma::vec& score,
                               arma::vec& beta,
                               int& status) {
  status = 0;
  bool ok = arma::solve(beta, gram, score,
                        arma::solve_opts::likely_sympd +
                          arma::solve_opts::no_approx);
  if (ok && beta.is_finite()) return true;

  status = 1;
  const double scale = std::max(1.0, arma::norm(gram, "inf"));
  beta = arma::pinv(gram, scale * 1e-10) * score;
  if (beta.is_finite()) return true;

  status = 2;
  beta.zeros(score.n_elem);
  return false;
}

inline double kkt_residual(const arma::mat& x_centered,
                           const arma::vec& residual,
                           const arma::vec& normalized_weights,
                           const arma::vec& beta,
                           const arma::vec& center,
                           const double lambda1,
                           const double lambda2) {
  double maximum = 0.0;
  for (arma::uword j = 0; j < beta.n_elem; ++j) {
    const double score = arma::dot(x_centered.col(j),
                                   normalized_weights % residual) +
      lambda2 * (center[j] - beta[j]);
    double violation = 0.0;
    if (std::abs(beta[j]) <= 1e-12) {
      violation = std::max(0.0, std::abs(score) - lambda1);
    } else {
      violation = std::abs(score - lambda1 * (beta[j] > 0.0 ? 1.0 : -1.0));
    }
    if (violation > maximum) maximum = violation;
  }
  return maximum;
}

} // namespace gwrs

#endif
