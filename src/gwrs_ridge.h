#ifndef GWRS_RIDGE_H
#define GWRS_RIDGE_H

#include "gwrs_common.h"

namespace gwrs {

// X and y are already weighted by sqrt(w / sum(w)) after weighted centering.
// Cache the smaller Gram system once per target, then reuse it across lambda.
// Only pure, positive Ridge uses this solver. No unpenalized approximation,
// diagonal jitter, pseudoinverse, or fallback estimator is introduced.
class RidgeSystem {
  bool dual_ = false;
  bool ready_ = false;
  arma::mat gram_;
  arma::vec rhs_;

public:
  bool prepare(const arma::mat& weighted_x, const arma::vec& weighted_y) {
    ready_ = false;
    if (weighted_x.n_rows != weighted_y.n_elem ||
        !weighted_x.is_finite() || !weighted_y.is_finite()) return false;
    dual_ = weighted_x.n_cols > weighted_x.n_rows;
    if (dual_) {
      gram_ = weighted_x * weighted_x.t();
      rhs_ = weighted_y;
    } else {
      gram_ = weighted_x.t() * weighted_x;
      rhs_ = weighted_x.t() * weighted_y;
    }
    ready_ = gram_.is_finite() && rhs_.is_finite();
    return ready_;
  }

  bool solve(const arma::mat& weighted_x, const double lambda,
             arma::vec& beta) const {
    if (!ready_ || !std::isfinite(lambda) || lambda <= 0.0) return false;
    arma::mat system = gram_;
    system.diag() += lambda;
    if (!system.is_finite()) return false;
    arma::mat lower;
    if (!arma::chol(lower, system, "lower")) return false;
    arma::vec intermediate;
    arma::vec solution;
    const auto options = arma::solve_opts::fast + arma::solve_opts::no_approx;
    if (!arma::solve(intermediate, arma::trimatl(lower), rhs_, options) ||
        !arma::solve(solution, arma::trimatu(lower.t()), intermediate, options) ||
        !solution.is_finite()) return false;
    if (dual_) beta = weighted_x.t() * solution;
    else beta = solution;
    return beta.is_finite();
  }
};

} // namespace gwrs
#endif
