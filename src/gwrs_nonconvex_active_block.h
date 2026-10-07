#ifndef GWRS_NONCONVEX_ACTIVE_BLOCK_H
#define GWRS_NONCONVEX_ACTIVE_BLOCK_H

#include "gwrs_nonconvex_block.h"

namespace gwrs {
constexpr arma::uword kActiveBlockMin = 3;
constexpr arma::uword kActiveBlockMax = 32;
// Uses the existing pair schedule, never adds a sweep or relaxes a stop.
static_assert(kPairEverySweeps == 5, "Reviewed block schedule changed");

// 0 not scheduled/disabled; 1 size; 2 current region; 3 rank/QR;
// 4 solve; 5 reserved (v13-only candidate veto); 6 nondecrease;
// 7 accepted; 8 invalid input. Status numbers are not renumbered.
struct ActiveBlockResult {
  int status = 0;
  arma::uword size = 0;
  double gain = 0.0, delta = 0.0, rcond = arma::datum::nan;
};

// Owned by one target of one worker: x and weights are constant for its
// lifetime. Only the last successful factorization is retained (<=32 columns).
// Coefficients, residuals, lambda/gamma, proposals and acceptance are not cached.
struct ActiveBlockQRCache {
  arma::uvec index;
  arma::mat q, r;
  double rcond = arma::datum::nan;
  bool valid = false;
  std::size_t hits = 0, factorizations = 0;
};

inline ActiveBlockResult apply_nonconvex_active_block(
    const arma::mat& x, const arma::vec& weights, arma::vec& beta,
    arma::vec& residual, double lambda, double gamma, int penalty,
    ActiveBlockQRCache* cache = nullptr) {
  ActiveBlockResult result;
  if (x.n_rows != weights.n_elem || x.n_rows != residual.n_elem ||
      x.n_cols != beta.n_elem || !weights.is_finite() || arma::any(weights < 0) ||
      !beta.is_finite() || !residual.is_finite() || !std::isfinite(lambda) ||
      lambda <= 0 || !std::isfinite(gamma) || !valid_nonconvex_penalty(penalty)) {
    result.status = 8; return result;
  }
  std::vector<arma::uword> selected;
  selected.reserve(kActiveBlockMax);
  for (arma::uword j = 0; j < beta.n_elem; ++j) {
    if (std::abs(beta[j]) <= 1e-12) continue; // Existing solver activity threshold.
    ++result.size;
    if (selected.size() < kActiveBlockMax) selected.push_back(j);
  }
  if (result.size < kActiveBlockMin || result.size > kActiveBlockMax ||
      x.n_rows < result.size) { result.status = 1; return result; }
  const double eps = std::numeric_limits<double>::epsilon();
  const double edge = gamma * lambda;
  const double margin = 64.0 * eps * (1.0 + std::abs(edge));
  const arma::uvec index(selected);
  const arma::vec current = beta.elem(index);
  if (!std::isfinite(edge) || arma::any(arma::abs(current) <= edge + margin)) {
    result.status = 2; return result;
  }

  const arma::mat block = x.cols(index);
  if (!block.is_finite()) { result.status = 8; return result; }
  const arma::vec sqrt_weights = arma::sqrt(weights);
  // Other coefficients, including tiny nonzero ones, are held fixed; not deleted.
  const arma::vec rhs = sqrt_weights % (residual + block * current);
  arma::mat local_q, local_r;
  arma::mat& q = cache ? cache->q : local_q;
  arma::mat& r = cache ? cache->r : local_r;
  const bool hit = cache && cache->valid && cache->index.n_elem == index.n_elem &&
    arma::all(cache->index == index);
  if (hit) {
    result.rcond = cache->rcond;
    ++cache->hits;
  } else {
    if (cache) { cache->valid = false; ++cache->factorizations; }
    arma::mat z = block;
    z.each_col() %= sqrt_weights;
    if (!arma::qr_econ(q, r, z) || !q.is_finite() || !r.is_finite()) {
      result.status = 3; return result;
    }
    result.rcond = arma::rcond(r);
    const double rank_floor = 64.0 * eps * std::max(z.n_rows, z.n_cols);
    if (!std::isfinite(result.rcond) || result.rcond <= rank_floor) {
      result.status = 3; return result;
    }
    if (cache) {
      cache->index = index; cache->rcond = result.rcond; cache->valid = true;
    }
  }
  arma::vec candidate;
  if (!arma::solve(candidate, arma::trimatu(r), q.t() * rhs,
                   arma::solve_opts::no_approx) || !candidate.is_finite()) {
    result.status = 4; return result;
  }
  // The current block lies in the flat tails, where each SCAD/MCP penalty
  // attains its global upper bound. Conditional weighted least squares cannot
  // increase the loss in exact arithmetic, and candidate penalties cannot
  // exceed those current bounds, even across knots or signs. Thus retaining
  // the current sign/region is not necessary for descent. The finite-precision
  // full original-objective guard below remains mandatory. This is only an
  // internal proposal; scalar/screening/convergence checks still follow.
  arma::vec next_beta = beta;
  next_beta.elem(index) = candidate;
  const arma::vec next_residual = residual - block * (candidate - current);
  const double before = pair_full_objective(residual, weights, beta, lambda, gamma, penalty);
  const double after = pair_full_objective(next_residual, weights, next_beta, lambda, gamma, penalty);
  const double decrease_floor = 64.0 * eps * (1.0 + std::abs(before));
  if (!std::isfinite(before) || !std::isfinite(after) ||
      !next_residual.is_finite() || !(before - after > decrease_floor)) {
    result.status = 6; return result;
  }
  result.gain = before - after;
  result.delta = arma::abs(candidate - current).max();
  result.status = 7;
  // Publish only after every guard has passed. Rejection is non-mutating.
  beta = std::move(next_beta);
  residual = next_residual;
  return result;
}
} // namespace gwrs
#endif
