#ifndef GWRS_NONCONVEX_BLOCK_H
#define GWRS_NONCONVEX_BLOCK_H

#include "gwrs_nonconvex.h"
#include <algorithm>
#include <limits>
#include <vector>

namespace gwrs {

// Frozen numerical scheduling constants, mirrored in v10/config/solver.json.
constexpr unsigned int kPairShortlist = 8;
constexpr int kPairEverySweeps = 5;
constexpr double kPairMinimumCorrelation = 0.99;
constexpr double kPairGuardEpsilon = 64.0;

struct PairSegment {
  double lower, upper, curvature, linear;
};

struct PairCandidate {
  double first, second;
  bool finite;
  unsigned int evaluated;
};

inline std::vector<PairSegment> pair_segments(const double lambda,
                                             const double gamma,
                                             const int penalty) {
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<PairSegment> pieces;
  for (int sign : {-1, 1}) {
    auto add = [&](const double lo, const double hi, const double curvature,
                   const double linear) {
      pieces.push_back({sign < 0 ? -hi : lo, sign < 0 ? -lo : hi,
                        curvature, sign * linear});
    };
    if (penalty == kMcpPenalty) {
      add(0.0, gamma * lambda, -1.0 / gamma, lambda);
    } else {
      add(0.0, lambda, 0.0, lambda);
      add(lambda, gamma * lambda, -1.0 / (gamma - 1.0),
          gamma * lambda / (gamma - 1.0));
    }
    add(gamma * lambda, inf, 0.0, 0.0);
  }
  return pieces;
}

// Enumerate interior minima on signed quadratic pieces and minima on every
// finite boundary. On a boundary the original exact scalar solver applies.
// An indefinite interior cannot be a minimum; a singular flat interior has
// an equally good boundary point. No ridge jitter or penalty rescaling.
// Near-singular positive systems use long-double determinant arithmetic;
// unresolved cases retain the current point for the external objective guard.
inline PairCandidate nonconvex_pair_candidate(
    const double q1, const double q12, const double q2,
    const double z1, const double z2, const double lambda,
    const double gamma, const int penalty, const double current1,
    const double current2) {
  PairCandidate result{current1, current2, false, 0};
  if (!valid_nonconvex_penalty(penalty) || !std::isfinite(lambda) ||
      lambda <= 0.0 || !std::isfinite(gamma) ||
      gamma <= (penalty == kMcpPenalty ? 1.0 : 2.0) ||
      !std::isfinite(q1) || !std::isfinite(q2) || !std::isfinite(q12) ||
      !std::isfinite(z1) || !std::isfinite(z2) ||
      !std::isfinite(current1) || !std::isfinite(current2) ||
      q1 <= kNumericalEpsilon || q2 <= kNumericalEpsilon) return result;
  using LD = long double;
  const LD determinant = LD(q1) * q2 - LD(q12) * q12;
  const LD scale = LD(q1) * q2 + LD(q12) * q12;
  if (determinant < -64 * std::numeric_limits<double>::epsilon() * scale)
    return result;
  auto value = [&](double u, double v) -> LD {
    return LD(0.5) * q1 * u * u + LD(q12) * u * v + LD(0.5) * q2 * v * v -
      LD(z1) * u - LD(z2) * v +
      nonconvex_penalty_value(u, lambda, gamma, penalty) +
      nonconvex_penalty_value(v, lambda, gamma, penalty);
  };
  LD best = value(current1, current2);
  if (!std::isfinite(best)) return result;
  result.finite = true;
  auto consider = [&](const double u, const double v) {
    if (!std::isfinite(u) || !std::isfinite(v)) return;
    const LD objective = value(u, v);
    ++result.evaluated;
    // Keep current/earlier candidates on exact ties; no gratuitous swapping.
    if (std::isfinite(objective) && objective < best) {
      best = objective; result.first = u; result.second = v;
    }
  };
  std::vector<double> knots{0.0, -gamma * lambda, gamma * lambda};
  if (penalty == kScadPenalty) {
    knots.push_back(-lambda); knots.push_back(lambda);
  }
  consider(0.0, 0.0);
  for (const double knot : knots) {
    consider(knot, nonconvex_coordinate_minimum(z2 - q12 * knot, q2,
                                               lambda, gamma, penalty));
    consider(nonconvex_coordinate_minimum(z1 - q12 * knot, q1,
                                          lambda, gamma, penalty), knot);
  }
  const auto pieces = pair_segments(lambda, gamma, penalty);
  for (const auto& a : pieces) for (const auto& b : pieces) {
    const LD h1 = LD(q1) + a.curvature, h2 = LD(q2) + b.curvature;
    const LD det = h1 * h2 - LD(q12) * q12;
    const LD det_scale = std::abs(h1 * h2) + LD(q12) * q12;
    if (h1 <= 0 || h2 <= 0 || det <= 64 * std::numeric_limits<LD>::epsilon() * det_scale)
      continue;
    const LD s1 = LD(z1) - a.linear, s2 = LD(z2) - b.linear;
    const double u = static_cast<double>((h2 * s1 - LD(q12) * s2) / det);
    const double v = static_cast<double>((h1 * s2 - LD(q12) * s1) / det);
    if (u >= a.lower && u <= a.upper && v >= b.lower && v <= b.upper)
      consider(u, v);
  }
  return result;
}

inline double pair_full_objective(const arma::vec& residual,
                                  const arma::vec& weights,
                                  const arma::vec& beta, const double lambda,
                                  const double gamma, const int penalty) {
  double objective = 0.5 * arma::dot(weights % residual, residual);
  for (arma::uword j = 0; j < beta.n_elem; ++j)
    objective += nonconvex_penalty_value(beta[j], lambda, gamma, penalty);
  return objective;
}

// Transactional update: a rejected/nonfinite proposal changes neither state.
inline bool apply_nonconvex_pair(const arma::mat& x, const arma::vec& weights,
                                 arma::vec& beta, arma::vec& residual,
                                 const arma::uword j, const arma::uword k,
                                 const PairCandidate& candidate,
                                 const double lambda, const double gamma,
                                 const int penalty, double& gain,
                                 double& maximum_delta) {
  gain = 0.0; maximum_delta = 0.0;
  if (!candidate.finite || !std::isfinite(candidate.first) ||
      !std::isfinite(candidate.second) || j == k || j >= beta.n_elem ||
      k >= beta.n_elem) return false;
  const double d1 = candidate.first - beta[j], d2 = candidate.second - beta[k];
  arma::vec proposed_beta = beta;
  proposed_beta[j] = candidate.first; proposed_beta[k] = candidate.second;
  const arma::vec proposed_residual = residual - x.col(j) * d1 - x.col(k) * d2;
  if (!proposed_residual.is_finite() || !proposed_beta.is_finite()) return false;
  const double before = pair_full_objective(residual, weights, beta, lambda, gamma, penalty);
  const double after = pair_full_objective(proposed_residual, weights, proposed_beta, lambda, gamma, penalty);
  const double margin = kPairGuardEpsilon * std::numeric_limits<double>::epsilon() *
    (1.0 + std::abs(before));
  if (!std::isfinite(before) || !std::isfinite(after) || !(after < before - margin)) return false;
  beta = proposed_beta; residual = proposed_residual;
  gain = before - after; maximum_delta = std::max(std::abs(d1), std::abs(d2));
  return true;
}

// Bound the work to O(p log L + k L^2), L <= 8, without allocating p-by-p.
// This changes scheduling only; excluded coordinates still get scalar updates
// and the existing full inactive-coordinate improvement scan.
inline bool select_nonconvex_pair(const arma::mat& x, const arma::vec& weights,
                                  const arma::vec& curvature,
                                  const arma::vec& beta, const arma::vec& deltas,
                                  arma::uword& first, arma::uword& second,
                                  double& cross_curvature) {
  std::vector<arma::uword> shortlist;
  auto better = [&](arma::uword a, arma::uword b) {
    return deltas[a] > deltas[b] || (deltas[a] == deltas[b] && a < b);
  };
  for (arma::uword j = 0; j < beta.n_elem; ++j) {
    if (std::abs(beta[j]) <= 1e-12 || curvature[j] <= kNumericalEpsilon) continue;
    shortlist.insert(std::lower_bound(shortlist.begin(), shortlist.end(), j, better), j);
    if (shortlist.size() > kPairShortlist) shortlist.pop_back();
  }
  double best = kPairMinimumCorrelation;
  bool found = false;
  for (unsigned int a = 0; a < shortlist.size(); ++a)
    for (unsigned int b = a + 1; b < shortlist.size(); ++b) {
      const arma::uword j = std::min(shortlist[a], shortlist[b]);
      const arma::uword k = std::max(shortlist[a], shortlist[b]);
      const double cross = arma::dot(x.col(j), weights % x.col(k));
      const double rho = std::abs(cross / std::sqrt(curvature[j] * curvature[k]));
      if (std::isfinite(rho) && (rho > best || (!found && rho == best) ||
          (found && rho == best && (j < first || (j == first && k < second))))) {
        found = true; best = rho; first = j; second = k; cross_curvature = cross;
      }
    }
  return found;
}

} // namespace gwrs
#endif
