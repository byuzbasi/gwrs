#ifndef GWRS_NONCONVEX_H
#define GWRS_NONCONVEX_H

#include "gwrs_common.h"

#include <array>

namespace gwrs {

constexpr int kMcpPenalty = 0;
constexpr int kScadPenalty = 1;

inline bool valid_nonconvex_penalty(const int penalty) {
  return penalty == kMcpPenalty || penalty == kScadPenalty;
}

inline double nonconvex_penalty_value(const double coefficient,
                                      const double lambda,
                                      const double gamma,
                                      const int penalty) {
  const double value = std::abs(coefficient);
  if (lambda <= 0.0 || value <= 0.0) return 0.0;

  if (penalty == kMcpPenalty) {
    if (value <= gamma * lambda) {
      return lambda * value - value * value / (2.0 * gamma);
    }
    return 0.5 * gamma * lambda * lambda;
  }

  if (value <= lambda) return lambda * value;
  if (value <= gamma * lambda) {
    return (-value * value + 2.0 * gamma * lambda * value -
            lambda * lambda) / (2.0 * (gamma - 1.0));
  }
  return 0.5 * (gamma + 1.0) * lambda * lambda;
}

inline double nonconvex_penalty_derivative(const double absolute_coefficient,
                                           const double lambda,
                                           const double gamma,
                                           const int penalty) {
  if (lambda <= 0.0) return 0.0;
  if (penalty == kMcpPenalty) {
    if (absolute_coefficient >= gamma * lambda) return 0.0;
    return std::max(0.0, lambda - absolute_coefficient / gamma);
  }
  if (absolute_coefficient <= lambda) return lambda;
  if (absolute_coefficient < gamma * lambda) {
    return (gamma * lambda - absolute_coefficient) / (gamma - 1.0);
  }
  return 0.0;
}

inline double nonconvex_penalty_second_derivative(
    const double absolute_coefficient,
    const double lambda,
    const double gamma,
    const int penalty) {
  if (lambda <= 0.0) return 0.0;
  if (penalty == kMcpPenalty) {
    return (absolute_coefficient > 0.0 &&
            absolute_coefficient < gamma * lambda) ? -1.0 / gamma : 0.0;
  }
  return (absolute_coefficient > lambda &&
          absolute_coefficient < gamma * lambda) ?
    -1.0 / (gamma - 1.0) : 0.0;
}

inline bool nonconvex_penalty_knot(const double absolute_coefficient,
                                   const double lambda,
                                   const double gamma,
                                   const int penalty) {
  if (lambda <= 0.0) return false;
  const double tolerance = 1e-8 * (1.0 + gamma * lambda);
  if (penalty == kMcpPenalty) {
    return std::abs(absolute_coefficient - gamma * lambda) <= tolerance;
  }
  return std::abs(absolute_coefficient - lambda) <= tolerance ||
    std::abs(absolute_coefficient - gamma * lambda) <= tolerance;
}

inline double coordinate_objective(const double coefficient,
                                   const double score,
                                   const double curvature,
                                   const double lambda,
                                   const double gamma,
                                   const int penalty) {
  return 0.5 * curvature * coefficient * coefficient -
    score * coefficient +
    nonconvex_penalty_value(coefficient, lambda, gamma, penalty);
}

inline void consider_coordinate_value(const double candidate,
                                      const double objective,
                                      double& best,
                                      double& best_objective) {
  if (!std::isfinite(candidate)) return;
  const double scale = 1.0 + std::abs(best_objective);
  if (objective < best_objective - 1e-14 * scale ||
      (std::abs(objective - best_objective) <= 1e-14 * scale &&
       std::abs(candidate) < std::abs(best))) {
    best = candidate;
    best_objective = objective;
  }
}

inline void consider_coordinate_candidate(const double candidate,
                                           const double score,
                                           const double curvature,
                                           const double lambda,
                                           const double gamma,
                                           const int penalty,
                                           double& best,
                                           double& best_objective) {
  if (!std::isfinite(candidate)) return;
  consider_coordinate_value(candidate, coordinate_objective(
    candidate, score, curvature, lambda, gamma, penalty), best, best_objective);
}

// One coordinate at one target/lambda: curvature and penalty settings stay fixed.
// Scores, candidates, comparison order and tie rules are never cached.
struct NonconvexCoordinateCache {
  double upper, denominator, middle_offset;
  double first_quadratic, first_penalty, upper_quadratic, upper_penalty;
  NonconvexCoordinateCache(double curvature, double lambda, double gamma, int penalty)
    : upper(gamma * lambda),
      denominator(curvature - (penalty == kMcpPenalty ? 1.0 / gamma : 1.0 / (gamma - 1.0))),
      middle_offset(penalty == kScadPenalty ? gamma * lambda / (gamma - 1.0) : 0.0),
      first_quadratic(0.5 * curvature * lambda * lambda),
      first_penalty(nonconvex_penalty_value(lambda, lambda, gamma, penalty)),
      upper_quadratic(0.5 * curvature * upper * upper),
      upper_penalty(nonconvex_penalty_value(upper, lambda, gamma, penalty)) {}
};

// Exact global minimizer of a one-coordinate SCAD/MCP subproblem:
//   0.5 * curvature * b^2 - score * b + p_lambda,gamma(|b|).
// Piecewise stationary points and every penalty knot are compared directly.
// This remains valid when local GWR curvature is not one and when a concave
// penalty segment makes a familiar closed-form denominator nonpositive.
inline double nonconvex_coordinate_minimum(const double score,
                                           const double curvature,
                                           const double lambda,
                                           const double gamma,
                                           const int penalty,
                                           const NonconvexCoordinateCache* cache = nullptr) {
  if (!std::isfinite(score) || !std::isfinite(curvature) ||
      curvature <= kNumericalEpsilon) return 0.0;
  if (lambda <= 0.0) return score / curvature;

  const double sign = score < 0.0 ? -1.0 : 1.0;
  const double absolute_score = std::abs(score);
  double best = 0.0;
  double best_objective = 0.0;
  auto consider_absolute = [&](const double value) {
    if (!std::isfinite(value) || value < 0.0) return;
    consider_coordinate_candidate(
      sign * value, score, curvature, lambda, gamma, penalty,
      best, best_objective
    );
  };
  auto consider_knot = [&](const double value, const double quadratic,
                           const double penalty_value) {
    if (!std::isfinite(value) || value < 0.0) return;
    const double candidate = sign * value;
    // Preserve the original quadratic - linear + penalty evaluation order.
    consider_coordinate_value(candidate, quadratic - score * candidate + penalty_value,
                              best, best_objective);
  };

  if (penalty == kMcpPenalty) {
    const double upper = cache ? cache->upper : gamma * lambda;
    if (cache) consider_knot(upper, cache->upper_quadratic, cache->upper_penalty);
    else consider_absolute(upper);
    const double denominator = cache ? cache->denominator : curvature - 1.0 / gamma;
    if (std::abs(denominator) > kNumericalEpsilon) {
      const double stationary = (absolute_score - lambda) / denominator;
      if (stationary >= 0.0 && stationary <= upper) {
        consider_absolute(stationary);
      }
    }
    const double unpenalized = absolute_score / curvature;
    if (unpenalized >= upper) consider_absolute(unpenalized);
    return best;
  }

  const double first_knot = lambda;
  const double second_knot = cache ? cache->upper : gamma * lambda;
  if (cache) {
    consider_knot(first_knot, cache->first_quadratic, cache->first_penalty);
    consider_knot(second_knot, cache->upper_quadratic, cache->upper_penalty);
  } else {
    consider_absolute(first_knot);
    consider_absolute(second_knot);
  }

  const double first_stationary = (absolute_score - lambda) / curvature;
  if (first_stationary >= 0.0 && first_stationary <= first_knot) {
    consider_absolute(first_stationary);
  }

  const double denominator = cache ? cache->denominator : curvature - 1.0 / (gamma - 1.0);
  if (std::abs(denominator) > kNumericalEpsilon) {
    const double middle_stationary =
      (absolute_score - (cache ? cache->middle_offset : gamma * lambda / (gamma - 1.0))) / denominator;
    if (middle_stationary >= first_knot &&
        middle_stationary <= second_knot) {
      consider_absolute(middle_stationary);
    }
  }

  const double unpenalized = absolute_score / curvature;
  if (unpenalized >= second_knot) consider_absolute(unpenalized);
  return best;
}

inline double nonconvex_coordinate_improvement(const double coefficient,
                                               const double score,
                                               const double curvature,
                                               const double lambda,
                                               const double gamma,
                                               const int penalty,
                                               const NonconvexCoordinateCache* cache = nullptr) {
  const double updated = nonconvex_coordinate_minimum(
    score, curvature, lambda, gamma, penalty, cache
  );
  const double current_objective = coordinate_objective(
    coefficient, score, curvature, lambda, gamma, penalty
  );
  const double updated_objective = coordinate_objective(
    updated, score, curvature, lambda, gamma, penalty
  );
  return std::max(0.0, current_objective - updated_objective);
}

inline double nonconvex_stationarity_violation(const double coefficient,
                                               const double gradient_score,
                                               const double lambda,
                                               const double gamma,
                                               const int penalty) {
  if (std::abs(coefficient) <= 1e-12) {
    return std::max(0.0, std::abs(gradient_score) - lambda);
  }
  const double derivative = nonconvex_penalty_derivative(
    std::abs(coefficient), lambda, gamma, penalty
  );
  const double signed_derivative = coefficient > 0.0 ? derivative : -derivative;
  return std::abs(gradient_score - signed_derivative);
}

} // namespace gwrs

#endif
