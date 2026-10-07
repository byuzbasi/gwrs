#include "gwrs_common.h"

#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

namespace {

struct KtddPrepareWorker : public RcppParallel::Worker {
  const arma::mat& design;
  const arma::vec& y;
  const arma::umat& neighbor_index;
  const arma::mat& neighbor_distance;
  const int kernel;
  const bool adaptive;
  const double fixed_bandwidth;
  const double gamma;
  const arma::sp_mat& coefficient_laplacian;
  const bool cache_factor;

  arma::mat& normalized_weights;
  arma::mat& local_score;
  arma::cube& factor;
  arma::ivec& factor_type;
  arma::mat& coefficient_preconditioner;
  arma::mat& initial_coefficients;
  arma::vec& weight_sum;
  arma::vec& bandwidth_used;
  arma::ivec& status;

  KtddPrepareWorker(const arma::mat& design,
                    const arma::vec& y,
                    const arma::umat& neighbor_index,
                    const arma::mat& neighbor_distance,
                    const int kernel,
                    const bool adaptive,
                    const double fixed_bandwidth,
                    const double gamma,
                    const arma::sp_mat& coefficient_laplacian,
                    const bool cache_factor,
                    arma::mat& normalized_weights,
                    arma::mat& local_score,
                    arma::cube& factor,
                    arma::ivec& factor_type,
                    arma::mat& coefficient_preconditioner,
                    arma::mat& initial_coefficients,
                    arma::vec& weight_sum,
                    arma::vec& bandwidth_used,
                    arma::ivec& status)
    : design(design), y(y), neighbor_index(neighbor_index),
      neighbor_distance(neighbor_distance), kernel(kernel),
      adaptive(adaptive), fixed_bandwidth(fixed_bandwidth), gamma(gamma),
      coefficient_laplacian(coefficient_laplacian),
      cache_factor(cache_factor),
      normalized_weights(normalized_weights), local_score(local_score),
      factor(factor), factor_type(factor_type),
      coefficient_preconditioner(coefficient_preconditioner),
      initial_coefficients(initial_coefficients), weight_sum(weight_sum),
      bandwidth_used(bandwidth_used), status(status) {}

  void operator()(std::size_t begin, std::size_t end) {
    const arma::uword k = neighbor_index.n_rows;
    const arma::uword p1 = design.n_cols;
    arma::mat local_design(k, p1, arma::fill::zeros);
    arma::mat weighted_design(k, p1, arma::fill::zeros);
    arma::mat gram(p1, p1, arma::fill::zeros);
    arma::mat system(p1, p1, arma::fill::zeros);
    arma::mat local_factor(p1, p1, arma::fill::zeros);
    arma::vec local_y(k, arma::fill::zeros);
    arma::vec weights(k, arma::fill::zeros);
    arma::vec sqrt_weights(k, arma::fill::zeros);
    arma::vec score(p1, arma::fill::zeros);
    arma::vec beta(p1, arma::fill::zeros);

    for (std::size_t target = begin; target < end; ++target) {
      bool bad_index = false;
      for (arma::uword row = 0; row < k; ++row) {
        const arma::uword encoded = neighbor_index(row, target);
        if (encoded == 0 || encoded > design.n_rows) {
          bad_index = true;
          break;
        }
        const arma::uword source = encoded - 1;
        local_y[row] = y[source];
        local_design.row(row) = design.row(source);
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

      weights /= local_weight_sum;
      normalized_weights.col(target) = weights;
      sqrt_weights = arma::sqrt(weights);
      weighted_design = local_design;
      weighted_design.each_col() %= sqrt_weights;
      gram = weighted_design.t() * weighted_design;
      score = local_design.t() * (weights % local_y);
      local_score.col(target) = score;

      int local_status = 0;
      gwrs::solve_local_system(gram, score, beta, local_status);
      initial_coefficients.col(target) = beta;
      status[target] = local_status;

      system = gram;
      const double graph_degree = coefficient_laplacian(target, target);
      system.diag() += gamma * graph_degree;
      if (cache_factor) {
        if (arma::chol(local_factor, system, "lower")) {
          factor.slice(target) = local_factor;
          factor_type[target] = 0;
        } else {
          factor.slice(target) = arma::pinv(system);
          factor_type[target] = 1;
        }
      } else {
        coefficient_preconditioner.col(target) = system.diag();
        coefficient_preconditioner.col(target).transform([](double value) {
          return std::max(value, gwrs::kNumericalEpsilon);
        });
        factor_type[target] = -1;
      }
    }
  }
};

struct KtddCoefficientWorker : public RcppParallel::Worker {
  const arma::mat& design;
  const arma::umat& neighbor_index;
  const arma::mat& normalized_weights;
  const arma::sp_mat& coefficient_laplacian;
  const double gamma;
  const arma::vec& process;
  const arma::mat& old_coefficients;
  const arma::mat& local_score;
  const arma::cube& factor;
  const arma::ivec& factor_type;
  arma::mat& new_coefficients;

  KtddCoefficientWorker(const arma::mat& design,
                        const arma::umat& neighbor_index,
                        const arma::mat& normalized_weights,
                        const arma::sp_mat& coefficient_laplacian,
                        const double gamma,
                        const arma::vec& process,
                        const arma::mat& old_coefficients,
                        const arma::mat& local_score,
                        const arma::cube& factor,
                        const arma::ivec& factor_type,
                        arma::mat& new_coefficients)
    : design(design), neighbor_index(neighbor_index),
      normalized_weights(normalized_weights),
      coefficient_laplacian(coefficient_laplacian), gamma(gamma),
      process(process), old_coefficients(old_coefficients),
      local_score(local_score), factor(factor), factor_type(factor_type),
      new_coefficients(new_coefficients) {}

  void operator()(std::size_t begin, std::size_t end) {
    const arma::uword k = neighbor_index.n_rows;
    const arma::uword p1 = design.n_cols;
    arma::vec rhs(p1, arma::fill::zeros);
    arma::vec temporary(p1, arma::fill::zeros);

    for (std::size_t target = begin; target < end; ++target) {
      rhs = local_score.col(target);
      for (arma::uword row = 0; row < k; ++row) {
        const arma::uword source = neighbor_index(row, target) - 1;
        rhs -= normalized_weights(row, target) * process[source] *
          design.row(source).t();
      }
      if (gamma > 0.0) {
        for (arma::sp_mat::const_col_iterator entry =
               coefficient_laplacian.begin_col(target);
             entry != coefficient_laplacian.end_col(target); ++entry) {
          if (entry.row() != target) {
            rhs -= gamma * (*entry) *
              old_coefficients.col(entry.row());
          }
        }
      }

      if (factor_type[target] == 0) {
        temporary = arma::solve(arma::trimatl(factor.slice(target)), rhs,
                                arma::solve_opts::fast);
        new_coefficients.col(target) = arma::solve(
          arma::trimatu(factor.slice(target).t()), temporary,
          arma::solve_opts::fast);
      } else {
        new_coefficients.col(target) = factor.slice(target) * rhs;
      }
    }
  }
};

struct KtddCoefficientSystemWorker : public RcppParallel::Worker {
  const arma::mat& design;
  const arma::umat& neighbor_index;
  const arma::mat& normalized_weights;
  const arma::sp_mat& coefficient_laplacian;
  const double gamma;
  const arma::mat& value;
  arma::mat& result;

  KtddCoefficientSystemWorker(
      const arma::mat& design,
      const arma::umat& neighbor_index,
      const arma::mat& normalized_weights,
      const arma::sp_mat& coefficient_laplacian,
      const double gamma,
      const arma::mat& value,
      arma::mat& result)
    : design(design), neighbor_index(neighbor_index),
      normalized_weights(normalized_weights),
      coefficient_laplacian(coefficient_laplacian), gamma(gamma),
      value(value), result(result) {}

  void operator()(std::size_t begin, std::size_t end) {
    const arma::uword k = neighbor_index.n_rows;
    const arma::uword p1 = design.n_cols;
    arma::vec local(p1, arma::fill::zeros);

    for (std::size_t target = begin; target < end; ++target) {
      local.zeros();
      for (arma::uword row = 0; row < k; ++row) {
        const arma::uword source = neighbor_index(row, target) - 1;
        const double projection = arma::dot(
          design.row(source), value.col(target).t());
        local += normalized_weights(row, target) * projection *
          design.row(source).t();
      }
      if (gamma > 0.0) {
        for (arma::sp_mat::const_col_iterator entry =
               coefficient_laplacian.begin_col(target);
             entry != coefficient_laplacian.end_col(target); ++entry) {
          local += gamma * (*entry) * value.col(entry.row());
        }
      }
      result.col(target) = local;
    }
  }
};

struct KtddCoefficientRhsWorker : public RcppParallel::Worker {
  const arma::mat& design;
  const arma::umat& neighbor_index;
  const arma::mat& normalized_weights;
  const arma::vec& process;
  const arma::mat& local_score;
  arma::mat& rhs;

  KtddCoefficientRhsWorker(const arma::mat& design,
                           const arma::umat& neighbor_index,
                           const arma::mat& normalized_weights,
                           const arma::vec& process,
                           const arma::mat& local_score,
                           arma::mat& rhs)
    : design(design), neighbor_index(neighbor_index),
      normalized_weights(normalized_weights), process(process),
      local_score(local_score), rhs(rhs) {}

  void operator()(std::size_t begin, std::size_t end) {
    const arma::uword k = neighbor_index.n_rows;
    for (std::size_t target = begin; target < end; ++target) {
      rhs.col(target) = local_score.col(target);
      for (arma::uword row = 0; row < k; ++row) {
        const arma::uword source = neighbor_index(row, target) - 1;
        rhs.col(target) -= normalized_weights(row, target) * process[source] *
          design.row(source).t();
      }
    }
  }
};

void build_coefficient_rhs(const arma::mat& design,
                           const arma::umat& neighbor_index,
                           const arma::mat& normalized_weights,
                           const arma::vec& process,
                           const arma::mat& local_score,
                           arma::mat& rhs,
                           const int n_threads,
                           const int grain_size) {
  KtddCoefficientRhsWorker worker(
    design, neighbor_index, normalized_weights, process, local_score, rhs);
  if (n_threads == 1 || rhs.n_cols < 2) {
    worker(0, rhs.n_cols);
  } else {
    RcppParallel::parallelFor(
      0, rhs.n_cols, worker,
      static_cast<std::size_t>(std::max(1, grain_size)), n_threads);
  }
}

void apply_coefficient_system(
    const arma::mat& design,
    const arma::umat& neighbor_index,
    const arma::mat& normalized_weights,
    const arma::sp_mat& coefficient_laplacian,
    const double gamma,
    const arma::mat& value,
    arma::mat& result,
    const int n_threads,
    const int grain_size) {
  KtddCoefficientSystemWorker worker(
    design, neighbor_index, normalized_weights, coefficient_laplacian,
    gamma, value, result);
  if (n_threads == 1 || value.n_cols < 2) {
    worker(0, value.n_cols);
  } else {
    RcppParallel::parallelFor(
      0, value.n_cols, worker,
      static_cast<std::size_t>(std::max(1, grain_size)), n_threads);
  }
}

bool solve_coefficients_pcg(
    const arma::mat& design,
    const arma::umat& neighbor_index,
    const arma::mat& normalized_weights,
    const arma::sp_mat& coefficient_laplacian,
    const double gamma,
    const arma::mat& preconditioner,
    const arma::mat& rhs,
    arma::mat& solution,
    const double tolerance,
    const int max_iterations,
    int& iterations,
    double& final_residual,
    const int n_threads,
    const int grain_size) {
  arma::mat system_value(solution.n_rows, solution.n_cols,
                         arma::fill::zeros);
  apply_coefficient_system(
    design, neighbor_index, normalized_weights, coefficient_laplacian,
    gamma, solution, system_value, n_threads, grain_size);
  arma::mat residual = rhs - system_value;
  const double stopping = tolerance * (1.0 + arma::norm(rhs, "fro"));
  final_residual = arma::norm(residual, "fro");
  if (final_residual <= stopping) {
    iterations = 0;
    return true;
  }

  arma::mat scaled_residual = residual / preconditioner;
  arma::mat direction = scaled_residual;
  double residual_inner = arma::accu(residual % scaled_residual);
  for (iterations = 0; iterations < max_iterations; ++iterations) {
    apply_coefficient_system(
      design, neighbor_index, normalized_weights, coefficient_laplacian,
      gamma, direction, system_value, n_threads, grain_size);
    const double denominator = arma::accu(direction % system_value);
    if (!std::isfinite(denominator) ||
        denominator <= std::numeric_limits<double>::min()) break;
    const double step = residual_inner / denominator;
    solution += step * direction;
    residual -= step * system_value;
    final_residual = arma::norm(residual, "fro");
    if (final_residual <= stopping) {
      ++iterations;
      return true;
    }
    scaled_residual = residual / preconditioner;
    const double updated_inner = arma::accu(residual % scaled_residual);
    if (!std::isfinite(updated_inner) ||
        residual_inner <= std::numeric_limits<double>::min()) break;
    direction = scaled_residual +
      (updated_inner / residual_inner) * direction;
    residual_inner = updated_inner;
  }
  return false;
}

struct KtddProcessRhsWorker : public RcppParallel::Worker {
  const arma::mat& design;
  const arma::vec& y;
  const arma::mat& normalized_weights;
  const arma::uvec& reverse_offsets;
  const arma::uvec& reverse_edges;
  const arma::uword neighbors_per_target;
  const arma::mat& coefficients;
  arma::vec& contribution_sum;

  KtddProcessRhsWorker(const arma::mat& design,
                       const arma::vec& y,
                       const arma::mat& normalized_weights,
                       const arma::uvec& reverse_offsets,
                       const arma::uvec& reverse_edges,
                       const arma::uword neighbors_per_target,
                       const arma::mat& coefficients,
                       arma::vec& contribution_sum)
    : design(design), y(y), normalized_weights(normalized_weights),
      reverse_offsets(reverse_offsets), reverse_edges(reverse_edges),
      neighbors_per_target(neighbors_per_target),
      coefficients(coefficients), contribution_sum(contribution_sum) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t source = begin; source < end; ++source) {
      double total = 0.0;
      for (arma::uword position = reverse_offsets[source];
           position < reverse_offsets[source + 1]; ++position) {
        const arma::uword edge = reverse_edges[position];
        const arma::uword target = edge / neighbors_per_target;
        const arma::uword row = edge - target * neighbors_per_target;
        total += normalized_weights(row, target) *
          (y[source] - arma::dot(design.row(source),
                                 coefficients.col(target).t()));
      }
      contribution_sum[source] = total;
    }
  }
};

void build_reverse_neighbor_index(const arma::umat& neighbor_index,
                                  arma::uvec& reverse_offsets,
                                  arma::uvec& reverse_edges) {
  const arma::uword n = neighbor_index.n_cols;
  const arma::uword k = neighbor_index.n_rows;
  reverse_offsets.zeros(n + 1);
  for (arma::uword target = 0; target < n; ++target) {
    for (arma::uword row = 0; row < k; ++row) {
      const arma::uword source = neighbor_index(row, target) - 1;
      ++reverse_offsets[source + 1];
    }
  }
  for (arma::uword source = 0; source < n; ++source) {
    reverse_offsets[source + 1] += reverse_offsets[source];
  }

  reverse_edges.set_size(n * k);
  arma::uvec cursor = reverse_offsets.head(n);
  for (arma::uword target = 0; target < n; ++target) {
    for (arma::uword row = 0; row < k; ++row) {
      const arma::uword source = neighbor_index(row, target) - 1;
      reverse_edges[cursor[source]++] = target * k + row;
    }
  }
}

void build_diffusion_row_index(const arma::sp_mat& diffusion_operator,
                               arma::uvec& row_offsets,
                               arma::uvec& row_columns,
                               arma::vec& row_values) {
  const arma::uword rows = diffusion_operator.n_rows;
  row_offsets.zeros(rows + 1);
  for (arma::sp_mat::const_iterator entry = diffusion_operator.begin();
       entry != diffusion_operator.end(); ++entry) {
    ++row_offsets[entry.row() + 1];
  }
  for (arma::uword row = 0; row < rows; ++row) {
    row_offsets[row + 1] += row_offsets[row];
  }

  row_columns.set_size(diffusion_operator.n_nonzero);
  row_values.set_size(diffusion_operator.n_nonzero);
  arma::uvec cursor = row_offsets.head(rows);
  for (arma::uword column = 0; column < diffusion_operator.n_cols; ++column) {
    for (arma::sp_mat::const_col_iterator entry =
           diffusion_operator.begin_col(column);
         entry != diffusion_operator.end_col(column); ++entry) {
      const arma::uword position = cursor[entry.row()]++;
      row_columns[position] = column;
      row_values[position] = *entry;
    }
  }
}

struct DiffusionForwardWorker : public RcppParallel::Worker {
  const arma::uvec& row_offsets;
  const arma::uvec& row_columns;
  const arma::vec& row_values;
  const arma::vec& value;
  arma::vec& product;

  DiffusionForwardWorker(const arma::uvec& row_offsets,
                         const arma::uvec& row_columns,
                         const arma::vec& row_values,
                         const arma::vec& value,
                         arma::vec& product)
    : row_offsets(row_offsets), row_columns(row_columns),
      row_values(row_values), value(value), product(product) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t row = begin; row < end; ++row) {
      double total = 0.0;
      for (arma::uword position = row_offsets[row];
           position < row_offsets[row + 1]; ++position) {
        total += row_values[position] * value[row_columns[position]];
      }
      product[row] = total;
    }
  }
};

struct DiffusionTransposeWorker : public RcppParallel::Worker {
  const arma::sp_mat& diffusion_operator;
  const arma::vec& observation_weight;
  const arma::vec& value;
  const arma::vec& forward_product;
  const double lambda_diffusion;
  arma::vec& result;

  DiffusionTransposeWorker(const arma::sp_mat& diffusion_operator,
                           const arma::vec& observation_weight,
                           const arma::vec& value,
                           const arma::vec& forward_product,
                           const double lambda_diffusion,
                           arma::vec& result)
    : diffusion_operator(diffusion_operator),
      observation_weight(observation_weight), value(value),
      forward_product(forward_product), lambda_diffusion(lambda_diffusion),
      result(result) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t column = begin; column < end; ++column) {
      double total = 0.0;
      for (arma::sp_mat::const_col_iterator entry =
             diffusion_operator.begin_col(column);
           entry != diffusion_operator.end_col(column); ++entry) {
        total += (*entry) * forward_product[entry.row()];
      }
      result[column] = observation_weight[column] * value[column] +
        lambda_diffusion * total;
    }
  }
};

void apply_process_system(const arma::sp_mat& diffusion_operator,
                          const arma::uvec& diffusion_row_offsets,
                          const arma::uvec& diffusion_row_columns,
                          const arma::vec& diffusion_row_values,
                          const arma::vec& observation_weight,
                          const double lambda_diffusion,
                          const arma::vec& value,
                          arma::vec& forward_product,
                          arma::vec& result,
                          const int n_threads,
                          const int grain_size) {
  if (n_threads == 1 || value.n_elem < 2) {
    result = observation_weight % value + lambda_diffusion *
      diffusion_operator.t() * (diffusion_operator * value);
    return;
  }

  DiffusionForwardWorker forward_worker(
    diffusion_row_offsets, diffusion_row_columns, diffusion_row_values,
    value, forward_product);
  RcppParallel::parallelFor(
    0, diffusion_operator.n_rows, forward_worker,
    static_cast<std::size_t>(std::max(1, grain_size)), n_threads);
  DiffusionTransposeWorker transpose_worker(
    diffusion_operator, observation_weight, value, forward_product,
    lambda_diffusion, result);
  RcppParallel::parallelFor(
    0, diffusion_operator.n_cols, transpose_worker,
    static_cast<std::size_t>(std::max(1, grain_size)), n_threads);
}

bool solve_process_pcg(const arma::sp_mat& diffusion_operator,
                       const arma::uvec& diffusion_row_offsets,
                       const arma::uvec& diffusion_row_columns,
                       const arma::vec& diffusion_row_values,
                       const arma::vec& observation_weight,
                       const arma::vec& preconditioner,
                       const double lambda_diffusion,
                       const arma::vec& rhs,
                       arma::vec& solution,
                       const double tolerance,
                       const int max_iterations,
                       int& iterations,
                       double& final_residual,
                       const int n_threads,
                       const int grain_size) {
  arma::vec forward_product(diffusion_operator.n_rows, arma::fill::zeros);
  arma::vec system_value(solution.n_elem, arma::fill::zeros);
  apply_process_system(
    diffusion_operator, diffusion_row_offsets, diffusion_row_columns,
    diffusion_row_values, observation_weight, lambda_diffusion, solution,
    forward_product, system_value, n_threads, grain_size);
  arma::vec residual = rhs - system_value;
  const double stopping = tolerance * (1.0 + arma::norm(rhs, 2));
  final_residual = arma::norm(residual, 2);
  if (final_residual <= stopping) {
    iterations = 0;
    return true;
  }

  arma::vec scaled_residual = residual / preconditioner;
  arma::vec direction = scaled_residual;
  double residual_inner = arma::dot(residual, scaled_residual);
  for (iterations = 0; iterations < max_iterations; ++iterations) {
    apply_process_system(
      diffusion_operator, diffusion_row_offsets, diffusion_row_columns,
      diffusion_row_values, observation_weight, lambda_diffusion, direction,
      forward_product, system_value, n_threads, grain_size);
    const double denominator = arma::dot(direction, system_value);
    if (!std::isfinite(denominator) ||
        denominator <= std::numeric_limits<double>::min()) break;
    const double step = residual_inner / denominator;
    solution += step * direction;
    residual -= step * system_value;
    final_residual = arma::norm(residual, 2);
    if (final_residual <= stopping) {
      ++iterations;
      return true;
    }
    scaled_residual = residual / preconditioner;
    const double updated_inner = arma::dot(residual, scaled_residual);
    if (!std::isfinite(updated_inner) ||
        residual_inner <= std::numeric_limits<double>::min()) break;
    direction = scaled_residual +
      (updated_inner / residual_inner) * direction;
    residual_inner = updated_inner;
  }
  return false;
}

void validate_ktdd_inputs(const arma::mat& x,
                          const arma::vec& y,
                          const arma::umat& neighbor_index,
                          const arma::mat& neighbor_distance,
                          const arma::sp_mat& diffusion_operator,
                          const arma::vec& source,
                          const arma::sp_mat& coefficient_laplacian) {
  const arma::uword n = x.n_rows;
  if (y.n_elem != n || neighbor_index.n_cols != n ||
      neighbor_distance.n_cols != n ||
      neighbor_index.n_rows != neighbor_distance.n_rows) {
    Rcpp::stop("model and neighbor dimensions are incompatible");
  }
  if (diffusion_operator.n_cols != n ||
      source.n_elem != diffusion_operator.n_rows) {
    Rcpp::stop("diffusion_operator, source, and data are incompatible");
  }
  if (coefficient_laplacian.n_rows != n ||
      coefficient_laplacian.n_cols != n) {
    Rcpp::stop("coefficient_laplacian must be n by n");
  }
  if (!x.is_finite() || !y.is_finite() ||
      !neighbor_distance.is_finite() || !source.is_finite()) {
    Rcpp::stop("dense numerical inputs must be finite");
  }
}

} // anonymous namespace

// [[Rcpp::export]]
int cpp_arma_uword_bits() {
  return static_cast<int>(8 * sizeof(arma::uword));
}

// [[Rcpp::export]]
arma::sp_mat cpp_graph_laplacian(const arma::umat& neighbor_index,
                                 const arma::mat& neighbor_distance,
                                 const int kernel,
                                 const bool adaptive,
                                 const double fixed_bandwidth) {
  if (neighbor_index.n_cols < 1 || neighbor_index.n_rows < 1 ||
      neighbor_index.n_rows != neighbor_distance.n_rows ||
      neighbor_index.n_cols != neighbor_distance.n_cols) {
    Rcpp::stop("neighbor matrices are malformed");
  }
  const arma::uword n = neighbor_index.n_cols;
  const arma::uword k = neighbor_index.n_rows;
  arma::umat locations(2, n * k, arma::fill::zeros);
  arma::vec values(n * k, arma::fill::zeros);
  arma::vec weights(k, arma::fill::zeros);
  arma::uword used = 0;
  for (arma::uword target = 0; target < n; ++target) {
    double bandwidth_used = 0.0;
    gwrs::compute_kernel_weights(neighbor_distance, target, kernel, adaptive,
                                 fixed_bandwidth, weights, bandwidth_used);
    for (arma::uword row = 0; row < k; ++row) {
      const arma::uword encoded = neighbor_index(row, target);
      if (encoded == 0 || encoded > n) {
        Rcpp::stop("neighbor index is outside 1,...,n");
      }
      const arma::uword source = encoded - 1;
      if (source == target || weights[row] <= gwrs::kNumericalEpsilon) continue;
      locations(0, used) = target;
      locations(1, used) = source;
      values[used] = weights[row];
      ++used;
    }
  }
  locations.resize(2, used);
  values.resize(used);
  arma::sp_mat directed(locations, values, n, n, true, true);
  arma::sp_mat adjacency = 0.5 * (directed + directed.t());
  adjacency.diag().zeros();
  arma::vec degree = arma::vec(arma::sum(adjacency, 1));
  arma::sp_mat laplacian = -adjacency;
  laplacian.diag() += degree;
  return laplacian;
}

// [[Rcpp::export]]
Rcpp::List cpp_gwr_ktdd_fit(
    const arma::mat& x,
    const arma::vec& y,
    const arma::umat& neighbor_index,
    const arma::mat& neighbor_distance,
    const int kernel,
    const bool adaptive,
    const double fixed_bandwidth,
    const arma::sp_mat& diffusion_operator,
    const arma::vec& source,
    const arma::sp_mat& coefficient_laplacian,
    const double lambda_diffusion,
    const double gamma,
    const int coefficient_solver,
    const double coefficient_tolerance,
    const int coefficient_max_iterations,
    const double tolerance,
    const int max_iterations,
    const double linear_tolerance,
    const int linear_max_iterations,
    const bool center_process,
    const int n_threads,
    const int grain_size) {
  validate_ktdd_inputs(x, y, neighbor_index, neighbor_distance,
                       diffusion_operator, source, coefficient_laplacian);
  if (kernel < 0 || kernel > 4 ||
      (!adaptive && (!std::isfinite(fixed_bandwidth) ||
                     fixed_bandwidth <= 0.0))) {
    Rcpp::stop("kernel or bandwidth is invalid");
  }
  if (!std::isfinite(lambda_diffusion) || lambda_diffusion <= 0.0 ||
      !std::isfinite(gamma) || gamma < 0.0 ||
      (coefficient_solver != 0 && coefficient_solver != 1) ||
      !std::isfinite(coefficient_tolerance) ||
      coefficient_tolerance <= 0.0 || coefficient_max_iterations < 1 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 ||
      !std::isfinite(linear_tolerance) || linear_tolerance <= 0.0 ||
      max_iterations < 1 || linear_max_iterations < 1) {
    Rcpp::stop("KTDD controls are outside their valid ranges");
  }

  const arma::uword n = x.n_rows;
  const arma::uword p1 = x.n_cols + 1;
  const arma::uword k = neighbor_index.n_rows;
  arma::mat design(n, p1, arma::fill::ones);
  design.cols(1, p1 - 1) = x;
  const bool use_factorized_solver = coefficient_solver == 0;
  arma::mat normalized_weights(k, n, arma::fill::zeros);
  arma::mat local_score(p1, n, arma::fill::zeros);
  arma::cube factor;
  arma::mat coefficient_preconditioner;
  if (use_factorized_solver) {
    factor.zeros(p1, p1, n);
  } else {
    coefficient_preconditioner.zeros(p1, n);
  }
  arma::ivec factor_type(n, arma::fill::value(-1));
  arma::mat coefficients(p1, n, arma::fill::zeros);
  arma::vec weight_sum(n, arma::fill::value(arma::datum::nan));
  arma::vec bandwidth_used(n, arma::fill::value(arma::datum::nan));
  arma::ivec status(n, arma::fill::zeros);

  KtddPrepareWorker prepare(
    design, y, neighbor_index, neighbor_distance, kernel, adaptive,
    fixed_bandwidth, gamma, coefficient_laplacian, use_factorized_solver,
    normalized_weights, local_score, factor, factor_type,
    coefficient_preconditioner, coefficients, weight_sum, bandwidth_used,
    status);
  if (n_threads == 1 || n < 2) {
    prepare(0, n);
  } else {
    RcppParallel::parallelFor(
      0, n, prepare,
      static_cast<std::size_t>(std::max(1, grain_size)), n_threads);
  }
  if (arma::any(status == 3) || arma::any(status == 4)) {
    Rcpp::stop("KTDD preparation failed because of weights or neighbor indices");
  }

  arma::uvec reverse_offsets;
  arma::uvec reverse_edges;
  build_reverse_neighbor_index(neighbor_index, reverse_offsets, reverse_edges);
  arma::uvec diffusion_row_offsets;
  arma::uvec diffusion_row_columns;
  arma::vec diffusion_row_values;
  if (n_threads != 1 && n > 1) {
    build_diffusion_row_index(
      diffusion_operator, diffusion_row_offsets, diffusion_row_columns,
      diffusion_row_values);
  }
  arma::vec observation_weight(n, arma::fill::zeros);
  for (arma::uword source = 0; source < n; ++source) {
    for (arma::uword position = reverse_offsets[source];
         position < reverse_offsets[source + 1]; ++position) {
      const arma::uword edge = reverse_edges[position];
      const arma::uword target = edge / k;
      const arma::uword row = edge - target * k;
      observation_weight[source] += normalized_weights(row, target);
    }
  }
  arma::vec process_diagonal = observation_weight;
  for (arma::sp_mat::const_iterator entry = diffusion_operator.begin();
       entry != diffusion_operator.end(); ++entry) {
    process_diagonal[entry.col()] += lambda_diffusion * (*entry) * (*entry);
  }
  process_diagonal.transform([](double value) {
    return std::max(value, gwrs::kNumericalEpsilon);
  });

  arma::vec process(n, arma::fill::zeros);
  arma::vec updated_process(n, arma::fill::zeros);
  arma::mat updated_coefficients(p1, n, arma::fill::zeros);
  arma::mat coefficient_rhs;
  if (!use_factorized_solver) coefficient_rhs.zeros(p1, n);
  arma::vec contribution_sum(n, arma::fill::zeros);
  const arma::vec physical_rhs = lambda_diffusion *
    diffusion_operator.t() * source;

  bool did_converge = false;
  bool linear_converged = false;
  bool coefficient_linear_converged = use_factorized_solver;
  int iteration = 0;
  int linear_iterations = 0;
  int coefficient_linear_iterations = 0;
  double linear_residual = arma::datum::inf;
  double coefficient_linear_residual = use_factorized_solver ?
    0.0 : arma::datum::inf;
  double coefficient_change = arma::datum::inf;
  double process_change = arma::datum::inf;

  for (; iteration < max_iterations; ++iteration) {
    if (use_factorized_solver) {
      KtddCoefficientWorker coefficient_worker(
        design, neighbor_index, normalized_weights, coefficient_laplacian,
        gamma, process, coefficients, local_score, factor, factor_type,
        updated_coefficients);
      if (n_threads == 1 || n < 2) {
        coefficient_worker(0, n);
      } else {
        RcppParallel::parallelFor(
          0, n, coefficient_worker,
          static_cast<std::size_t>(std::max(1, grain_size)), n_threads);
      }
      coefficient_linear_converged = true;
      coefficient_linear_iterations = 0;
      coefficient_linear_residual = 0.0;
    } else {
      build_coefficient_rhs(
        design, neighbor_index, normalized_weights, process, local_score,
        coefficient_rhs, n_threads, grain_size);
      updated_coefficients = coefficients;
      coefficient_linear_converged = solve_coefficients_pcg(
        design, neighbor_index, normalized_weights, coefficient_laplacian,
        gamma, coefficient_preconditioner, coefficient_rhs,
        updated_coefficients, coefficient_tolerance,
        coefficient_max_iterations, coefficient_linear_iterations,
        coefficient_linear_residual, n_threads, grain_size);
    }

    KtddProcessRhsWorker process_rhs_worker(
      design, y, normalized_weights, reverse_offsets, reverse_edges, k,
      updated_coefficients, contribution_sum);
    if (n_threads == 1 || n < 2) {
      process_rhs_worker(0, n);
    } else {
      RcppParallel::parallelFor(
        0, n, process_rhs_worker,
        static_cast<std::size_t>(std::max(1, grain_size)), n_threads);
    }

    updated_process = process;
    linear_converged = solve_process_pcg(
      diffusion_operator, diffusion_row_offsets, diffusion_row_columns,
      diffusion_row_values, observation_weight, process_diagonal,
      lambda_diffusion, contribution_sum + physical_rhs, updated_process,
      linear_tolerance, linear_max_iterations, linear_iterations,
      linear_residual, n_threads, grain_size);

    if (center_process) {
      const double process_shift = arma::mean(updated_process);
      updated_process -= process_shift;
      updated_coefficients.row(0) += process_shift;
    }

    coefficient_change = arma::abs(updated_coefficients - coefficients).max() /
      (1.0 + arma::abs(coefficients).max());
    process_change = arma::abs(updated_process - process).max() /
      (1.0 + arma::abs(process).max());
    coefficients.swap(updated_coefficients);
    process.swap(updated_process);
    if (coefficient_linear_converged && linear_converged &&
        std::max(coefficient_change, process_change) <= tolerance) {
      did_converge = true;
      ++iteration;
      break;
    }
  }

  arma::vec fitted(n, arma::fill::zeros);
  for (arma::uword target = 0; target < n; ++target) {
    fitted[target] = arma::dot(design.row(target),
                               coefficients.col(target).t()) +
      process[target];
  }
  const arma::vec residuals = y - fitted;

  double local_loss = 0.0;
  for (arma::uword target = 0; target < n; ++target) {
    for (arma::uword row = 0; row < k; ++row) {
      const arma::uword source_index = neighbor_index(row, target) - 1;
      const double local_residual = y[source_index] - process[source_index] -
        arma::dot(design.row(source_index), coefficients.col(target).t());
      local_loss += 0.5 * normalized_weights(row, target) *
        local_residual * local_residual;
    }
  }
  const arma::vec physical_residual = diffusion_operator * process - source;
  const double diffusion_loss = 0.5 * lambda_diffusion *
    arma::dot(physical_residual, physical_residual);
  double smoothness_loss = 0.0;
  for (arma::uword column = 0; column < p1; ++column) {
    const arma::vec surface = coefficients.row(column).t();
    smoothness_loss += 0.5 * gamma *
      arma::dot(surface, coefficient_laplacian * surface);
  }

  return Rcpp::List::create(
    Rcpp::_["coefficients"] = coefficients,
    Rcpp::_["process"] = process,
    Rcpp::_["fitted"] = fitted,
    Rcpp::_["residuals"] = residuals,
    Rcpp::_["objective"] = local_loss + diffusion_loss + smoothness_loss,
    Rcpp::_["loss_components"] = Rcpp::NumericVector::create(
      Rcpp::_["local"] = local_loss,
      Rcpp::_["diffusion"] = diffusion_loss,
      Rcpp::_["coefficient_smoothness"] = smoothness_loss),
    Rcpp::_["iterations"] = iteration,
    Rcpp::_["converged"] = did_converge,
    Rcpp::_["coefficient_change"] = coefficient_change,
    Rcpp::_["process_change"] = process_change,
    Rcpp::_["linear_converged"] = linear_converged,
    Rcpp::_["linear_iterations"] = linear_iterations,
    Rcpp::_["linear_residual"] = linear_residual,
    Rcpp::_["coefficient_linear_converged"] =
      coefficient_linear_converged,
    Rcpp::_["coefficient_linear_iterations"] =
      coefficient_linear_iterations,
    Rcpp::_["coefficient_linear_residual"] =
      coefficient_linear_residual,
    Rcpp::_["coefficient_solver"] = coefficient_solver,
    Rcpp::_["process_centered"] = center_process,
    Rcpp::_["factor_type"] = factor_type,
    Rcpp::_["status"] = status,
    Rcpp::_["weight_sum"] = weight_sum,
    Rcpp::_["bandwidth"] = bandwidth_used,
    Rcpp::_["observation_weight"] = observation_weight
  );
}
