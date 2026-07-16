#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

namespace {

struct GroupStats {
  arma::uvec count;
  arma::mat sum;
  arma::cube cross;
};

inline double corr_coeff(const double lambda, const unsigned int n) {
  return lambda / (1.0 + lambda * static_cast<double>(n));
}

arma::uvec checked_indices(const Rcpp::IntegerVector& selected,
                           const arma::uword n_rows) {
  if (selected.size() == 0) {
    Rcpp::stop("selected must contain at least one row");
  }
  arma::uvec out(selected.size());
  std::vector<unsigned char> seen(n_rows, 0);
  for (R_xlen_t i = 0; i < selected.size(); ++i) {
    const int value = selected[i];
    if (value < 1 || static_cast<arma::uword>(value) > n_rows) {
      Rcpp::stop("selected contains an out-of-range 1-based row index");
    }
    const arma::uword idx = static_cast<arma::uword>(value - 1);
    if (seen[idx]) {
      Rcpp::stop("selected contains duplicate row indices");
    }
    seen[idx] = 1;
    out[i] = idx;
  }
  return out;
}

arma::uvec checked_groups(const Rcpp::IntegerVector& group,
                          const arma::uword n_rows,
                          arma::uword& n_groups) {
  if (static_cast<arma::uword>(group.size()) != n_rows) {
    Rcpp::stop("group must have one entry per row of X");
  }
  arma::uvec out(n_rows);
  int max_group = 0;
  for (arma::uword i = 0; i < n_rows; ++i) {
    const int value = group[i];
    if (value < 1) {
      Rcpp::stop("group labels must be positive 1-based integers");
    }
    out[i] = static_cast<arma::uword>(value - 1);
    if (value > max_group) max_group = value;
  }
  n_groups = static_cast<arma::uword>(max_group);
  return out;
}

GroupStats build_stats(const arma::mat& X,
                       const arma::uvec& groups,
                       const arma::uvec& selected,
                       const arma::uword n_groups) {
  const arma::uword p = X.n_cols;
  GroupStats stats{
    arma::uvec(n_groups, arma::fill::zeros),
    arma::mat(n_groups, p, arma::fill::zeros),
    arma::cube(p, p, n_groups, arma::fill::zeros)
  };
  for (arma::uword j = 0; j < selected.n_elem; ++j) {
    const arma::uword row = selected[j];
    const arma::uword g = groups[row];
    const arma::rowvec x = X.row(row);
    stats.count[g] += 1;
    stats.sum.row(g) += x;
    stats.cross.slice(g) += x.t() * x;
  }
  return stats;
}

arma::mat info_from_stats(const GroupStats& stats, const double lambda) {
  const arma::uword p = stats.sum.n_cols;
  arma::mat info(p, p, arma::fill::zeros);
  for (arma::uword g = 0; g < stats.count.n_elem; ++g) {
    if (stats.count[g] == 0) continue;
    info += stats.cross.slice(g);
    const double coeff = corr_coeff(lambda, stats.count[g]);
    if (coeff > 0.0) {
      const arma::vec s = stats.sum.row(g).t();
      info -= coeff * (s * s.t());
    }
  }
  return arma::symmatu(info);
}

bool safe_logdet(const arma::mat& input, double& value, double& jitter_used) {
  arma::mat info = arma::symmatu(input);
  const double base = std::max(1.0, std::abs(arma::trace(info)) /
                                      static_cast<double>(info.n_rows));
  jitter_used = 0.0;
  for (int attempt = 0; attempt < 8; ++attempt) {
    arma::mat chol_upper;
    if (arma::chol(chol_upper, info)) {
      value = 2.0 * arma::accu(arma::log(chol_upper.diag()));
      return std::isfinite(value);
    }
    jitter_used = (attempt == 0 ? 1e-12 * base : jitter_used * 10.0);
    info = arma::symmatu(input) + jitter_used * arma::eye(input.n_rows, input.n_cols);
  }
  value = -std::numeric_limits<double>::infinity();
  return false;
}

inline void add_term(arma::mat& U,
                     arma::vec& coeff,
                     int& k,
                     const arma::rowvec& vector,
                     const double weight) {
  if (std::abs(weight) <= 1e-15) return;
  U.col(k) = vector.t();
  coeff[k] = weight;
  ++k;
}

int fill_swap_terms(arma::mat& U,
                    arma::vec& coeff,
                    const arma::mat& X,
                    const arma::uvec& groups,
                    const GroupStats& stats,
                    const arma::uword row_out,
                    const arma::uword row_in,
                    const double lambda) {
  int k = 0;
  const arma::uword group_out = groups[row_out];
  const arma::uword group_in = groups[row_in];
  const arma::rowvec x_out = X.row(row_out);
  const arma::rowvec x_in = X.row(row_in);

  add_term(U, coeff, k, x_out, -1.0);
  add_term(U, coeff, k, x_in, 1.0);

  if (group_out == group_in) {
    const arma::rowvec old_sum = stats.sum.row(group_out);
    const arma::rowvec new_sum = old_sum - x_out + x_in;
    const double c = corr_coeff(lambda, stats.count[group_out]);
    add_term(U, coeff, k, old_sum, c);
    add_term(U, coeff, k, new_sum, -c);
    return k;
  }

  const arma::rowvec old_sum_out = stats.sum.row(group_out);
  const arma::rowvec new_sum_out = old_sum_out - x_out;
  const double c_old_out = corr_coeff(lambda, stats.count[group_out]);
  const double c_new_out = corr_coeff(lambda, stats.count[group_out] - 1);
  add_term(U, coeff, k, old_sum_out, c_old_out);
  add_term(U, coeff, k, new_sum_out, -c_new_out);

  const arma::rowvec old_sum_in = stats.sum.row(group_in);
  const arma::rowvec new_sum_in = old_sum_in + x_in;
  const double c_old_in = corr_coeff(lambda, stats.count[group_in]);
  const double c_new_in = corr_coeff(lambda, stats.count[group_in] + 1);
  add_term(U, coeff, k, old_sum_in, c_old_in);
  add_term(U, coeff, k, new_sum_in, -c_new_in);
  return k;
}

void apply_swap(GroupStats& stats,
                const arma::mat& X,
                const arma::uvec& groups,
                const arma::uword row_out,
                const arma::uword row_in) {
  const arma::uword group_out = groups[row_out];
  const arma::uword group_in = groups[row_in];
  const arma::rowvec x_out = X.row(row_out);
  const arma::rowvec x_in = X.row(row_in);

  stats.count[group_out] -= 1;
  stats.sum.row(group_out) -= x_out;
  stats.cross.slice(group_out) -= x_out.t() * x_out;

  stats.count[group_in] += 1;
  stats.sum.row(group_in) += x_in;
  stats.cross.slice(group_in) += x_in.t() * x_in;
}

}  // namespace

// Random-intercept information matrix, apart from the common 1/sigma_e^2
// factor. lambda is sigma_a^2 / sigma_e^2. Row indices are 1-based.
// [[Rcpp::export]]
arma::mat weighted_info_cpp(const arma::mat& X,
                            const Rcpp::IntegerVector& group,
                            const Rcpp::IntegerVector& selected,
                            const double lambda) {
  if (!std::isfinite(lambda) || lambda < 0.0) {
    Rcpp::stop("lambda must be finite and non-negative");
  }
  arma::uword n_groups = 0;
  const arma::uvec groups = checked_groups(group, X.n_rows, n_groups);
  const arma::uvec selected_zero = checked_indices(selected, X.n_rows);
  const GroupStats stats = build_stats(X, groups, selected_zero, n_groups);
  return info_from_stats(stats, lambda);
}

// [[Rcpp::export]]
Rcpp::List weighted_logdet_cpp(const arma::mat& X,
                               const Rcpp::IntegerVector& group,
                               const Rcpp::IntegerVector& selected,
                               const double lambda) {
  const arma::mat info = weighted_info_cpp(X, group, selected, lambda);
  double value = 0.0;
  double jitter = 0.0;
  const bool ok = safe_logdet(info, value, jitter);
  return Rcpp::List::create(
    Rcpp::Named("logdet") = value,
    Rcpp::Named("ok") = ok,
    Rcpp::Named("jitter") = jitter
  );
}

// Exhaustive best-swap exchange for the locally weighted LMM criterion.
// The result is a one-row-exchange local optimum when converged is TRUE.
// [[Rcpp::export]]
Rcpp::List weighted_exchange_cpp(const arma::mat& X,
                                 const Rcpp::IntegerVector& group,
                                 const Rcpp::IntegerVector& initial_selected,
                                 const double lambda,
                                 const int max_exchanges = 40,
                                 const double tolerance = 1e-9) {
  if (!std::isfinite(lambda) || lambda < 0.0) {
    Rcpp::stop("lambda must be finite and non-negative");
  }
  if (max_exchanges < 1) Rcpp::stop("max_exchanges must be at least 1");
  if (!std::isfinite(tolerance) || tolerance <= 0.0) {
    Rcpp::stop("tolerance must be finite and positive");
  }

  arma::uword n_groups = 0;
  const arma::uvec groups = checked_groups(group, X.n_rows, n_groups);
  arma::uvec selected = checked_indices(initial_selected, X.n_rows);
  GroupStats stats = build_stats(X, groups, selected, n_groups);
  arma::mat info = info_from_stats(stats, lambda);
  double current_logdet = 0.0;
  double current_jitter = 0.0;
  if (!safe_logdet(info, current_logdet, current_jitter)) {
    Rcpp::stop("initial design has a singular weighted information matrix");
  }

  std::vector<unsigned char> is_selected(X.n_rows, 0);
  for (arma::uword j = 0; j < selected.n_elem; ++j) {
    is_selected[selected[j]] = 1;
  }

  const arma::uword p = X.n_cols;
  arma::mat U(p, 6, arma::fill::zeros);
  arma::vec coeff(6, arma::fill::zeros);
  int exchanges = 0;
  bool converged = false;
  bool numerical_stop = false;
  double last_gain = 0.0;

  for (int iteration = 0; iteration < max_exchanges; ++iteration) {
    arma::mat info_inverse;
    if (!arma::inv_sympd(info_inverse, info)) {
      const double scale = std::max(1.0, std::abs(arma::trace(info)) /
                                            static_cast<double>(p));
      const arma::mat stabilized = info + 1e-10 * scale * arma::eye(p, p);
      if (!arma::inv_sympd(info_inverse, stabilized)) {
        Rcpp::stop("could not invert the current weighted information matrix");
      }
    }

    double best_gain = tolerance;
    arma::uword best_position = 0;
    arma::uword best_row_in = 0;
    bool found = false;

    for (arma::uword pos = 0; pos < selected.n_elem; ++pos) {
      const arma::uword row_out = selected[pos];
      for (arma::uword row_in = 0; row_in < X.n_rows; ++row_in) {
        if (is_selected[row_in]) continue;
        const int k = fill_swap_terms(U, coeff, X, groups, stats,
                                      row_out, row_in, lambda);
        if (k == 0) continue;
        const arma::mat U_k = U.cols(0, static_cast<arma::uword>(k - 1));
        const arma::vec d_k = coeff.head(static_cast<arma::uword>(k));
        const arma::mat gram = U_k.t() * info_inverse * U_k;
        const arma::mat small = arma::eye(k, k) + arma::diagmat(d_k) * gram;
        double gain = 0.0;
        double sign = 0.0;
        arma::log_det(gain, sign, small);
        if (sign > 0.0 && std::isfinite(gain) && gain > best_gain) {
          best_gain = gain;
          best_position = pos;
          best_row_in = row_in;
          found = true;
        }
      }
    }

    if (!found) {
      converged = true;
      break;
    }

    const arma::uword row_out = selected[best_position];
    const int k = fill_swap_terms(U, coeff, X, groups, stats,
                                  row_out, best_row_in, lambda);
    const arma::mat U_k = U.cols(0, static_cast<arma::uword>(k - 1));
    const arma::vec d_k = coeff.head(static_cast<arma::uword>(k));
    const arma::mat proposed_info = arma::symmatu(
      info + U_k * arma::diagmat(d_k) * U_k.t()
    );
    double proposed_logdet = 0.0;
    double proposed_jitter = 0.0;
    if (!safe_logdet(proposed_info, proposed_logdet, proposed_jitter) ||
        proposed_logdet - current_logdet <= tolerance) {
      numerical_stop = true;
      break;
    }

    apply_swap(stats, X, groups, row_out, best_row_in);
    is_selected[row_out] = 0;
    is_selected[best_row_in] = 1;
    selected[best_position] = best_row_in;
    info = info_from_stats(stats, lambda);
    double rebuilt_jitter = 0.0;
    double rebuilt_logdet = 0.0;
    if (!safe_logdet(info, rebuilt_logdet, rebuilt_jitter)) {
      Rcpp::stop("weighted information matrix became singular after an exchange");
    }
    last_gain = rebuilt_logdet - current_logdet;
    current_logdet = rebuilt_logdet;
    current_jitter = rebuilt_jitter;
    ++exchanges;
  }

  arma::uvec selected_sorted = arma::sort(selected) + 1;
  Rcpp::IntegerVector selected_out(selected_sorted.n_elem);
  for (arma::uword i = 0; i < selected_sorted.n_elem; ++i) {
    selected_out[i] = static_cast<int>(selected_sorted[i]);
  }
  return Rcpp::List::create(
    Rcpp::Named("selected") = selected_out,
    Rcpp::Named("logdet") = current_logdet,
    Rcpp::Named("exchanges") = exchanges,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("numerical_stop") = numerical_stop,
    Rcpp::Named("last_gain") = last_gain,
    Rcpp::Named("jitter") = current_jitter
  );
}
