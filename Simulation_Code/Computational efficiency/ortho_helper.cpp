#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Helper to get residuals from OLS: y = X*beta + e
arma::vec get_ols_residuals(const arma::mat& X, const arma::vec& y) {
    arma::vec beta = arma::solve(X, y);
    return y - X * beta;
}

// [[Rcpp::export]]
arma::mat create_structured_orthogonal_terms_cpp(arma::mat main_effects) {
    int N = main_effects.n_rows;
    int m = main_effects.n_cols;
    
    if (m == 0) return arma::mat(N, 0);

    // 1. Center the main effects
    arma::mat X_centered = main_effects;
    arma::rowvec means = arma::mean(X_centered, 0);
    X_centered.each_row() -= means;

    // 2. Quadratic Terms Orthogonalization
    arma::mat quad_terms_ortho(N, m);
    arma::vec ones = arma::ones<vec>(N);

    for (int i = 0; i < m; ++i) {
        arma::vec raw_squared = square(main_effects.col(i));
        
        // Basis: [1, linear_i]
        arma::mat basis(N, 2);
        basis.col(0) = ones;
        basis.col(1) = X_centered.col(i);
        
        quad_terms_ortho.col(i) = get_ols_residuals(basis, raw_squared);
    }

    // 3. Interaction Terms Orthogonalization
    int n_int = (m * (m - 1)) / 2;
    if (n_int == 0) {
        return quad_terms_ortho;
    }

    arma::mat interaction_terms_ortho(N, n_int);
    int current_col_idx = 0;

    for (int i = 0; i < m - 1; ++i) {
        for (int j = i + 1; j < m; ++j) {
            arma::vec raw_interaction = main_effects.col(i) % main_effects.col(j);
            
            // Basis: [1, linear_i, linear_j, quad_i, quad_j]
            arma::mat basis(N, 5);
            basis.col(0) = ones;
            basis.col(1) = X_centered.col(i);
            basis.col(2) = X_centered.col(j);
            basis.col(3) = quad_terms_ortho.col(i);
            basis.col(4) = quad_terms_ortho.col(j);
            
            interaction_terms_ortho.col(current_col_idx) = get_ols_residuals(basis, raw_interaction);
            current_col_idx++;
        }
    }

    // Combine output: [Quad, Inter]
    return join_horiz(quad_terms_ortho, interaction_terms_ortho);
}
