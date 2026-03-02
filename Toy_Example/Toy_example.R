# HOSS Toy Example: Blind vs Informed and Raw vs Ortho
# Symbols match the paper (X, N, n_sub)
rm(list=ls())
library(ggplot2)
library(gridExtra)
library(Rcpp)
library(RcppArmadillo)
library(MASS)

# ==============================================================================
# 1. C++ Integration (OSS / Leverage)
# ==============================================================================
sourceCpp(code='
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::vec lev(const arma::mat& X) {
  arma::mat U, V;
  arma::vec s;
  arma::svd_econ(U, s, V, X, "left");
  arma::vec row_crossprods(U.n_rows);
  for (arma::uword i = 0; i < U.n_rows; ++i) {
    arma::rowvec row = U.row(i);
    row_crossprods(i) = arma::dot(row, row);
  }
  return row_crossprods / X.n_cols;
}

// [[Rcpp::export]]
arma::vec bottom_k(arma::vec x, unsigned int k) {
  arma::vec x2 = x; 
  arma::vec ind(k); 
  std::nth_element(x.begin(), x.begin() + k - 1, x.end()); 
  for(int ii=0, i=0; i<int(x.n_elem) && ii<int(k); i++){
    if(x2[i] <= x[k-1])  ind[ii++] = i; 
  }
  return ind;
}

// [[Rcpp::export]]
arma::vec top_k(arma::vec x, unsigned int k) {
  return bottom_k(-x,k);
}
  
// [[Rcpp::export]]
arma::vec Dscr_cpp(arma::mat X, arma::vec xa, arma::mat y, double ya, double tPow) {
  int n=X.n_rows;
  int p=X.n_cols;
  arma::vec B = zeros<vec>(n);
  for(int i=0; i<n; i++){
    B(i) = pow(accu(X.row(i)==y)+p-xa(i)/2-ya/2,tPow); 
  }
  return B;
}

// [[Rcpp::export]]
arma::uvec OAJ2_cpp(arma::mat x, int k, double tPow=2){
  int n=x.n_rows;
  arma::uvec candi=linspace<uvec>(1,n,n);
  arma::uvec ind=linspace<uvec>(1,k,k);
  arma::vec L=sum(pow(x,2),1);
  arma::vec xa=L;
  uword mm=L.index_max();
  ind(0)=candi(mm);
  candi.shed_row(mm);
  L.shed_row(mm);
  
  arma::mat sx=sign(x);
  double r=log(n/k)/log(k);
  for(int i=1; i<k; i++){
    if(i==1)
      L=Dscr_cpp(sx.rows(candi-1),xa.elem(candi-1),sx.row(ind(i-1)-1),xa(ind(i-1)-1),tPow);
    else
      L=L+Dscr_cpp(sx.rows(candi-1),xa.elem(candi-1),sx.row(ind(i-1)-1),xa(ind(i-1)-1),tPow);
    
    mm=L.index_min();
    ind(i)=candi(mm);
    candi.shed_row(mm);
    L.shed_row(mm);
    
    int nc=floor(n/pow(i,r));
    if((i>1) & (L.n_elem>double(nc))){
      arma::uvec tt=arma::conv_to<arma::uvec>::from(bottom_k(L,nc));
      L=L.elem(tt);
      candi=candi.elem(tt);
    }
  }
  return ind;
}
')

# ==============================================================================
# 2. R Helper Functions
# ==============================================================================

iboss=function(x,k){
  ind=NULL
  m = ncol(x)
  r = rep(floor(k/2/m),m)
  if(sum(r)<k/2) r[1:((k-2*sum(r))/2)] = r[1:((k-2*sum(r))/2)]+1
  candi=1:nrow(x)
  for (i in 1:m)
  {
    xi = x[candi,i]
    j1 = top_k(xi,r[i])+1
    j2 = bottom_k(xi,r[i])+1
    j = unique(c(j1,j2))
    if(length(j)<2*r[i]) {jj=(1:length(candi))[-j];j=c(j,jj[1:(2*r[i]-length(j))])}
    ind = c(ind,candi[j])
    candi=setdiff(candi,ind)
  }
  return(ind)
}

scalex=function(a){
  if(max(a)==min(a)) return(rep(0,length(a)))
  2*(a-min(a))/(max(a)-min(a))-1
}

standardize_matrix <- function(X) {
  X <- as.matrix(X)
  apply(X, 2, function(col) {
    if (sd(col) == 0) return(col)
    return((col - mean(col)) / sd(col))
  })
}

# Calculate Log-Determinant of Information Matrix
calc_score <- function(indices, X_ground_truth) {
  X_sub <- cbind(1, X_ground_truth)[indices, ] 
  info_mat <- t(X_sub) %*% X_sub
  val <- det(info_mat)
  if(val <= 1e-12) return(-999) 
  return(log(val))
}

calc_kappa <- function(indices, X_ground_truth) {
  X_sub <- cbind(1, X_ground_truth)[indices, ]
  kappa(X_sub, exact = FALSE)
}

run_lev_method <- function(X_design, X_eval, N, n_sub, n_loop=100) {
  lev_probs <- lev(standardize_matrix(X_design))
  scores_lev <- numeric(n_loop)
  
  for(i in 1:n_loop) {
    idx_l <- sample(N, n_sub, prob=lev_probs)
    scores_lev[i] <- calc_score(idx_l, X_eval)
  }
  return(mean(scores_lev))
}

run_lev_method_kappa <- function(X_design, X_eval, N, n_sub, n_loop=100) {
  lev_probs <- lev(standardize_matrix(X_design))
  scores_lev <- numeric(n_loop)
  kappas_lev <- numeric(n_loop)
  for(i in 1:n_loop) {
    idx_l <- sample(N, n_sub, prob=lev_probs)
    scores_lev[i] <- calc_score(idx_l, X_eval)
    kappas_lev[i] <- calc_kappa(idx_l, X_eval)
  }
  return(list(score=mean(scores_lev), kappa=mean(kappas_lev)))
}

plot_result <- function(x1, x2, indices, method_name, score, kappa_val = NA_real_) {
  df <- data.frame(x1=x1, x2=x2, Type="Unselected")
  df$Type[indices] <- "Selected"
  df <- df[order(df$Type), ] 
  label_text <- if (is.na(kappa_val)) {
    paste0("D-criterion: ", round(score, 2))
  } else {
    paste0("D-criterion: ", round(score, 2), "\nCondition number: ", round(kappa_val, 1))
  }
  
  ggplot(df, aes(x=x1, y=x2, color=Type, size=Type)) +
    geom_point(alpha=0.6) +
    scale_color_manual(values=c("Red", "Grey90")) +
    scale_size_manual(values=c(2.5, 0.8)) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size=11, face="bold"),
          axis.title = element_blank()) +
    labs(title = method_name, subtitle = label_text)
}

# ==============================================================================
# 3. Main Simulation
# ==============================================================================

run_simulation <- function() {
  
  N <- 1000
  n_sub <- 40
  
  # --------------------------------------------------------------------------
  # PART 1: Blind vs Informed (Trade-off)
  # --------------------------------------------------------------------------
  set.seed(123)
  x1 <- runif(N, -1, 1)
  x2 <- runif(N, -1, 1)
  
  X_blind <- cbind(x1, x2)
  X_informed <- cbind(x1, x2, x1^2, x2^2, x1*x2)
  
  # Deterministic Sampling (IBOSS, OSS)
  idx_blind_iboss <- iboss(X_blind, n_sub)
  idx_blind_oss   <- OAJ2_cpp(apply(X_blind, 2, scalex), n_sub)
  
  idx_info_iboss  <- iboss(X_informed, n_sub)
  idx_info_oss    <- OAJ2_cpp(apply(X_informed, 2, scalex), n_sub)
  
  # Stochastic Sampling (LEV, 100 loops)
  lev_blind_simple  <- run_lev_method(X_blind, X_blind, N, n_sub)
  lev_blind_complex <- run_lev_method(X_blind, X_informed, N, n_sub)
  lev_info_simple   <- run_lev_method(X_informed, X_blind, N, n_sub)
  lev_info_complex  <- run_lev_method(X_informed, X_informed, N, n_sub)
  
  # Calculate Scores
  res_iboss_blind <- c(calc_score(idx_blind_iboss, X_blind), calc_score(idx_blind_iboss, X_informed))
  res_oss_blind   <- c(calc_score(idx_blind_oss, X_blind), calc_score(idx_blind_oss, X_informed))
  res_iboss_info  <- c(calc_score(idx_info_iboss, X_blind), calc_score(idx_info_iboss, X_informed))
  res_oss_info    <- c(calc_score(idx_info_oss, X_blind), calc_score(idx_info_oss, X_informed))
  
  df_part1 <- data.frame(
    Method = c("IBOSS", "OSS", "LEV (Avg)"),
    Blind_on_Simple = c(res_iboss_blind[1], res_oss_blind[1], lev_blind_simple),
    Blind_on_Complex = c(res_iboss_blind[2], res_oss_blind[2], lev_blind_complex),
    Informed_on_Simple = c(res_iboss_info[1], res_oss_info[1], lev_info_simple),
    Informed_on_Complex = c(res_iboss_info[2], res_oss_info[2], lev_info_complex)
  )
  
  # --------------------------------------------------------------------------
  # PART 2: Raw vs Ortho (Strategy)
  # Setting: Non-centered data [1, 3] creates correlation between X and X^2
  # --------------------------------------------------------------------------
  set.seed(999)
  x1_s <- runif(N, 1, 3) 
  x2_s <- runif(N, 1, 3)
  
  # Raw Design Matrix
  X_raw <- cbind(x1_s, x2_s, x1_s^2, x2_s^2, x1_s*x2_s)
  
  # Orthogonalized Design Matrix (Proposed Strategy)
  r_x1sq <- lm(x1_s^2 ~ x1_s)$residuals
  r_x2sq <- lm(x2_s^2 ~ x2_s)$residuals
  r_int  <- lm(x1_s*x2_s ~ x1_s + x2_s + r_x1sq + r_x2sq)$residuals
  X_ortho <- cbind(x1_s, x2_s, r_x1sq, r_x2sq, r_int)
  
  # Evaluation ground truth is always X_raw
  X_eval_p2 <- X_raw 
  
  # Deterministic Sampling (IBOSS, OSS)
  idx_raw_iboss <- iboss(X_raw, n_sub)
  idx_raw_oss   <- OAJ2_cpp(apply(X_raw, 2, scalex), n_sub)
  
  idx_orth_iboss <- iboss(X_ortho, n_sub)
  idx_orth_oss   <- OAJ2_cpp(apply(X_ortho, 2, scalex), n_sub)
  
  get_metrics <- function(idx, X_true) {
    c(calc_score(idx, X_true), calc_kappa(idx, X_true))
  }
  
  lev_raw_p2  <- run_lev_method_kappa(X_raw, X_eval_p2, N, n_sub)
  lev_orth_p2 <- run_lev_method_kappa(X_ortho, X_eval_p2, N, n_sub)
  
  m_raw_iboss <- get_metrics(idx_raw_iboss, X_eval_p2)
  m_raw_oss   <- get_metrics(idx_raw_oss, X_eval_p2)
  m_orth_iboss <- get_metrics(idx_orth_iboss, X_eval_p2)
  m_orth_oss   <- get_metrics(idx_orth_oss, X_eval_p2)
  
  df_part2 <- data.frame(
    Method = c("IBOSS", "OSS", "LEV (Avg)"),
    Raw_Score = c(m_raw_iboss[1], m_raw_oss[1], lev_raw_p2$score),
    Ortho_Score = c(m_orth_iboss[1], m_orth_oss[1], lev_orth_p2$score),
    Raw_Kappa = c(m_raw_iboss[2], m_raw_oss[2], lev_raw_p2$kappa),
    Ortho_Kappa = c(m_orth_iboss[2], m_orth_oss[2], lev_orth_p2$kappa)
  )
  
  # Plots for Part 2 (IBOSS and OSS only)
  p2_list <- list()
  p2_list[[1]] <- plot_result(x1_s, x2_s, idx_raw_iboss,  "IBOSS (Raw)", m_raw_iboss[1], m_raw_iboss[2])
  p2_list[[2]] <- plot_result(x1_s, x2_s, idx_orth_iboss, "IBOSS (Orthogonalized)", m_orth_iboss[1], m_orth_iboss[2])
  p2_list[[3]] <- plot_result(x1_s, x2_s, idx_raw_oss,    "OSS (Raw)", m_raw_oss[1], m_raw_oss[2])
  p2_list[[4]] <- plot_result(x1_s, x2_s, idx_orth_oss,   "OSS (Orthogonalized)", m_orth_oss[1], m_orth_oss[2])
  
  return(list(t1=df_part1, t2=df_part2, p2=p2_list))
}

# ==============================================================================
# 4. Execution and Output
# ==============================================================================

res <- run_simulation()

cat("\n--- Part 1: Blind vs Informed ---\n")
print(res$t1, digits=4)

cat("\n--- Part 2: Raw vs Ortho ---\n")
print(res$t2, digits=4)

# Create Output PNG
png("simulation_result.png", width = 12, height = 10, units = "in", res = 300)
grid.arrange(grobs=res$p2, ncol=2)
dev.off()
