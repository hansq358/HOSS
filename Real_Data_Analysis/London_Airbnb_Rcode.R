# HOSS Real Data Analysis: London Airbnb
# ==============================================================================
# SIMULATION (N=1000) & ANALYSIS
# ==============================================================================
rm(list=ls())
library(MASS)
library(caret)
library(AlgDesign)
library(lme4)         
library(dplyr)       
library(Rcpp)
library(RcppArmadillo)
library(stringr)
library(ggplot2)

# --- 1. Embed C++ Core Algorithms ---
sourceCpp(code='
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::vec lev(const arma::mat& X) {
  if (X.n_rows == 0) return arma::vec();
  arma::mat U, V;
  arma::vec s;
  // Use "std" instead of "econ" if n < p, but here n is usually > p. 
  // Safety check for svd failure
  bool success = arma::svd_econ(U, s, V, X, "left");
  if (!success) return arma::ones<arma::vec>(X.n_rows) / X.n_rows;
  
  arma::vec row_crossprods(U.n_rows);
  for (arma::uword i = 0; i < U.n_rows; ++i) {
    row_crossprods(i) = arma::dot(U.row(i), U.row(i));
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

# --- Helper Functions ---
iboss <- function(x, k) {
  ind <- NULL
  m <- ncol(x)
  r <- rep(floor(k/2/m), m)
  if(sum(r) < k/2) r[1:((k - 2*sum(r))/2)] <- r[1:((k - 2*sum(r))/2)] + 1
  candi <- 1:nrow(x)
  for (i in 1:m) {
    xi <- x[candi, i]
    j1 <- top_k(xi, r[i]) + 1
    j2 <- bottom_k(xi, r[i]) + 1
    j <- unique(c(j1, j2))
    if(length(j) < 2*r[i]) {
      jj <- (1:length(candi))[-j]
      j <- c(j, jj[1:(2*r[i] - length(j))])
    }
    ind <- c(ind, candi[j])
    candi <- setdiff(candi, ind)
  }
  return(ind)
}

scalex <- function(a) {
  if(max(a) == min(a)) return(rep(0, length(a)))
  2 * (a - min(a)) / (max(a) - min(a)) - 1
}

safe_scale <- function(x) {
  apply(x, 2, function(col) {
    if (sd(col) == 0) return(rep(0, length(col)))
    return(scale(col))
  })
}

create_scaled_poly_contrasts <- function(factor_column, factor_name) {
  if (!is.factor(factor_column)) {
    factor_column <- factor(as.character(factor_column))
  }
  num_levels <- nlevels(factor_column)
  if (num_levels <= 1) {
    return(matrix(nrow = length(factor_column), ncol = 0))
  }
  contrasts_matrix <- model.matrix(~ factor_column, 
                                   contrasts.arg = list(factor_column = contr.poly(num_levels)))[, -1, drop = FALSE]
  colnames(contrasts_matrix) <- paste(factor_name, colnames(contrasts_matrix), sep = ".")
  return(contrasts_matrix)
}

# --- Orthogonalization Generator ---
create_structured_orthogonal_terms <- function(main_effects_matrix) {
  X_mat <- as.matrix(main_effects_matrix)
  N <- nrow(X_mat)
  m <- ncol(X_mat)
  col_names <- colnames(X_mat)
  X_center <- scale(X_mat, center = TRUE, scale = FALSE)
  
  quad_terms <- matrix(NA, nrow=N, ncol=m)
  for(i in 1:m) {
    quad_raw <- X_mat[,i]^2
    quad_terms[,i] <- lm.fit(cbind(1, X_center[,i]), quad_raw)$residuals
  }
  colnames(quad_terms) <- paste0(col_names, "_sq_orth")
  
  n_int <- choose(m, 2)
  if(n_int > 0) {
    int_terms <- matrix(NA, nrow=N, ncol=n_int)
    curr <- 1
    cnames <- c()
    for(i in 1:(m-1)) {
      for(j in (i+1):m) {
        int_raw <- X_mat[,i] * X_mat[,j]
        basis <- cbind(1, X_center[,i], X_center[,j], quad_terms[,i], quad_terms[,j])
        int_terms[,curr] <- lm.fit(basis, int_raw)$residuals
        cnames <- c(cnames, paste0(col_names[i], "x", col_names[j], "_orth"))
        curr <- curr + 1
      }
    }
    colnames(int_terms) <- cnames
    return(cbind(quad_terms, int_terms))
  } else {
    return(quad_terms)
  }
}

# --- Raw Complex Terms Generator ---
create_raw_complex_terms <- function(main_effects_matrix) {
  X_mat <- as.matrix(main_effects_matrix)
  N <- nrow(X_mat)
  m <- ncol(X_mat)
  col_names <- colnames(X_mat)
  
  quad_terms <- X_mat^2
  colnames(quad_terms) <- paste0(col_names, "_sq_RAW")
  
  n_int <- choose(m, 2)
  if(n_int > 0) {
    int_terms <- matrix(NA, nrow=N, ncol=n_int)
    cnames <- c()
    curr <- 1
    for(i in 1:(m-1)) {
      for(j in (i+1):m) {
        int_terms[,curr] <- X_mat[,i] * X_mat[,j]
        cnames <- c(cnames, paste0(col_names[i], "x", col_names[j], "_RAW"))
        curr <- curr + 1
      }
    }
    colnames(int_terms) <- cnames
    return(cbind(quad_terms, int_terms))
  } else {
    colnames(quad_terms) <- paste0(col_names, "_sq_RAW")
    return(quad_terms)
  }
}

# MSE Calculation
calc_mse_manual <- function(idx, data_df, X_full_matrix, truth) {
  if(length(idx) < 50) return(NA)
  
  mod <- tryCatch({
    suppressMessages(
      lmer(y ~ . - Group + (1 | Group), data = data_df[idx, ], 
           control = lmerControl(optimizer="bobyqa", calc.derivs=FALSE,
                                 check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))
    )
  }, error = function(e) { return(NULL) })
  
  if(is.null(mod)) return(NA)
  
  beta_hat <- fixef(mod)
  all_vars <- c("(Intercept)", colnames(X_full_matrix))
  beta_full <- setNames(rep(0, length(all_vars)), all_vars)
  beta_full[names(beta_hat)] <- beta_hat
  
  X_full_with_int <- cbind(1, X_full_matrix)
  pred_fixed <- as.vector(X_full_with_int %*% beta_full)
  
  ran_effs <- ranef(mod)$Group
  re_map <- setNames(ran_effs[,1], rownames(ran_effs))
  groups_in_full <- as.character(data_df$Group)
  pred_random <- re_map[groups_in_full]
  pred_random[is.na(pred_random)] <- 0
  
  preds <- pred_fixed + pred_random
  return(mean((preds - truth)^2))
}

# ==============================================================================
# 3. MAIN SIMULATION PROCESS
# ==============================================================================

# 1. Load Data
cat("Loading Cleaned Data...\n")
X_data <- read.csv("London_Airbnb.csv", stringsAsFactors = FALSE)

# 2. Define Variables
y <- log(X_data$price)
var_min_nights <- X_data$minimum_nights
var_reviews    <- X_data$number_of_reviews
var_host_count <- X_data$calculated_host_listings_count
var_pay        <- X_data$pay_2024
var_crime      <- X_data$crime_rate

# 3. Construct Main Effects Matrix
X_cont_raw <- cbind(var_min_nights, var_reviews, var_host_count)
X_cont_scaled <- scale(X_cont_raw)
colnames(X_cont_scaled) <- c("MinNights", "Reviews", "HostCount")

X_room_poly <- create_scaled_poly_contrasts(factor(X_data$room_type), "RoomType")

X_group_vars <- cbind(var_pay, var_crime)
X_group_vars_scaled <- scale(X_group_vars)
colnames(X_group_vars_scaled) <- c("Pay", "Crime")

X_MAIN_EFFECTS <- cbind(X_cont_scaled, X_room_poly, X_group_vars_scaled)

# 4. Construct Sampling & Modeling Matrices
cat("Generating Sampling Matrix...\n")
X_complex_ortho <- create_structured_orthogonal_terms(X_MAIN_EFFECTS)
X_SAMPLING_FULL <- cbind(X_MAIN_EFFECTS, X_complex_ortho)
X_within_sampling <- X_SAMPLING_FULL

cat("Generating Modeling Matrix...\n")
X_complex_raw <- create_raw_complex_terms(X_MAIN_EFFECTS)
X_MODELING_FULL <- cbind(X_MAIN_EFFECTS, X_complex_raw)

colnames(X_MODELING_FULL) <- paste0("V", 1:ncol(X_MODELING_FULL))
DF_Modeling <- as.data.frame(X_MODELING_FULL)
DF_Modeling$y <- y
DF_Modeling$Group <- factor(X_data$borough) 

X_MODELING_MATRIX <- as.matrix(X_MODELING_FULL)

# Group Sampling Matrix
Group_Info <- X_data %>% distinct(borough, pay_2024, crime_rate) %>% arrange(borough)
X_group_sampling <- scale(Group_Info[, c("pay_2024", "crime_rate")])

# Group Stats
C_i <- as.vector(X_data %>% group_by(borough) %>% count() %>% pull(n))
C_cum <- c(0, cumsum(C_i))
R <- length(C_i)

# 5. Simulation Loop
n_sub_target <- 1000 
nloop <- 50 
ratio <- 0.9

group_strats_base <- c("LEV", "Dopt", "IBOSS", "UNIF")
within_strats <- c("GOSS", "GIBOSS", "GLEV", "GUNIF")
col_names <- c("FullGroup", group_strats_base) 

cat("\n==================================================\n")
cat("Running Simulation for Optimal N:", n_sub_target, "\n")
cat("==================================================\n")

m <- ceiling(n_sub_target / R) 
results_matrix_list <- list()
for(ws in within_strats) {
  results_matrix_list[[ws]] <- matrix(NA, nrow = nloop, ncol = length(col_names))
  colnames(results_matrix_list[[ws]]) <- col_names
}

last_goss_indices <- NULL

for (k in 1:nloop) {
  if(k %% 10 == 0) cat("  Iter:", k, "/", nloop, "\n")
  set.seed(k * 2025 + n_sub_target) 
  
  # --- Phase 1: Full Group ---
  indices_full <- list(GOSS=c(), GIBOSS=c(), GLEV=c(), GUNIF=c())
  for (i in 1:R) {
    row_start <- C_cum[i] + 1; row_end <- C_cum[i+1]
    curr_n <- C_i[i]; curr_m <- min(m, curr_n)
    if(curr_m == 0) next
    
    X_sub <- X_within_sampling[row_start:row_end, , drop=FALSE]
    X_sub_scaled_goss <- apply(X_sub, 2, scalex)
    
    # Safe scaling for IBOSS/LEV to prevent NaN
    X_sub_std <- safe_scale(X_sub)
    
    sel <- OAJ2_cpp(X_sub_scaled_goss, curr_m, tPow=2); indices_full$GOSS <- c(indices_full$GOSS, row_start + sel - 1)
    sel <- iboss(X_sub_std, curr_m); indices_full$GIBOSS <- c(indices_full$GIBOSS, row_start + sel - 1)
    
    # Calculate prob with safe check
    prob <- lev(X_sub_std)
    if(length(prob) != curr_n || any(is.na(prob))) prob <- rep(1/curr_n, curr_n) # Fallback
    sel <- sample.int(curr_n, curr_m, replace=FALSE, prob=prob); indices_full$GLEV <- c(indices_full$GLEV, row_start + sel - 1)
    
    sel <- sample.int(curr_n, curr_m, replace=FALSE); indices_full$GUNIF <- c(indices_full$GUNIF, row_start + sel - 1)
  }
  
  if (k == nloop) last_goss_indices <- indices_full$GOSS
  
  for(ws in within_strats) {
    results_matrix_list[[ws]][k, "FullGroup"] <- calc_mse_manual(indices_full[[ws]], DF_Modeling, X_MODELING_MATRIX, y)
  }
  
  # --- Phase 2: Group Selection ---
  n_groups_sel <- floor(R * ratio)
  
  for (gs in group_strats_base) {
    selected_groups <- switch(gs,
                              LEV = { prob_g <- lev(scale(X_group_sampling)); sort(sample.int(R, n_groups_sel, replace=FALSE, prob=prob_g)) },
                              Dopt = { res <- optFederov(~., data = as.data.frame(X_group_sampling), nTrials = n_groups_sel, criterion = "D", evaluateI=FALSE); sort(res$rows) },
                              IBOSS = { sort(iboss(scale(X_group_sampling), n_groups_sel)) },
                              UNIF = { sort(sample.int(R, n_groups_sel, replace=FALSE)) }
    )
    
    indices_sub <- list(GOSS=c(), GIBOSS=c(), GLEV=c(), GUNIF=c())
    
    for (g_idx in selected_groups) {
      row_start <- C_cum[g_idx] + 1; row_end <- C_cum[g_idx+1]
      curr_n <- C_i[g_idx]; curr_m <- min(m, curr_n)
      if(curr_m == 0) next
      
      X_sub <- X_within_sampling[row_start:row_end, , drop=FALSE]
      X_sub_scaled_goss <- apply(X_sub, 2, scalex)
      X_sub_std <- safe_scale(X_sub)
      
      sel <- OAJ2_cpp(X_sub_scaled_goss, curr_m, tPow=2); indices_sub$GOSS <- c(indices_sub$GOSS, row_start + sel - 1)
      sel <- iboss(X_sub_std, curr_m); indices_sub$GIBOSS <- c(indices_sub$GIBOSS, row_start + sel - 1)
      
      prob <- lev(X_sub_std)
      if(length(prob) != curr_n || any(is.na(prob))) prob <- rep(1/curr_n, curr_n)
      sel <- sample.int(curr_n, curr_m, replace=FALSE, prob=prob); indices_sub$GLEV <- c(indices_sub$GLEV, row_start + sel - 1)
      
      sel <- sample.int(curr_n, curr_m, replace=FALSE); indices_sub$GUNIF <- c(indices_sub$GUNIF, row_start + sel - 1)
    }
    
    for(ws in within_strats) {
      results_matrix_list[[ws]][k, gs] <- calc_mse_manual(indices_sub[[ws]], DF_Modeling, X_MODELING_MATRIX, y)
    }
  } 
} 

# 6. Results
final_mse_table <- matrix(0, nrow=length(within_strats), ncol=length(col_names))
rownames(final_mse_table) <- within_strats
colnames(final_mse_table) <- col_names

for(ws in within_strats) {
  final_mse_table[ws, ] <- colMeans(results_matrix_list[[ws]], na.rm=TRUE)
}

cat("\n=== MSE Results (N=1000) ===\n")
print(final_mse_table)
write.csv(final_mse_table, "LondonAirbnb_MSE_N1000.csv")

# ==============================================================================
# 4. SIGNIFICANT FACTOR VISUALIZATION & EXPORT
# ==============================================================================
cat("\n=== Generating Visualization ===\n")

# 1. Refit Model using Last GOSS Sample
mod_goss <- lmer(y ~ . - Group + (1 | Group), data = DF_Modeling[last_goss_indices, ],
                 control = lmerControl(optimizer="bobyqa", 
                                       check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)))

# 2. Extract Coefficients
summ <- summary(mod_goss)
coefs <- summ$coefficients
coef_df <- data.frame(
  Term = rownames(coefs),
  Estimate = coefs[, "Estimate"], 
  tValue = coefs[, "t value"]
)

# 3. Map V-names back to Real Names
real_names_vec <- c(colnames(X_MAIN_EFFECTS), colnames(X_complex_raw))
v_names_vec    <- paste0("V", 1:length(real_names_vec))
name_map       <- setNames(real_names_vec, v_names_vec)
coef_df$RealName <- name_map[coef_df$Term]

# 4. Name Cleaning Function
clean_term_names <- function(raw_names) {
  new_names <- raw_names
  new_names <- gsub("_RAW", "", new_names)
  new_names <- gsub("RoomType\\.factor_column\\.L", "RoomType [Linear]", new_names)
  new_names <- gsub("RoomType\\.factor_column\\.Q", "RoomType [Quad]", new_names)
  new_names <- gsub("RoomType\\.factor_column\\.C", "RoomType [Cubic]", new_names)
  new_names <- gsub("RoomType\\.factor_column", "RoomType", new_names)
  new_names <- gsub("x", " × ", new_names)
  new_names <- gsub("_sq", "²", new_names)
  new_names <- gsub("MinNights", "Min. Stay", new_names)
  new_names <- gsub("HostCount", "Host Listings", new_names)
  new_names <- gsub("Reviews", "Review Count", new_names)
  new_names <- gsub("Pay", "Avg. Income", new_names)
  new_names <- gsub("Crime", "Crime Rate", new_names)
  return(new_names)
}

coef_df$CleanTerm <- clean_term_names(coef_df$RealName)

# 5. Filter Top 10 Significant Factors
top_10 <- coef_df %>% 
  filter(!is.na(CleanTerm)) %>% 
  arrange(desc(abs(tValue))) %>% 
  head(10)

# 6. Create Plot
p <- ggplot(top_10, aes(x = reorder(CleanTerm, abs(tValue)), y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  geom_segment(aes(x = CleanTerm, xend = CleanTerm, y = 0, yend = Estimate), 
               color = "grey60", size = 1.2) +
  geom_point(aes(color = abs(tValue)), size = 6) +
  coord_flip() +
  scale_color_gradient(low = "#5D9CEC", high = "#E04F5F", name = "|t-value|") +
  labs(
    title = "Top 10 Significant Predictors of Airbnb Price",
    subtitle = paste0("GOSS Subsampling (N=", n_sub_target, ")"),
    x = "", 
    y = "Coefficient Estimate (Log Price)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

# 7. Export High-Quality PNG
ggsave("Top10_Factors_GOSS.png", plot = p, width = 10, height = 7, dpi = 300)
cat("Plot saved as 'Top10_Factors_GOSS.png'. Analysis Complete.\n")
