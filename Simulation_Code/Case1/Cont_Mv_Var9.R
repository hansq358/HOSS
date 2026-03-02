# HOSS Simulation Case 1: Continuous MVN
# Symbols match the paper (N, R, C_i, N_G, N_I, sigma_a2, sigma_e2, B)
# --- 0. PRE-REQUISITES ---

rm(list=ls())
library(MASS)
library(caret)
library(AlgDesign)
Rcpp::sourceCpp('myoss.cpp')

# --- 1. HELPER FUNCTIONS ---

find.sigma=function(X, y, beta, C_i, R){
  
  eta.hat<-y-cbind(1, X)%*%beta
  
  b <- c()
  Ue <- 0
  for (i in 1:length(eta.hat)) {
    b[i] <- sum((eta.hat[i]*rep(1,length(eta.hat)) - eta.hat)^2)/2
    Ue <- Ue + b[i]
  }
  
  if(length(which(C_i==0))!=0){
    C_i <- C_i[-which(C_i==0)]
  }
  R.1 <- length(C_i)
  
  e1 <- c()
  for(j in 1:C_i[1]) {
    e1[j] <- sum((eta.hat[1:C_i[1]][j]*rep(1,C_i[1]) - eta.hat[1:C_i[1]])^2)/(2*C_i[1])
  }
  
  
  e <- c()
  e[1] <- sum(e1)
  
  for (i in 2:R.1) {
    e2 <- c()
    for (j in 1:C_i[i]) {
      e2[j] <- sum((eta.hat[(sum(C_i[1:(i-1)]) + j)]*rep(1,C_i[i]) - eta.hat[(sum(C_i[1:(i-1)]) + 1):sum(C_i[1:i])])^2)/(2*C_i[i])
    }
    e[i] <- sum(e2)
  }
  
  
  Ua <- sum(e)
  NF <- sum(C_i)
  NS <- sum(C_i^2)
  sigma_e2_hat <- Ua/(NF-R) 
  sigma_a2_hat <- Ue/(NF^2-NS) -  Ua*(NF^2 - NF)/((NF^2 - NS)*(NF-R)) 
  return(list(sigma_a2_hat, sigma_e2_hat))
}

find.beta=function(X, y, sigma_a2, sigma_e2, C_i, R, p){
  if(length(which(C_i==0))!=0){
    C_i <- C_i[-which(C_i==0)]
  }
  R.1 <- length(C_i)
  sn <- c(0, cumsum(C_i))
  
  XVX <- matrix(0, p+1, p+1)
  XVY <- matrix(0, p+1, 1)
  for (i in 1:R.1){
    gamma <- (C_i[i]*sigma_a2)/(sigma_e2+C_i[i]*sigma_a2)
    inv.V <- (1/sigma_e2)*(diag(1,(C_i[i]))-(gamma/(C_i[i]))*matrix(1,(C_i[i]),(C_i[i])))
    if(C_i[i] ==1){
      XVX <- XVX + t(cbind(1,t(X[(sn[i]+1):(sn[i+1]),])))%*%inv.V%*%cbind(1,t(X[(sn[i]+1):(sn[i+1]),]))
      XVY <- XVY + t(cbind(1,t(X[(sn[i]+1):(sn[i+1]),])))%*%inv.V%*%y[(sn[i]+1):(sn[i+1])]
    }else{
      XVX <- XVX + t(cbind(1,X[(sn[i]+1):(sn[i+1]),]))%*%inv.V%*%cbind(1,X[(sn[i]+1):(sn[i+1]),])
      XVY <- XVY + t(cbind(1,X[(sn[i]+1):(sn[i+1]),]))%*%inv.V%*%y[(sn[i]+1):(sn[i+1])]
    }
  }
  
  bt<-solve(XVX)%*%XVY
  
  return(bt)
}

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

Est_hat=function(X, y, beta, sigma_a2, sigma_e2, C_i, R, p){
  
  beta0 <- as.matrix(lm(y ~ X)$coefficients)
  sigma.hat0 <- find.sigma(X, y, beta0, C_i, R)
  beta1 <- find.beta(X, y, sigma.hat0[[1]],sigma.hat0[[2]], C_i, R, p)
  sigma.hat1 <- find.sigma(X, y, beta1, C_i, R)
  beta2 <- find.beta(X, y, sigma.hat1[[1]], sigma.hat1[[2]], C_i, R, p)
  bt.mse <- sum((beta2[-1]-beta)^2)
  
  return(list(bt.mse))
}

scalex=function(a){
  2*(a-min(a))/(max(a)-min(a))-1
}

standardize_matrix <- function(X) {
  X <- as.matrix(X)
  apply(X, 2, function(col) {
    if (sd(col) == 0) return(col)
    return((col - mean(col)) / sd(col))
  })
}

# Create Scaled Polynomial Contrasts for Categorical Variables
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

# Create Orthogonalized High-Order Terms (Gram-Schmidt)
create_structured_orthogonal_terms <- function(main_effects_matrix) {
  
  if (!is.matrix(main_effects_matrix)) {
    main_effects_matrix <- as.matrix(main_effects_matrix)
  }
  
  N <- nrow(main_effects_matrix)
  m <- ncol(main_effects_matrix)
  
  if (m == 0) {
    return(matrix(nrow = N, ncol = 0))
  }
  
  col_names <- colnames(main_effects_matrix)
  if (is.null(col_names)) {
    col_names <- paste0("V", 1:m)
  }
  
  X_main_centered <- scale(main_effects_matrix, center = TRUE, scale = FALSE)
  
  n_quad <- m
  quad_terms_ortho <- matrix(NA, nrow = N, ncol = n_quad)
  
  for (i in 1:m) {
    raw_squared_term <- main_effects_matrix[, i]^2
    linear_term <- X_main_centered[, i]
    fit_quad <- lm.fit(x = cbind(1, linear_term), y = raw_squared_term)
    quad_terms_ortho[, i] <- fit_quad$residuals
  }
  colnames(quad_terms_ortho) <- paste0(col_names, "_quad_ortho")
  
  n_int <- if (m > 1) choose(m, 2) else 0
  if (n_int == 0) {
    return(quad_terms_ortho)
  }
  
  interaction_terms_ortho <- matrix(NA, nrow = N, ncol = n_int)
  interaction_col_names <- character(n_int) 
  current_col_idx <- 1 
  
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      raw_interaction <- main_effects_matrix[, i] * main_effects_matrix[, j]
      
      basis_matrix <- cbind(
        1,
        X_main_centered[, i],
        X_main_centered[, j],
        quad_terms_ortho[, i],
        quad_terms_ortho[, j]
      )
      
      fit_int <- lm.fit(x = basis_matrix, y = raw_interaction)
      
      interaction_terms_ortho[, current_col_idx] <- fit_int$residuals
      interaction_col_names[current_col_idx] <- paste0(col_names[i], "x", col_names[j], "_int_ortho")
      current_col_idx <- current_col_idx + 1
    }
  }
  
  colnames(interaction_terms_ortho) <- interaction_col_names
  final_orthogonal_terms <- cbind(quad_terms_ortho, interaction_terms_ortho)
  
  return(final_orthogonal_terms)
}

# --- 2. MAIN SIMULATION FUNCTION (All Continuous Strategy) ---

mysimu = function(R, sigma_a2, sigma_e2, B, n_values) {
  
  # --- Step 2.1: Factor Structure (All Continuous) ---
  N_G <- 5
  N_I <- 5
  
  # Active Terms Config (Strong Signals)
  N_G_active_complex <- 3
  N_I_active_complex <- 5
  
  Sigma_G <- diag(0.5, N_G, N_G) + matrix(0.5, N_G, N_G)
  Sigma_I <- diag(0.5, N_I, N_I) + matrix(0.5, N_I, N_I)
  
  # --- Step 2.2: Dimensions ---
  n_main_cols <- N_G + N_I
  n_total_vars <- n_main_cols
  
  # Calculate number of complex terms
  n_quad <- n_total_vars
  n_int <- choose(n_total_vars, 2)
  n_complex_total <- n_quad + n_int
  total_fxx_cols <- n_main_cols + n_complex_total
  
  cat("--- Model Configuration (All Continuous Strategy) ---\n")
  cat("Group-level factors (N_G): ", N_G, " continuous.\n", sep="")
  cat("Individual-level factors (N_I): ", N_I, " continuous.\n", sep="")
  cat("Total main effect columns (signal=1.0):", n_main_cols, "\n")
  cat("Total strong signal complex terms (signal=0.5):", N_G_active_complex + N_I_active_complex, "\n")
  cat("Total weak signal complex terms (signal=0.05):", n_complex_total - (N_G_active_complex + N_I_active_complex), "\n")
  cat("Total model columns (Evaluation Space):", total_fxx_cols, "\n\n")
  
  # --- Setup ---
  lrs <- length(n_values)
  ratios_to_test <- c(0.9, 0.95)
  group_strats_base <- c("LEV", "Dopt", "IBOSS", "OSS", "Unif")
  expanded_group_strats <- c("FullGroup", as.vector(outer(group_strats_base, ratios_to_test, paste, sep = "_")))
  within_strats <- c("GOSS", "GUNIF", "GLEV", "GIBOSS")
  res_names <- apply(expand.grid(expanded_group_strats, within_strats), 1, function(x) paste(x[1], x[2], sep = "_"))
  results_list <- lapply(res_names, function(x) matrix(NA, nrow = 1, ncol = B * lrs))
  names(results_list) <- paste0(res_names, ".mat")
  
  itr <- 0
  for (n in n_values) {
    m <- ceiling(n / R)
    for (k in 1:B) {
      if (k %% 10 == 0) cat("n=", n, ", k=", k, "\n")
      itr <- itr + 1
      
      # >>> Seed A: HTC & Structure <<<
      set.seed(k * n * 123)
      
      # --- Structure Generation ---
      if (R %% 5 != 0) {
        stop("R must be a multiple of 5.")
      }
      C_i = sample(c(rep(300, R / 5), rep(400, R / 5), rep(500, R / 5),
                    rep(600, R / 5), rep(700, R / 5)))
      C_cum <- c(0, cumsum(C_i))
      a_i <- rnorm(R, mean = 0, sd = sqrt(sigma_a2))
      
      # --- Group-level covariates (raw) ---
      Z_group_raw <- mvrnorm(R, rep(0, N_G), Sigma_G)
      colnames(Z_group_raw) <- paste0("z_g", 1:N_G)
      
      # --- PATH 1: Matrix A (Group Sampling - Orthogonalized) ---
      # We use the helper to create orthogonal Quad and Int terms for SAMPLING only
      Z_group_orth <- create_structured_orthogonal_terms(Z_group_raw)
      X_group_select <- cbind(Z_group_raw, Z_group_orth)
      X_group_select_std <- standardize_matrix(X_group_select)
      
      # >>> Seed B: Individual-level covariates <<<
      set.seed(k * n * 123 + 99999)
      
      # --- Individual-level covariates (raw) ---
      X_ind <- matrix(NA, nrow = max(C_cum), ncol = N_I)
      for (i in 1:R) { 
        X_ind[((C_cum[i] + 1):C_cum[i + 1]), ] <- mvrnorm(C_i[i], rep(-1 + i / (R/2), N_I), Sigma_I) 
      }
      colnames(X_ind) <- paste0("x_i", 1:N_I)
      
      # --- PATH 1: Matrix B (Within Sampling - Orthogonalized) ---
      # We use the helper to create orthogonal Quad and Int terms for SAMPLING only
      X_ind_orth <- create_structured_orthogonal_terms(X_ind)
      X_within_select <- cbind(X_ind, X_ind_orth)
      X_within_select_std <- standardize_matrix(X_within_select)
      
      # --- PATH 2: X_full (evaluation, raw terms) ---
      # Expand group-level covariates to population level
      Z_group_expanded <- matrix(0, nrow = max(C_cum), ncol = N_G)
      for (i in 1:R) Z_group_expanded[((C_cum[i] + 1):C_cum[i + 1]), ] <- rep(Z_group_raw[i, ], each = C_i[i])
      colnames(Z_group_expanded) <- colnames(Z_group_raw)
      
      # Combine Main Effects
      X_main <- cbind(Z_group_expanded, X_ind)
      cont_names <- colnames(X_main)
      
      # Generate RAW Quadratic Terms (for Model Truth)
      X_quad_raw <- X_main^2
      colnames(X_quad_raw) <- paste0(cont_names, "_sq")
      
      # Generate RAW Interaction Terms (for Model Truth)
      X_int_list <- list()
      if (n_total_vars > 1) {
        for (i in 1:(n_total_vars - 1)) {
          for (j in (i + 1):n_total_vars) {
            term <- X_main[, i] * X_main[, j]
            X_int_list[[paste0(cont_names[i], "x", cont_names[j])]] <- term
          }
        }
      }
      X_int_raw <- do.call(cbind, X_int_list)
      
      # Final Evaluation Matrix
      X_full <- cbind(X_main, X_quad_raw, X_int_raw)
      
      # --- 5. Beta & Response ---
      beta <- rep(0.05, ncol(X_full))
      # 5.1: Main Effects = 1.0
      beta[1:n_main_cols] <- 1.0 
      
      # 5.2: Identify Strong Signals (Complex Terms)
      complex_col_names <- colnames(X_full)[(n_main_cols + 1):ncol(X_full)]
      
      # Filter for pure group-level complex terms
      htc_pure_complex <- grep("z_g", complex_col_names, value=TRUE)
      htc_pure_complex <- grep("x_i", htc_pure_complex, value=TRUE, invert=TRUE)
      
      # Filter for pure individual-level complex terms
      etc_pure_complex <- grep("x_i", complex_col_names, value=TRUE)
      etc_pure_complex <- grep("z_g", etc_pure_complex, value=TRUE, invert=TRUE)
      
      # Select specific terms to be Strong (0.5)
      active_htc_names <- htc_pure_complex[1:min(length(htc_pure_complex), N_G_active_complex)]
      active_etc_names <- etc_pure_complex[1:min(length(etc_pure_complex), N_I_active_complex)]
      
      idx_htc_active <- match(active_htc_names, colnames(X_full))
      idx_etc_active <- match(active_etc_names, colnames(X_full))
      
      beta[idx_htc_active] <- 0.5
      beta[idx_etc_active] <- 0.5
      
      epsilon <- rnorm(max(C_cum), mean = 0, sd = sqrt(sigma_e2))
      a_i_rep <- rep(a_i, C_i)
      y <- 1 + X_full %*% beta + a_i_rep + epsilon
      
      # --- 6. LOOPS ---
      
      # --- Full Group ---
      selected_groups_full <- 1:R
      all_indices_full <- list(
        GOSS = unlist(sapply(selected_groups_full, function(g_idx) {
          rows <- (C_cum[g_idx] + 1):C_cum[g_idx + 1]; size <- min(m, C_i[g_idx]); if (size == 0) return(NULL)
          rows[OAJ2_cpp(apply(X_within_select[rows, , drop = FALSE], 2, scalex), size, tPow = 2)]
        })),
        GUNIF = unlist(sapply(selected_groups_full, function(g_idx) {
          rows <- (C_cum[g_idx] + 1):C_cum[g_idx + 1]; size <- min(m, C_i[g_idx]); if (size == 0) return(NULL)
          rows[sample.int(C_i[g_idx], size, replace = FALSE)]
        })),
        GLEV = unlist(sapply(selected_groups_full, function(g_idx) {
          rows <- (C_cum[g_idx] + 1):C_cum[g_idx + 1]; size <- min(m, C_i[g_idx]); if (size == 0) return(NULL)
          rows[sample.int(C_i[g_idx], size, replace = FALSE, prob = lev(X_within_select_std[rows, , drop = FALSE]))]
        })),
        GIBOSS = unlist(sapply(selected_groups_full, function(g_idx) {
          rows <- (C_cum[g_idx] + 1):C_cum[g_idx + 1]; size <- min(m, C_i[g_idx]); if (size == 0) return(NULL)
          rows[iboss(X_within_select_std[rows, , drop = FALSE], size)]
        }))
      )
      
      for (w_strat_name in names(all_indices_full)) {
        res_mat_name <- paste0("FullGroup_", w_strat_name, ".mat")
        indices <- all_indices_full[[w_strat_name]]
        if (!is.null(indices) && length(indices) > total_fxx_cols + 5) {
          C_i_sub <- sapply(1:R, function(i) sum(indices > C_cum[i] & indices <= C_cum[i + 1]))
          est_result <- Est_hat(X = X_full[indices, , drop = FALSE], y = y[indices], 
                                beta = beta, sigma_a2, sigma_e2, C_i_sub, R, p = total_fxx_cols)
          results_list[[res_mat_name]][, itr] <- est_result[[1]]
        }
      }
      
      # --- Sub-Group Selection ---
      for (g_strat_base in group_strats_base) {
        for (ratio in ratios_to_test) {
          n_groups_to_select <- floor(R * ratio)
          
          selected_groups <- switch(g_strat_base,
                                    LEV   = sort(sample.int(R, size = n_groups_to_select, replace = FALSE, prob = lev(X_group_select_std))),
                                    Dopt  = sort(optFederov(~., data = as.data.frame(X_group_select_std), nTrials = n_groups_to_select, criterion = "D")$rows),
                                    IBOSS = sort(iboss(X_group_select_std, n_groups_to_select)),
                                    OSS   = sort(OAJ2_cpp(apply(X_group_select, 2, scalex), n_groups_to_select, tPow = 2)),
                                    Unif  = sort(sample.int(R, size = n_groups_to_select, replace = FALSE))
          )
          
          if (length(selected_groups) == 0) next
          
          all_indices <- list(
            GOSS   = unlist(sapply(selected_groups, function(g_idx) {
              rows <- (C_cum[g_idx] + 1):C_cum[g_idx + 1]; size <- min(m, C_i[g_idx]); if (size == 0) return(NULL)
              rows[OAJ2_cpp(apply(X_within_select[rows, , drop = FALSE], 2, scalex), size, tPow = 2)]
            })),
            GUNIF  = unlist(sapply(selected_groups, function(g_idx) {
              rows <- (C_cum[g_idx] + 1):C_cum[g_idx + 1]; size <- min(m, C_i[g_idx]); if (size == 0) return(NULL)
              rows[sample.int(C_i[g_idx], size, replace = FALSE)]
            })),
            GLEV   = unlist(sapply(selected_groups, function(g_idx) {
              rows <- (C_cum[g_idx] + 1):C_cum[g_idx + 1]; size <- min(m, C_i[g_idx]); if (size == 0) return(NULL)
              rows[sample.int(C_i[g_idx], size, replace = FALSE, prob = lev(X_within_select_std[rows, , drop = FALSE]))]
            })),
            GIBOSS = unlist(sapply(selected_groups, function(g_idx) {
              rows <- (C_cum[g_idx] + 1):C_cum[g_idx + 1]; size <- min(m, C_i[g_idx]); if (size == 0) return(NULL)
              rows[iboss(X_within_select_std[rows, , drop = FALSE], size)]
            }))
          )
          
          current_g_strat_name <- paste(g_strat_base, ratio, sep = "_")
          for (w_strat_name in names(all_indices)) {
            res_mat_name <- paste0(current_g_strat_name, "_", w_strat_name, ".mat")
            indices <- all_indices[[w_strat_name]]
            if (!is.null(indices) && length(indices) > total_fxx_cols + 5) {
              C_i_sub <- sapply(1:R, function(i) sum(indices > C_cum[i] & indices <= C_cum[i + 1]))
              est_result <- Est_hat(X = X_full[indices, , drop = FALSE], y = y[indices], 
                                    beta = beta, sigma_a2, sigma_e2, C_i_sub, R, p = total_fxx_cols)
              results_list[[res_mat_name]][, itr] <- est_result[[1]]
            }
          }
        }
      }
    }
  }
  
  mse_list <- list()
  for (res_name in res_names) {
    res_mat <- results_list[[paste0(res_name, ".mat")]]
    mse_vec <- sapply(1:lrs, function(i) { loc <- ((i - 1) * B + 1):(i * B); mean(res_mat[, loc], na.rm = TRUE) })
    mse_list[[res_name]] <- mse_vec
  }
  rec1 <- do.call(cbind, mse_list); rownames(rec1) <- n_values; return(list(rec1))
}

# --- 3. CALL AND EXPORT ---

simulation_results <- mysimu(R = 100, 
                             sigma_a2 = 9, 
                             sigma_e2 = 9, 
                             B = 300, 
                             n_values = c(500, 1000, 1500, 2000))

# Print Results
print(simulation_results[[1]])

# Export Results to CSV
write.csv(simulation_results[[1]], file = "simulation_results_Case1_var9.csv")
