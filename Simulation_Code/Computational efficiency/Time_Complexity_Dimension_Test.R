rm(list=ls())
library(MASS)
library(AlgDesign)
library(Rcpp)
library(RcppArmadillo)

sourceCpp('myoss.cpp') 
sourceCpp('ortho_helper.cpp')

scalex = function(a) {
  2 * (a - min(a)) / (max(a) - min(a)) - 1
}

standardize_matrix <- function(X) {
  X <- as.matrix(X)
  apply(X, 2, function(col) {
    if (sd(col) == 0) return(col)
    (col - mean(col)) / sd(col)
  })
}

create_structured_orthogonal_terms_optimized <- function(main_effects_matrix) {
  if (!is.matrix(main_effects_matrix)) {
    main_effects_matrix <- as.matrix(main_effects_matrix)
  }
  return(create_structured_orthogonal_terms_cpp(main_effects_matrix))
}

make_unbalanced_group_sizes = function(N_total, R, sizes) {
  if (R %% length(sizes) != 0) stop("R must be a multiple of length(sizes).")
  base <- rep(sizes, each = R / length(sizes))
  Ni <- sample(base, size = R, replace = FALSE)
  if (sum(Ni) != N_total) stop("Sum of group sizes must equal N_total.")
  Ni
}

make_full_raw_features = function(X_main) {
  X_main <- as.matrix(X_main)
  n_i <- nrow(X_main)
  p <- ncol(X_main)
  X_quad <- X_main^2
  n_int <- p * (p - 1) / 2
  if (n_int == 0) return(cbind(X_main, X_quad))
  X_int <- matrix(0, n_i, n_int)
  idx <- 1
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      X_int[, idx] <- X_main[, i] * X_main[, j]
      idx <- idx + 1
    }
  }
  cbind(X_main, X_quad, X_int)
}

run_dimension_simulation = function(
  N_total = 1000000,
  R = 100,
  n = 1000,
  nloop = 100,
  group_ratio = 0.9,
  p_group = 5,
  p_within_list = c(5, 10, 15, 20),
  rho_group = 0.5,
  rho_within = 0.5
) {
  results <- data.frame()
  n_groups_to_select <- floor(R * group_ratio)
  m_sub <- ceiling(n / n_groups_to_select)
  stopifnot(n_groups_to_select >= 2)
  stopifnot(m_sub >= 1)
  cat(sprintf("N=%d, R=%d, n=%d, K=%d, m=%d, loops=%d\n", N_total, R, n, n_groups_to_select, m_sub, nloop))
  
  for (p_within in p_within_list) {
    cat(sprintf("p_group=%d, p_within=%d\n", p_group, p_within))

    time_hoss_total <- numeric(nloop)
    time_hoss_group_prep <- numeric(nloop)
    time_hoss_group_select <- numeric(nloop)
    time_hoss_within_prep <- numeric(nloop)
    time_hoss_within_sample <- numeric(nloop)

    time_goss_within_sample <- numeric(nloop)
    time_goss_total <- numeric(nloop)

    sigma_group <- matrix(rho_group, p_group, p_group)
    diag(sigma_group) <- 1
    sigma_within <- matrix(rho_within, p_within, p_within)
    diag(sigma_within) <- 1

    for (k in 1:nloop) {
      set.seed(1000 * k + 10 * p_within)

      Ni <- make_unbalanced_group_sizes(N_total, R, sizes = c(8000, 9000, 10000, 11000, 12000))
      SC <- c(0, cumsum(Ni))
      group_id <- rep(1:R, times = Ni)

      Z_group <- mvrnorm(R, rep(0, p_group), sigma_group)
      X_within_main <- mvrnorm(N_total, rep(0, p_within), sigma_within)
      mu <- matrix(rep(seq(-1, 1, length.out = R), p_within), nrow = R, ncol = p_within)
      X_within_main <- X_within_main + mu[group_id, , drop = FALSE]

      t0_hoss <- Sys.time()
      t_gp0 <- Sys.time()
      Z_group_ortho <- create_structured_orthogonal_terms_optimized(Z_group)
      Z_group_tilde <- cbind(Z_group, Z_group_ortho)
      Z_group_std <- standardize_matrix(Z_group_tilde)
      time_hoss_group_prep[k] <- as.numeric(Sys.time() - t_gp0, units = "secs")

      t_gs0 <- Sys.time()
      selected_groups <- sort(optFederov(~., data = as.data.frame(Z_group_std), nTrials = n_groups_to_select, criterion = "D")$rows)
      time_hoss_group_select[k] <- as.numeric(Sys.time() - t_gs0, units = "secs")

      unlist(lapply(selected_groups, function(g_idx) {
        rows <- (SC[g_idx] + 1):SC[g_idx + 1]
        size <- min(m_sub, Ni[g_idx])
        if (size <= 0) return(NULL)

        t_wpi <- Sys.time()
        Xg_main <- X_within_main[rows, , drop = FALSE]
        Xg_ortho <- create_structured_orthogonal_terms_optimized(Xg_main)
        Xg_tilde <- cbind(Xg_main, Xg_ortho)
        time_hoss_within_prep[k] <<- time_hoss_within_prep[k] + as.numeric(Sys.time() - t_wpi, units = "secs")

        t_wsi <- Sys.time()
        group_data <- apply(Xg_tilde, 2, scalex)
        out <- rows[OAJ2_cpp(group_data, size, tPow = 2)]
        time_hoss_within_sample[k] <<- time_hoss_within_sample[k] + as.numeric(Sys.time() - t_wsi, units = "secs")
        out
      }))

      time_hoss_total[k] <- as.numeric(Sys.time() - t0_hoss, units = "secs")

      unlist(lapply(selected_groups, function(g_idx) {
        rows <- (SC[g_idx] + 1):SC[g_idx + 1]
        size <- min(m_sub, Ni[g_idx])
        if (size <= 0) return(NULL)

        Xg_full_raw <- make_full_raw_features(X_within_main[rows, , drop = FALSE])

        t_samp <- Sys.time()
        group_data <- apply(Xg_full_raw, 2, scalex)
        out <- rows[OAJ2_cpp(group_data, size, tPow = 2)]
        time_goss_within_sample[k] <<- time_goss_within_sample[k] + as.numeric(Sys.time() - t_samp, units = "secs")
        out
      }))
      time_goss_total[k] <- time_hoss_group_prep[k] + time_hoss_group_select[k] + time_goss_within_sample[k]
    }
    
    # HOSS: Avg_Time = group_prep + group_select + within_prep + within_sample
    hoss_avg_time <- mean(time_hoss_group_prep + time_hoss_group_select + 
                          time_hoss_within_prep + time_hoss_within_sample, na.rm = TRUE)
    
    results <- rbind(results, data.frame(
      N_total = N_total,
      R = R,
      n = n,
      K = n_groups_to_select,
      m = m_sub,
      p_group = p_group,
      p_within = p_within,
      rho_group = rho_group,
      rho_within = rho_within,
      Method = "HOSS",
      Avg_Time = hoss_avg_time,
      Avg_GroupPreprocess = mean(time_hoss_group_prep, na.rm = TRUE),
      Avg_GroupSelection = mean(time_hoss_group_select, na.rm = TRUE),
      Avg_WithinPreprocess = mean(time_hoss_within_prep, na.rm = TRUE),
      Avg_WithinSampling = mean(time_hoss_within_sample, na.rm = TRUE)
    ))

    # GOSS: Avg_Time = group_prep + group_select + within_sample 
    goss_avg_time <- mean(time_hoss_group_prep + time_hoss_group_select + 
                          time_goss_within_sample, na.rm = TRUE)
    
    results <- rbind(results, data.frame(
      N_total = N_total,
      R = R,
      n = n,
      K = n_groups_to_select,
      m = m_sub,
      p_group = p_group,
      p_within = p_within,
      rho_group = rho_group,
      rho_within = rho_within,
      Method = "GOSS",
      Avg_Time = goss_avg_time,
      Avg_GroupPreprocess = mean(time_hoss_group_prep, na.rm = TRUE),
      Avg_GroupSelection = mean(time_hoss_group_select, na.rm = TRUE),
      Avg_WithinPreprocess = NA,
      Avg_WithinSampling = mean(time_goss_within_sample, na.rm = TRUE)
    ))

    cat(sprintf("HOSS Avg_Time: %.2fs (group prep %.2fs + group select %.2fs + within prep %.2fs + within sample %.2fs)\n",
                hoss_avg_time,
                mean(time_hoss_group_prep, na.rm = TRUE),
                mean(time_hoss_group_select, na.rm = TRUE),
                mean(time_hoss_within_prep, na.rm = TRUE),
                mean(time_hoss_within_sample, na.rm = TRUE)))
    cat(sprintf("GOSS Avg_Time: %.2fs (group prep %.2fs + group select %.2fs + within sample %.2fs)\n",
                goss_avg_time,
                mean(time_hoss_group_prep, na.rm = TRUE),
                mean(time_hoss_group_select, na.rm = TRUE),
                mean(time_goss_within_sample, na.rm = TRUE)))
  }
  
  return(results)
}

final_res <- run_dimension_simulation()

print(final_res)

write.csv(final_res, "Time_Complexity_HOSS_Full_Results.csv", row.names = FALSE)

