argv0 <- commandArgs(trailingOnly = FALSE)
file_arg <- argv0[grep("^--file=", argv0)]
if (length(file_arg) > 0) {
  root <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  root <- getwd()
}
argv <- commandArgs(trailingOnly = TRUE)
mode <- NULL
if (length(argv) > 0 && !grepl("^--", argv[1])) {
  mode <- argv[1]
  args <- argv[-1]
} else {
  args <- argv
}

get_arg <- function(name, default = NULL) {
  eq <- paste0(name, "=")
  hit <- grep(paste0("^", eq), args)
  if (length(hit) > 0) return(sub(eq, "", args[hit[1]], fixed = TRUE))
  pos <- match(name, args)
  if (!is.na(pos) && pos < length(args)) return(args[pos + 1])
  default
}

csv_chr <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])
csv_num <- function(x) as.numeric(csv_chr(x))
csv_int <- function(x) as.integer(csv_num(x))

suppressPackageStartupMessages({
  library(MASS)
  library(AlgDesign)
  library(Rcpp)
  library(RcppArmadillo)
  library(parallel)
})

Rcpp::sourceCpp(file.path(root, "myoss.cpp"))

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
  return(list(sigma_a2_hat,sigma_e2_hat))
}

find.beta=function(X, y, sigma_a2, sigma_e2, C_i, R, p){
  if(length(which(C_i==0))!=0){
    C_i <- C_i[-which(C_i==0)]
  }
  R.1 <- length(C_i)
  C_cum <- c(0, cumsum(C_i))
  XVX <- matrix(0, p+1, p+1)
  XVY <- matrix(0, p+1, 1)
  for (i in 1:R.1){
    gamma <- (C_i[i]*sigma_a2)/(sigma_e2+C_i[i]*sigma_a2)
    inv.V <- (1/sigma_e2)*(diag(1,(C_i[i]))-(gamma/(C_i[i]))*matrix(1,(C_i[i]),(C_i[i])))
    if(C_i[i] ==1){
      XVX <- XVX + t(cbind(1,t(X[(C_cum[i]+1):(C_cum[i+1]),])))%*%inv.V%*%cbind(1,t(X[(C_cum[i]+1):(C_cum[i+1]),]))
      XVY <- XVY + t(cbind(1,t(X[(C_cum[i]+1):(C_cum[i+1]),])))%*%inv.V%*%y[(C_cum[i]+1):(C_cum[i+1])]
    }else{
      XVX <- XVX + t(cbind(1,X[(C_cum[i]+1):(C_cum[i+1]),]))%*%inv.V%*%cbind(1,X[(C_cum[i]+1):(C_cum[i+1]),])
      XVY <- XVY + t(cbind(1,X[(C_cum[i]+1):(C_cum[i+1]),]))%*%inv.V%*%y[(C_cum[i]+1):(C_cum[i+1])]
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

create_scaled_contrasts <- function(factor_column, factor_name, contrast_type) {
  if (contrast_type == "poly") return(create_scaled_poly_contrasts(factor_column, factor_name))
  if (!is.factor(factor_column)) factor_column <- factor(as.character(factor_column))
  num_levels <- nlevels(factor_column)
  if (num_levels <= 1) return(matrix(nrow = length(factor_column), ncol = 0))
  if (contrast_type == "helmert") {
    cm <- contr.helmert(num_levels)
  } else if (contrast_type == "sum") {
    cm <- qr.Q(qr(contr.sum(num_levels)))
  } else {
    mm <- model.matrix(~ factor_column)[, -1, drop = FALSE]
    colnames(mm) <- paste(factor_name, colnames(mm), sep = ".")
    return(mm)
  }
  contrasts_matrix <- model.matrix(~ factor_column, contrasts.arg = list(factor_column = cm))[, -1, drop = FALSE]
  colnames(contrasts_matrix) <- paste(factor_name, colnames(contrasts_matrix), sep = ".")
  contrasts_matrix
}

make_group_sizes <- function(R, N_total) {
  base <- sample(c(rep(300, R / 5), rep(400, R / 5), rep(500, R / 5), rep(600, R / 5), rep(700, R / 5)))
  if (is.na(N_total) || sum(base) == N_total) return(base)
  sizes <- pmax(2, floor(base * N_total / sum(base)))
  diff <- N_total - sum(sizes)
  pos <- 1
  while (diff != 0) {
    if (diff > 0) {
      sizes[pos] <- sizes[pos] + 1
      diff <- diff - 1
    } else if (sizes[pos] > 2) {
      sizes[pos] <- sizes[pos] - 1
      diff <- diff + 1
    }
    pos <- ifelse(pos == length(sizes), 1, pos + 1)
  }
  sample(sizes)
}

case_dims <- function(case) {
  if (case %in% c("case1", "case2")) list(mixed = FALSE, uniform = case == "case2", N_G_CONT = 5, N_I_CONT = 5) else list(mixed = TRUE, uniform = case == "case4", N_G_CONT = 3, N_I_CONT = 3)
}

make_raw_terms <- function(X_main, names_main) {
  X_quad_raw <- X_main^2
  colnames(X_quad_raw) <- paste0(names_main, "_sq")
  X_int_list <- list()
  if (ncol(X_main) > 1) {
    for (i in 1:(ncol(X_main) - 1)) {
      for (j in (i + 1):ncol(X_main)) {
        term <- X_main[, i] * X_main[, j]
        X_int_list[[paste0(names_main[i], "x", names_main[j])]] <- term
      }
    }
  }
  X_int_raw <- do.call(cbind, X_int_list)
  cbind(X_quad_raw, X_int_raw)
}

make_beta <- function(X_full, n_main_cols) {
  beta <- rep(0.05, ncol(X_full))
  beta[1:n_main_cols] <- 1.0
  complex_col_names <- colnames(X_full)[(n_main_cols + 1):ncol(X_full)]
  group_pool <- grep("z_g", complex_col_names, value=TRUE)
  group_pool <- grep("x_i", group_pool, value=TRUE, invert=TRUE)
  ind_pool <- grep("x_i", complex_col_names, value=TRUE)
  ind_pool <- grep("z_g", ind_pool, value=TRUE, invert=TRUE)
  idx_group_active <- match(group_pool[1:min(length(group_pool), 3)], colnames(X_full))
  idx_ind_active <- match(ind_pool[1:min(length(ind_pool), 5)], colnames(X_full))
  beta[idx_group_active] <- 0.5
  beta[idx_ind_active] <- 0.5
  beta
}

permute_factor_df <- function(df) {
  out <- df
  for (nm in names(df)) {
    vals <- sort(unique(df[[nm]]))
    map <- sample(vals)
    names(map) <- vals
    out[[nm]] <- unname(map[as.character(df[[nm]])])
  }
  out
}

make_population <- function(case, seed, R, N_total, sigma_a2, sigma_e2, contrast_type = "poly", label_permute = FALSE, label_seed = NULL, truth_model = "full") {
  set.seed(seed)
  cfg <- case_dims(case)
  C_i <- make_group_sizes(R, N_total)
  C_cum <- c(0, cumsum(C_i))
  N <- max(C_cum)
  a_i <- rnorm(R, mean = 0, sd = sqrt(sigma_a2))
  if (cfg$uniform) {
    Z_group_cont_raw <- matrix(runif(R * cfg$N_G_CONT, -1, 1), R, cfg$N_G_CONT)
  } else {
    Sigma_G <- diag(0.5, cfg$N_G_CONT, cfg$N_G_CONT) + matrix(0.5, cfg$N_G_CONT, cfg$N_G_CONT)
    Z_group_cont_raw <- mvrnorm(R, rep(0, cfg$N_G_CONT), Sigma_G)
  }
  colnames(Z_group_cont_raw) <- paste0("z_g", 1:cfg$N_G_CONT)
  if (cfg$uniform) {
    X_ind_cont <- matrix(NA, nrow = N, ncol = cfg$N_I_CONT)
    for (i in 1:R) X_ind_cont[((C_cum[i] + 1):C_cum[i + 1]), ] <- matrix(runif(C_i[i]*cfg$N_I_CONT, -2+i/(R/2), 0+i/(R/2)), C_i[i], cfg$N_I_CONT)
  } else {
    Sigma_I <- diag(0.5, cfg$N_I_CONT, cfg$N_I_CONT) + matrix(0.5, cfg$N_I_CONT, cfg$N_I_CONT)
    X_ind_cont <- matrix(NA, nrow = N, ncol = cfg$N_I_CONT)
    for (i in 1:R) X_ind_cont[((C_cum[i] + 1):C_cum[i + 1]), ] <- mvrnorm(C_i[i], rep(-1 + i / (R/2), cfg$N_I_CONT), Sigma_I)
  }
  colnames(X_ind_cont) <- paste0("x_i", 1:cfg$N_I_CONT)
  if (cfg$mixed) {
    G_CAT_FACTORS <- list(g_cat1 = 1:3, g_cat2 = 1:2)
    I_CAT_FACTORS <- list(i_cat1 = 1:2, i_cat2 = 1:3)
    Z_group_cat_raw <- as.data.frame(sapply(G_CAT_FACTORS, function(l) sample(l, R, replace = TRUE), simplify = FALSE))
    X_ind_cat_raw <- as.data.frame(sapply(I_CAT_FACTORS, function(l) sample(l, N, replace = TRUE), simplify = FALSE))
    if (label_permute) {
      if (!is.null(label_seed)) set.seed(label_seed)
      Z_group_cat_raw <- permute_factor_df(Z_group_cat_raw)
      X_ind_cat_raw <- permute_factor_df(X_ind_cat_raw)
    }
    list_group_contrasts <- mapply(function(col, nm) create_scaled_contrasts(col, nm, contrast_type), Z_group_cat_raw, names(Z_group_cat_raw), SIMPLIFY = FALSE)
    Z_group_cat_contrast <- do.call(cbind, list_group_contrasts)
    list_ind_contrasts <- mapply(function(col, nm) create_scaled_contrasts(col, nm, contrast_type), X_ind_cat_raw, names(X_ind_cat_raw), SIMPLIFY = FALSE)
    X_ind_cat_contrasts <- do.call(cbind, list_ind_contrasts)
  } else {
    Z_group_cat_raw <- data.frame()
    X_ind_cat_raw <- data.frame()
    Z_group_cat_contrast <- matrix(nrow = R, ncol = 0)
    X_ind_cat_contrasts <- matrix(nrow = N, ncol = 0)
  }
  Z_group_orth <- create_structured_orthogonal_terms(Z_group_cont_raw)
  X_group_select <- cbind(Z_group_cont_raw, Z_group_cat_contrast, Z_group_orth)
  X_group_select_std <- standardize_matrix(X_group_select)
  X_ind_orth <- create_structured_orthogonal_terms(X_ind_cont)
  X_within_select <- cbind(X_ind_cont, X_ind_cat_contrasts, X_ind_orth)
  X_within_select_std <- standardize_matrix(X_within_select)
  Z_group_cont_expanded <- matrix(0, nrow = N, ncol = ncol(Z_group_cont_raw))
  for (i in 1:R) Z_group_cont_expanded[((C_cum[i] + 1):C_cum[i + 1]), ] <- rep(Z_group_cont_raw[i, ], each = C_i[i])
  colnames(Z_group_cont_expanded) <- colnames(Z_group_cont_raw)
  if (cfg$mixed) {
    list_group_dummies <- list()
    for(cname in names(Z_group_cat_raw)) {
      fac <- factor(Z_group_cat_raw[[cname]])
      mm <- model.matrix(~ fac)[, -1, drop=FALSE]
      colnames(mm) <- paste0("g_", cname, "_lev", levels(fac)[-1])
      list_group_dummies[[cname]] <- mm
    }
    Z_group_cat_dummy <- do.call(cbind, list_group_dummies)
    Z_group_cat_dummy_expanded <- matrix(0, nrow = N, ncol = ncol(Z_group_cat_dummy))
    for (i in 1:R) Z_group_cat_dummy_expanded[((C_cum[i] + 1):C_cum[i + 1]), ] <- rep(Z_group_cat_dummy[i, ], each = C_i[i])
    colnames(Z_group_cat_dummy_expanded) <- colnames(Z_group_cat_dummy)
    list_ind_dummies <- list()
    for(cname in names(X_ind_cat_raw)) {
      fac <- factor(X_ind_cat_raw[[cname]])
      mm <- model.matrix(~ fac)[, -1, drop=FALSE]
      colnames(mm) <- paste0("i_", cname, "_lev", levels(fac)[-1])
      list_ind_dummies[[cname]] <- mm
    }
    X_ind_cat_dummies <- do.call(cbind, list_ind_dummies)
    X_main <- cbind(Z_group_cont_expanded, Z_group_cat_dummy_expanded, X_ind_cont, X_ind_cat_dummies)
  } else {
    X_main <- cbind(Z_group_cont_expanded, X_ind_cont)
  }
  cont_names <- colnames(cbind(Z_group_cont_expanded, X_ind_cont))
  X_cont_all <- cbind(Z_group_cont_expanded, X_ind_cont)
  X_complex <- make_raw_terms(X_cont_all, cont_names)
  X_full <- cbind(X_main, X_complex)
  n_main_cols <- ncol(X_main)
  beta_full <- make_beta(X_full, n_main_cols)
  X_quad_no <- X_cont_all^2
  colnames(X_quad_no) <- paste0(colnames(X_cont_all), "_sq")
  X_no_interactions <- cbind(X_main, X_quad_no)
  beta_no <- beta_full[match(colnames(X_no_interactions), colnames(X_full))]
  beta_main <- beta_full[1:n_main_cols]
  beta_truth <- beta_full
  X_truth <- X_full
  if (truth_model == "main") {
    X_truth <- X_main
    beta_truth <- beta_main
  }
  epsilon <- rnorm(N, mean = 0, sd = sqrt(sigma_e2))
  y <- 1 + X_truth %*% beta_truth + rep(a_i, C_i) + epsilon
  list(case = case, R = R, C_i = C_i, C_cum = C_cum, N = N, X_group_select = X_group_select, X_group_select_std = X_group_select_std, X_within_select = X_within_select, X_within_select_std = X_within_select_std, X_full = X_full, X_main = X_main, X_no_interactions = X_no_interactions, beta_full = beta_full, beta_main = beta_main, beta_no_interactions = beta_no, y = as.numeric(y), sigma_a2 = sigma_a2, sigma_e2 = sigma_e2, Z_group_cat_raw = Z_group_cat_raw, X_ind_cat_raw = X_ind_cat_raw)
}

method_grid <- function(ratios, within = c("GOSS", "GUNIF", "GLEV", "GIBOSS")) {
  group_strats_base <- c("LEV", "Dopt", "IBOSS", "OSS", "Unif")
  c(paste0("FullGroup_", within), as.vector(sapply(group_strats_base, function(g) as.vector(sapply(ratios, function(r) paste0(g, "_", r, "_", within))))))
}

methods_from_arg <- function(value, ratios) {
  if (value == "workflow") return(method_grid(ratios))
  if (value == "core") return(c("FullGroup_GOSS", "Dopt_0.9_GOSS", "OSS_0.9_GOSS", "Unif_0.9_GUNIF"))
  csv_chr(value)
}

select_groups <- function(pop, g_strat_base, ratio) {
  R <- pop$R
  n_groups_to_select <- floor(R * ratio)
  switch(g_strat_base,
         LEV   = sort(sample.int(R, size = n_groups_to_select, replace = FALSE, prob = lev(pop$X_group_select_std))),
         Dopt  = sort(optFederov(~., data = as.data.frame(pop$X_group_select_std), nTrials = n_groups_to_select, criterion = "D")$rows),
         IBOSS = sort(iboss(pop$X_group_select_std, n_groups_to_select)),
         OSS   = sort(OAJ2_cpp(apply(pop$X_group_select, 2, scalex), n_groups_to_select, tPow = 2)),
         Unif  = sort(sample.int(R, size = n_groups_to_select, replace = FALSE)))
}

select_within <- function(pop, selected_groups, n, w_strat_name) {
  m <- ceiling(n / pop$R)
  unlist(sapply(selected_groups, function(g_idx) {
    rows <- (pop$C_cum[g_idx] + 1):pop$C_cum[g_idx + 1]
    size <- min(m, pop$C_i[g_idx])
    if (size == 0) return(NULL)
    if (w_strat_name == "GOSS") return(rows[OAJ2_cpp(apply(pop$X_within_select[rows, , drop = FALSE], 2, scalex), size, tPow = 2)])
    if (w_strat_name == "GUNIF") return(rows[sample.int(pop$C_i[g_idx], size, replace = FALSE)])
    if (w_strat_name == "GLEV") return(rows[sample.int(pop$C_i[g_idx], size, replace = FALSE, prob = lev(pop$X_within_select_std[rows, , drop = FALSE]))])
    rows[iboss(pop$X_within_select_std[rows, , drop = FALSE], size)]
  }))
}

select_method <- function(pop, method, n) {
  parts <- strsplit(method, "_", fixed = TRUE)[[1]]
  if (parts[1] == "FullGroup") {
    selected_groups <- 1:pop$R
    w <- parts[2]
  } else {
    selected_groups <- select_groups(pop, parts[1], as.numeric(parts[2]))
    w <- parts[3]
  }
  select_within(pop, selected_groups, n, w)
}

fit_matrix <- function(pop, fit_model) {
  if (fit_model == "main") return(list(X = pop$X_main, beta = pop$beta_main))
  if (fit_model == "no_interactions") return(list(X = pop$X_no_interactions, beta = pop$beta_no_interactions))
  list(X = pop$X_full, beta = pop$beta_full)
}

Est_hat_rev <- function(pop, indices, fit_model) {
  fm <- fit_matrix(pop, fit_model)
  X_all <- fm$X
  beta <- fm$beta
  X <- X_all[indices, , drop = FALSE]
  y <- pop$y[indices]
  C_i_sub <- sapply(1:pop$R, function(i) sum(indices > pop$C_cum[i] & indices <= pop$C_cum[i + 1]))
  p <- ncol(X)
  if (length(indices) <= p + 5) return(data.frame(selected_n = length(indices), beta_mse = NA, mspe = NA, sigma_a2_hat = NA, sigma_e2_hat = NA, sigma_a2_bias = NA, sigma_e2_bias = NA, condition_number = NA))
  beta0 <- as.matrix(lm(y ~ X)$coefficients)
  sigma.hat0 <- find.sigma(X, y, beta0, C_i_sub, pop$R)
  beta1 <- find.beta(X, y, sigma.hat0[[1]],sigma.hat0[[2]], C_i_sub, pop$R, p)
  sigma.hat1 <- find.sigma(X, y, beta1, C_i_sub, pop$R)
  beta2 <- find.beta(X, y, sigma.hat1[[1]], sigma.hat1[[2]], C_i_sub, pop$R, p)
  sigma.hat2 <- find.sigma(X, y, beta2, C_i_sub, pop$R)
  bt.mse <- sum((beta2[-1]-beta)^2)
  resid <- as.numeric(y - cbind(1, X) %*% beta2)
  group_eff <- rep(0, pop$R)
  start <- 1
  for (g in 1:pop$R) {
    c <- C_i_sub[g]
    if (c > 0) {
      rg <- resid[start:(start + c - 1)]
      gamma <- (c * sigma.hat2[[1]]) / (sigma.hat2[[2]] + c * sigma.hat2[[1]])
      group_eff[g] <- gamma * mean(rg)
      start <- start + c
    }
  }
  group_id <- rep(1:pop$R, pop$C_i)
  yhat <- as.numeric(cbind(1, X_all) %*% beta2 + group_eff[group_id])
  data.frame(selected_n = length(indices), beta_mse = bt.mse, mspe = mean((yhat - pop$y)^2), sigma_a2_hat = sigma.hat2[[1]], sigma_e2_hat = sigma.hat2[[2]], sigma_a2_bias = sigma.hat2[[1]] - pop$sigma_a2, sigma_e2_bias = sigma.hat2[[2]] - pop$sigma_e2, condition_number = kappa(cbind(1, X), exact = FALSE))
}

task_prediction <- function(task) {
  pop <- make_population(task$case, task$seed, task$R, task$N_total, task$sigma_a2, task$sigma_e2, truth_model = task$truth_model)
  rows <- lapply(task$methods, function(method) {
    set.seed(task$seed + nchar(method) + task$n)
    indices <- select_method(pop, method, task$n)
    out <- Est_hat_rev(pop, indices, task$fit_model)
    cbind(data.frame(case = task$case, scenario = task$scenario, rep = task$rep, seed = task$seed, n_total = pop$N, n_sub = task$n, method = method, truth_model = task$truth_model, fit_model = task$fit_model), out)
  })
  do.call(rbind, rows)
}

collect <- function(tasks, worker, jobs) {
  jobs <- min(max(1, jobs), length(tasks))
  if (.Platform$OS.type == "windows") jobs <- 1
  out <- mclapply(tasks, worker, mc.cores = jobs)
  do.call(rbind, out)
}

summarize_df <- function(df, keys) {
  metrics <- intersect(c("selected_n", "beta_mse", "mspe", "sigma_a2_hat", "sigma_e2_hat", "sigma_a2_bias", "sigma_e2_bias", "condition_number", "jaccard_vs_poly", "jaccard_vs_unpermuted", "d_efficiency", "rho_max", "rho_mean", "time_total", "time_group_prep", "time_group_select", "time_within_prep", "time_within_sample"), names(df))
  groups <- split(df, interaction(df[keys], drop = TRUE), drop = TRUE)
  rows <- lapply(groups, function(g) {
    head <- g[1, keys, drop = FALSE]
    head$replications <- nrow(g)
    for (m in metrics) head[[m]] <- mean(g[[m]], na.rm = TRUE)
    head
  })
  do.call(rbind, rows)
}

write_outputs <- function(out_dir, prefix, raw, keys) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(raw, file.path(out_dir, paste0(prefix, "_raw.csv")), row.names = FALSE)
  write.csv(summarize_df(raw, keys), file.path(out_dir, paste0(prefix, "_summary.csv")), row.names = FALSE)
}

run_prediction <- function(case_override = NULL) {
  case <- if (is.null(case_override)) get_arg("--case", "case1") else case_override
  B <- as.integer(get_arg("--B", "300"))
  R <- as.integer(get_arg("--R", "100"))
  N_total <- as.integer(get_arg("--N-total", "50000"))
  n_values <- csv_int(get_arg("--n-values", "500,1000,1500,2000"))
  sigma_a2 <- as.numeric(get_arg("--sigma-a2", "0.5"))
  sigma_e2 <- as.numeric(get_arg("--sigma-e2", "9"))
  ratios <- csv_num(get_arg("--ratios", "0.9,0.95"))
  methods <- methods_from_arg(get_arg("--methods", "workflow"), ratios)
  jobs <- as.integer(get_arg("--jobs", Sys.getenv("NJOBS", "128")))
  seed0 <- as.integer(get_arg("--seed", "20260603"))
  out_dir <- get_arg("--output-dir", file.path(root, "results"))
  tasks <- list()
  idx <- 1
  for (n in n_values) for (rep in 1:B) {
    tasks[[idx]] <- list(case = case, rep = rep, seed = seed0 + rep * 1009 + n * 17, R = R, N_total = N_total, n = n, sigma_a2 = sigma_a2, sigma_e2 = sigma_e2, methods = methods, truth_model = "full", fit_model = "full", scenario = "default")
    idx <- idx + 1
  }
  raw <- collect(tasks, task_prediction, jobs)
  write_outputs(out_dir, paste0(case, "_mspe_variance"), raw, c("case", "n_sub", "method"))
}

run_prediction_case <- function(case_name) {
  run_prediction(case_name)
}

run_n_scan <- function(cases_override = NULL) {
  cases <- if (is.null(cases_override)) csv_chr(get_arg("--cases", "case1,case2,case3,case4")) else cases_override
  B <- as.integer(get_arg("--B", "300"))
  R <- as.integer(get_arg("--R", "100"))
  N_values <- csv_int(get_arg("--N-values", "50000,100000,200000,500000"))
  n_sub <- as.integer(get_arg("--n-sub", "1000"))
  sigma_a2 <- as.numeric(get_arg("--sigma-a2", "0.5"))
  sigma_e2 <- as.numeric(get_arg("--sigma-e2", "9"))
  methods <- methods_from_arg(get_arg("--methods", "core"), c(0.9, 0.95))
  jobs <- as.integer(get_arg("--jobs", Sys.getenv("NJOBS", "128")))
  seed0 <- as.integer(get_arg("--seed", "20260603"))
  out_dir <- get_arg("--output-dir", file.path(root, "results"))
  tasks <- list()
  idx <- 1
  for (case in cases) for (N_total in N_values) for (rep in 1:B) {
    tasks[[idx]] <- list(case = case, rep = rep, seed = seed0 + rep * 1009 + N_total + match(case, cases) * 100000, R = R, N_total = N_total, n = n_sub, sigma_a2 = sigma_a2, sigma_e2 = sigma_e2, methods = methods, truth_model = "full", fit_model = "full", scenario = "default")
    idx <- idx + 1
  }
  raw <- collect(tasks, task_prediction, jobs)
  write_outputs(out_dir, "mse_mspe_vs_N", raw, c("case", "n_total", "n_sub", "method"))
}

task_encoding <- function(task) {
  base_indices <- NULL
  rows <- list()
  pos <- 1
  for (enc in task$encodings) {
    pop <- make_population(task$case, task$seed, task$R, task$N_total, task$sigma_a2, task$sigma_e2, contrast_type = enc)
    set.seed(task$seed + nchar(enc))
    indices <- select_method(pop, task$method, task$n)
    if (enc == "poly") base_indices <- indices
    out <- Est_hat_rev(pop, indices, "full")
    rows[[pos]] <- cbind(data.frame(case = task$case, experiment = "encoding", rep = task$rep, seed = task$seed, n_total = pop$N, n_sub = task$n, method = task$method, encoding = enc, jaccard_vs_poly = NA, jaccard_vs_unpermuted = NA), out)
    attr(rows[[pos]], "indices") <- list(indices)
    pos <- pos + 1
  }
  for (i in seq_along(rows)) {
    indices <- attr(rows[[i]], "indices")[[1]]
    rows[[i]]$jaccard_vs_poly <- length(intersect(indices, base_indices)) / length(union(indices, base_indices))
    attr(rows[[i]], "indices") <- NULL
  }
  if (task$label_permutations > 0) {
    pop0 <- make_population(task$case, task$seed, task$R, task$N_total, task$sigma_a2, task$sigma_e2)
    set.seed(task$seed + 991)
    idx0 <- select_method(pop0, task$method, task$n)
    for (perm in 1:task$label_permutations) {
      popp <- make_population(task$case, task$seed, task$R, task$N_total, task$sigma_a2, task$sigma_e2, label_permute = TRUE, label_seed = task$seed + perm * 31)
      set.seed(task$seed + perm * 997)
      idxp <- select_method(popp, task$method, task$n)
      rows[[pos]] <- data.frame(case = task$case, experiment = "label_permutation", rep = task$rep, seed = task$seed, n_total = popp$N, n_sub = task$n, method = task$method, encoding = "poly", jaccard_vs_poly = NA, jaccard_vs_unpermuted = length(intersect(idxp, idx0)) / length(union(idxp, idx0)), selected_n = length(idxp), beta_mse = NA, mspe = NA, sigma_a2_hat = NA, sigma_e2_hat = NA, sigma_a2_bias = NA, sigma_e2_bias = NA, condition_number = NA)
      pos <- pos + 1
    }
  }
  do.call(rbind, rows)
}

run_encoding <- function() {
  case <- get_arg("--case", "case3")
  B <- as.integer(get_arg("--B", "300"))
  R <- as.integer(get_arg("--R", "100"))
  N_total <- as.integer(get_arg("--N-total", "50000"))
  n <- as.integer(get_arg("--n-sub", "1000"))
  sigma_a2 <- as.numeric(get_arg("--sigma-a2", "0.5"))
  sigma_e2 <- as.numeric(get_arg("--sigma-e2", "9"))
  method <- get_arg("--method", "Dopt_0.9_GOSS")
  encodings <- csv_chr(get_arg("--encodings", "poly,helmert,sum,dummy"))
  label_permutations <- as.integer(get_arg("--label-permutations", "100"))
  jobs <- as.integer(get_arg("--jobs", Sys.getenv("NJOBS", "128")))
  seed0 <- as.integer(get_arg("--seed", "20260603"))
  out_dir <- get_arg("--output-dir", file.path(root, "results"))
  tasks <- lapply(1:B, function(rep) list(case = case, rep = rep, seed = seed0 + rep * 1009, R = R, N_total = N_total, n = n, sigma_a2 = sigma_a2, sigma_e2 = sigma_e2, method = method, encodings = encodings, label_permutations = label_permutations))
  raw <- collect(tasks, task_encoding, jobs)
  write_outputs(out_dir, paste0(case, "_encoding_sensitivity"), raw, c("case", "experiment", "encoding", "method"))
}

run_misspecification <- function() {
  case <- get_arg("--case", "case4")
  B <- as.integer(get_arg("--B", "300"))
  R <- as.integer(get_arg("--R", "100"))
  N_total <- as.integer(get_arg("--N-total", "50000"))
  n <- as.integer(get_arg("--n-sub", "1000"))
  sigma_a2 <- as.numeric(get_arg("--sigma-a2", "0.5"))
  sigma_e2 <- as.numeric(get_arg("--sigma-e2", "9"))
  methods <- methods_from_arg(get_arg("--methods", "core"), c(0.9, 0.95))
  jobs <- as.integer(get_arg("--jobs", Sys.getenv("NJOBS", "128")))
  seed0 <- as.integer(get_arg("--seed", "20260603"))
  out_dir <- get_arg("--output-dir", file.path(root, "results"))
  scenarios <- data.frame(scenario = c("correct", "over_spec", "under_spec", "no_interactions"), truth_model = c("full", "main", "full", "full"), fit_model = c("full", "full", "main", "no_interactions"), stringsAsFactors = FALSE)
  tasks <- list()
  idx <- 1
  for (i in 1:nrow(scenarios)) for (rep in 1:B) {
    tasks[[idx]] <- list(case = case, rep = rep, seed = seed0 + rep * 1009 + i * 100000, R = R, N_total = N_total, n = n, sigma_a2 = sigma_a2, sigma_e2 = sigma_e2, methods = methods, truth_model = scenarios$truth_model[i], fit_model = scenarios$fit_model[i], scenario = scenarios$scenario[i])
    idx <- idx + 1
  }
  raw <- collect(tasks, task_prediction, jobs)
  write_outputs(out_dir, paste0(case, "_misspecification"), raw, c("case", "scenario", "method"))
}

logdet_info <- function(X) {
  M <- t(cbind(1, standardize_matrix(X))) %*% cbind(1, standardize_matrix(X))
  as.numeric(determinant(M + diag(1e-8, nrow(M)), logarithm = TRUE)$modulus)
}

task_sequential <- function(task) {
  pop <- make_population(task$case, task$seed, task$R, task$N_total, task$sigma_a2, task$sigma_e2)
  idx_seq <- select_method(pop, task$method, task$n)
  group_expanded <- do.call(rbind, lapply(1:pop$R, function(i) matrix(rep(pop$X_group_select[i, ], pop$C_i[i]), nrow = pop$C_i[i], byrow = TRUE)))
  joint_X <- cbind(group_expanded, pop$X_within_select)
  idx_joint <- sort(optFederov(~., data = as.data.frame(standardize_matrix(joint_X)), nTrials = length(idx_seq), criterion = "D")$rows)
  log_seq <- logdet_info(joint_X[idx_seq, , drop = FALSE])
  log_joint <- logdet_info(joint_X[idx_joint, , drop = FALSE])
  rho <- cancor(standardize_matrix(group_expanded[idx_seq, , drop = FALSE]), standardize_matrix(pop$X_within_select[idx_seq, , drop = FALSE]))$cor
  data.frame(case = task$case, experiment = "sequential_vs_joint", rep = task$rep, seed = task$seed, n_total = pop$N, n_sub = task$n, method = task$method, selected_n = length(idx_seq), joint_selected_n = length(idx_joint), logdet_sequential = log_seq, logdet_joint = log_joint, d_efficiency = exp((log_seq - log_joint) / ncol(cbind(1, joint_X))), rho_max = max(rho), rho_mean = mean(rho))
}

run_sequential <- function() {
  case <- get_arg("--case", "case3")
  B <- as.integer(get_arg("--B", "300"))
  R <- as.integer(get_arg("--R", "20"))
  N_total <- as.integer(get_arg("--N-total", "1000"))
  n <- as.integer(get_arg("--n-sub", "120"))
  sigma_a2 <- as.numeric(get_arg("--sigma-a2", "0.5"))
  sigma_e2 <- as.numeric(get_arg("--sigma-e2", "9"))
  method <- get_arg("--method", "Dopt_0.9_GOSS")
  jobs <- as.integer(get_arg("--jobs", Sys.getenv("NJOBS", "128")))
  seed0 <- as.integer(get_arg("--seed", "20260603"))
  out_dir <- get_arg("--output-dir", file.path(root, "results"))
  tasks <- lapply(1:B, function(rep) list(case = case, rep = rep, seed = seed0 + rep * 1009, R = R, N_total = N_total, n = n, sigma_a2 = sigma_a2, sigma_e2 = sigma_e2, method = method))
  raw <- collect(tasks, task_sequential, jobs)
  write_outputs(out_dir, paste0(case, "_sequential_vs_joint"), raw, c("case", "method", "n_total", "n_sub"))
}

make_full_raw_features <- function(X_main) {
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

task_dimension <- function(task) {
  set.seed(task$seed)
  n_groups_to_select <- floor(task$R * task$group_ratio)
  m_sub <- ceiling(task$n / n_groups_to_select)
  C_i <- make_group_sizes(task$R, task$N_total)
  C_cum <- c(0, cumsum(C_i))
  Sigma_G <- matrix(task$rho_group, task$p_group, task$p_group)
  diag(Sigma_G) <- 1
  Sigma_I <- matrix(task$rho_within, task$p_within, task$p_within)
  diag(Sigma_I) <- 1
  Z_group <- mvrnorm(task$R, rep(0, task$p_group), Sigma_G)
  tg0 <- Sys.time()
  Z_group_ortho <- create_structured_orthogonal_terms(Z_group)
  Z_group_tilde <- cbind(Z_group, Z_group_ortho)
  Z_group_std <- standardize_matrix(Z_group_tilde)
  tg1 <- Sys.time()
  selected_groups <- sort(optFederov(~., data = as.data.frame(Z_group_std), nTrials = n_groups_to_select, criterion = "D")$rows)
  tg2 <- Sys.time()
  X_within_main <- matrix(NA, nrow = sum(C_i), ncol = task$p_within)
  for (i in 1:task$R) X_within_main[(C_cum[i] + 1):C_cum[i + 1], ] <- mvrnorm(C_i[i], rep(-1 + i / (task$R/2), task$p_within), Sigma_I)
  time_within_prep <- 0
  time_within_sample <- 0
  time_goss_prep <- 0
  time_goss_sample <- 0
  for (g_idx in selected_groups) {
    rows <- (C_cum[g_idx] + 1):C_cum[g_idx + 1]
    size <- min(m_sub, C_i[g_idx])
    tw0 <- Sys.time()
    Xg_ortho <- create_structured_orthogonal_terms(X_within_main[rows, , drop = FALSE])
    Xg_tilde <- cbind(X_within_main[rows, , drop = FALSE], Xg_ortho)
    tw1 <- Sys.time()
    OAJ2_cpp(apply(Xg_tilde, 2, scalex), size, tPow = 2)
    tw2 <- Sys.time()
    tgoss0 <- Sys.time()
    Xg_full_raw <- make_full_raw_features(X_within_main[rows, , drop = FALSE])
    tgoss1 <- Sys.time()
    OAJ2_cpp(apply(Xg_full_raw, 2, scalex), size, tPow = 2)
    tgoss2 <- Sys.time()
    time_within_prep <- time_within_prep + as.numeric(tw1 - tw0, units = "secs")
    time_within_sample <- time_within_sample + as.numeric(tw2 - tw1, units = "secs")
    time_goss_prep <- time_goss_prep + as.numeric(tgoss1 - tgoss0, units = "secs")
    time_goss_sample <- time_goss_sample + as.numeric(tgoss2 - tgoss1, units = "secs")
  }
  group_prep <- as.numeric(tg1 - tg0, units = "secs")
  group_select <- as.numeric(tg2 - tg1, units = "secs")
  # Algorithm timings exclude simulated data generation. Both methods share the same group selection.
  hoss_total <- group_prep + group_select + time_within_prep + time_within_sample
  goss_total <- group_prep + group_select + time_goss_prep + time_goss_sample
  rbind(data.frame(experiment = "dimension_runtime", method = "HOSS", rep = task$rep, seed = task$seed, n_total = sum(C_i), n_sub = task$n, p_group = task$p_group, p_within = task$p_within, time_total = hoss_total, time_group_prep = group_prep, time_group_select = group_select, time_within_prep = time_within_prep, time_within_sample = time_within_sample), data.frame(experiment = "dimension_runtime", method = "GOSS", rep = task$rep, seed = task$seed, n_total = sum(C_i), n_sub = task$n, p_group = task$p_group, p_within = task$p_within, time_total = goss_total, time_group_prep = group_prep, time_group_select = group_select, time_within_prep = time_goss_prep, time_within_sample = time_goss_sample))
}

run_dimension <- function() {
  nloop <- as.integer(get_arg("--nloop", "100"))
  N_total <- as.integer(get_arg("--N-total", "1000000"))
  R <- as.integer(get_arg("--R", "100"))
  n <- as.integer(get_arg("--n-sub", "1000"))
  group_ratio <- as.numeric(get_arg("--group-ratio", "0.9"))
  p_group <- as.integer(get_arg("--p-group", "5"))
  p_within_list <- csv_int(get_arg("--p-within-list", "5,10,15,20,30,40,50"))
  rho_group <- as.numeric(get_arg("--rho-group", "0.5"))
  rho_within <- as.numeric(get_arg("--rho-within", "0.5"))
  jobs <- as.integer(get_arg("--jobs", Sys.getenv("NJOBS", "128")))
  seed0 <- as.integer(get_arg("--seed", "20260603"))
  out_dir <- get_arg("--output-dir", file.path(root, "results"))
  tasks <- list()
  idx <- 1
  for (p_within in p_within_list) for (rep in 1:nloop) {
    tasks[[idx]] <- list(rep = rep, seed = seed0 + rep * 1009 + p_within * 100000, N_total = N_total, R = R, n = n, group_ratio = group_ratio, p_group = p_group, p_within = p_within, rho_group = rho_group, rho_within = rho_within)
    idx <- idx + 1
  }
  raw <- collect(tasks, task_dimension, jobs)
  write_outputs(out_dir, "dimension_runtime", raw, c("experiment", "method", "p_within"))
}

if (!is.null(mode)) {
  if (mode == "prediction") run_prediction()
  if (mode == "n_scan") run_n_scan()
  if (mode == "encoding") run_encoding()
  if (mode == "misspecification") run_misspecification()
  if (mode == "sequential") run_sequential()
  if (mode == "dimension") run_dimension()
}
