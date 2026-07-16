Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

argv0 <- commandArgs(trailingOnly = FALSE)
file_arg <- argv0[grep("^--file=", argv0)]
if (length(file_arg) != 1L) stop("Run this file with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))

source(file.path(script_dir, "common.R"))
Rcpp::sourceCpp(file.path(script_dir, "weighted_joint_core.cpp"))

as_flag <- function(x) tolower(trimws(x)) %in% c("1", "true", "yes", "y")
analysis_version <- "2026-07-16-v1"

logdet_spd <- function(M) {
  R <- tryCatch(chol((M + t(M)) / 2), error = function(e) NULL)
  if (is.null(R)) stop("Information matrix is not positive definite.")
  2 * sum(log(diag(R)))
}

block_diagnostics <- function(M, p_group) {
  p_total <- nrow(M)
  if (p_group < 1L || p_group >= p_total) stop("Invalid group-block dimension.")
  idx_group <- seq_len(p_group)
  idx_within <- (p_group + 1L):p_total
  M_group <- M[idx_group, idx_group, drop = FALSE]
  M_within <- M[idx_within, idx_within, drop = FALSE]
  M_cross <- M[idx_group, idx_within, drop = FALSE]

  R_group <- chol((M_group + t(M_group)) / 2)
  R_within <- chol((M_within + t(M_within)) / 2)
  left_whitened <- forwardsolve(t(R_group), M_cross)
  whitened <- t(forwardsolve(t(R_within), t(left_whitened)))
  rho <- svd(whitened, nu = 0, nv = 0)$d
  rho <- pmin(1, pmax(0, rho))

  coupling_from_rho <- exp(sum(log(pmax(1e-15, 1 - rho^2))) / p_total)
  coupling_from_det <- exp(
    (logdet_spd(M) - logdet_spd(M_group) - logdet_spd(M_within)) / p_total
  )
  list(
    rho_max = max(rho),
    rho_mean = mean(rho),
    coupling_efficiency = coupling_from_rho,
    coupling_efficiency_check = coupling_from_det,
    coupling_check_error = abs(coupling_from_rho - coupling_from_det)
  )
}

expand_group_features <- function(pop) {
  do.call(rbind, lapply(seq_len(pop$R), function(i) {
    matrix(
      rep(pop$X_group_select[i, ], pop$C_i[i]),
      nrow = pop$C_i[i],
      byrow = TRUE
    )
  }))
}

choose_start <- function(candidates, X, group_id, lambda) {
  scores <- vapply(candidates, function(candidate) {
    weighted_logdet_cpp(X, group_id, as.integer(candidate$indices), lambda)$logdet
  }, numeric(1))
  best <- which.max(scores)
  list(
    indices = sort(as.integer(candidates[[best]]$indices)),
    source = candidates[[best]]$source,
    logdet = scores[best]
  )
}

legacy_efficiency <- function(joint_X, idx_sequential, idx_joint_unweighted) {
  p_total <- ncol(cbind(1, joint_X))
  log_seq <- logdet_info(joint_X[idx_sequential, , drop = FALSE])
  log_joint <- logdet_info(joint_X[idx_joint_unweighted, , drop = FALSE])
  exp((log_seq - log_joint) / p_total)
}

run_replication <- function(task) {
  set.seed(task$seed)
  pop <- make_population(
    task$case, task$seed, task$R, task$N_total,
    task$sigma_a2, task$sigma_e2
  )
  set.seed(task$seed + nchar(task$method) + task$n)
  idx_sequential <- sort(select_method(pop, task$method, task$n))

  group_expanded <- expand_group_features(pop)
  joint_X_raw <- cbind(group_expanded, pop$X_within_select)

  joint_X_standardized <- standardize_matrix(joint_X_raw)
  X_design <- cbind(`(Intercept)` = 1, joint_X_standardized)
  group_id <- rep(seq_len(pop$R), pop$C_i)
  p_group <- 1L + ncol(group_expanded)
  p_total <- ncol(X_design)

  set.seed(task$seed + 700001L)
  idx_unweighted <- sort(optFederov(
    ~.,
    data = as.data.frame(joint_X_standardized),
    nTrials = length(idx_sequential),
    criterion = "D"
  )$rows)
  legacy_eff <- legacy_efficiency(joint_X_raw, idx_sequential, idx_unweighted)

  rows <- vector("list", length(task$lambda_grid))
  previous_joint <- NULL
  for (j in seq_along(task$lambda_grid)) {
    lambda <- task$lambda_grid[j]
    candidates <- list(
      list(source = "unweighted_optFederov", indices = idx_unweighted),
      list(source = "sequential_HOSS", indices = idx_sequential)
    )
    if (!is.null(previous_joint)) {
      candidates[[length(candidates) + 1L]] <- list(
        source = "previous_lambda_warm_start",
        indices = previous_joint
      )
    }
    start <- choose_start(candidates, X_design, group_id, lambda)

    started <- proc.time()[[3]]
    exchange <- weighted_exchange_cpp(
      X_design,
      as.integer(group_id),
      as.integer(start$indices),
      lambda,
      task$max_exchanges,
      task$tolerance
    )
    elapsed <- proc.time()[[3]] - started
    idx_joint <- sort(as.integer(exchange$selected))
    previous_joint <- idx_joint

    M_sequential <- weighted_info_cpp(
      X_design, as.integer(group_id), as.integer(idx_sequential), lambda
    )
    M_joint <- weighted_info_cpp(
      X_design, as.integer(group_id), as.integer(idx_joint), lambda
    )
    logdet_sequential <- logdet_spd(M_sequential)
    logdet_joint <- logdet_spd(M_joint)
    d_efficiency <- exp((logdet_sequential - logdet_joint) / p_total)
    if (d_efficiency > 1 + 1e-7) {
      stop(sprintf(
        "Joint exchange failed to dominate the sequential start (rep=%d, lambda=%.8g).",
        task$rep, lambda
      ))
    }
    d_efficiency <- min(1, d_efficiency)

    block <- block_diagnostics(M_sequential, p_group)
    if (block$coupling_check_error > 1e-7) {
      stop(sprintf(
        "Canonical-correlation determinant check failed (rep=%d, lambda=%.8g).",
        task$rep, lambda
      ))
    }

    seq_counts <- tabulate(group_id[idx_sequential], nbins = pop$R)
    joint_counts <- tabulate(group_id[idx_joint], nbins = pop$R)
    rows[[j]] <- data.frame(
      analysis_version = analysis_version,
      case = task$case,
      experiment = "weighted_sequential_vs_joint",
      rep = task$rep,
      seed = task$seed,
      n_total = pop$N,
      groups_total = pop$R,
      n_requested_before_group_ratio = task$n,
      n_selected = length(idx_sequential),
      method = task$method,
      sigma_a2_data = task$sigma_a2,
      sigma_e2_data = task$sigma_e2,
      lambda = lambda,
      lambda_multiplier = if (task$true_lambda == 0) NA_real_ else lambda / task$true_lambda,
      assumed_icc = lambda / (1 + lambda),
      p_total = p_total,
      p_group_block = p_group,
      p_within_block = p_total - p_group,
      sequential_groups_selected = sum(seq_counts > 0),
      joint_groups_selected = sum(joint_counts > 0),
      joint_group_n_min_positive = min(joint_counts[joint_counts > 0]),
      joint_group_n_max = max(joint_counts),
      logdet_sequential = logdet_sequential,
      logdet_joint = logdet_joint,
      d_efficiency = d_efficiency,
      rho_max = block$rho_max,
      rho_mean = block$rho_mean,
      coupling_efficiency = block$coupling_efficiency,
      coupling_efficiency_check = block$coupling_efficiency_check,
      coupling_check_error = block$coupling_check_error,
      legacy_unweighted_efficiency = legacy_eff,
      exchange_start = start$source,
      optimizer_max_exchanges = task$max_exchanges,
      optimizer_tolerance = task$tolerance,
      exchange_count = as.integer(exchange$exchanges),
      exchange_converged = isTRUE(exchange$converged),
      exchange_numerical_stop = isTRUE(exchange$numerical_stop),
      exchange_last_gain = as.numeric(exchange$last_gain),
      exchange_jitter = as.numeric(exchange$jitter),
      exchange_elapsed_seconds = elapsed,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarize_results <- function(raw, true_lambda) {
  split_rows <- split(raw, raw$lambda)
  summary_rows <- lapply(split_rows, function(g) {
    q <- quantile(g$d_efficiency, c(0.05, 0.50, 0.95), na.rm = TRUE)
    data.frame(
      case = g$case[1],
      method = g$method[1],
      lambda = g$lambda[1],
      lambda_multiplier = if (true_lambda == 0) NA_real_ else g$lambda[1] / true_lambda,
      assumed_icc = g$assumed_icc[1],
      replications = nrow(g),
      d_efficiency_mean = mean(g$d_efficiency),
      d_efficiency_sd = sd(g$d_efficiency),
      d_efficiency_se = sd(g$d_efficiency) / sqrt(nrow(g)),
      d_efficiency_q05 = unname(q[1]),
      d_efficiency_median = unname(q[2]),
      d_efficiency_q95 = unname(q[3]),
      d_efficiency_min = min(g$d_efficiency),
      d_efficiency_max = max(g$d_efficiency),
      rho_max_mean = mean(g$rho_max),
      rho_mean_mean = mean(g$rho_mean),
      coupling_efficiency_mean = mean(g$coupling_efficiency),
      joint_groups_selected_mean = mean(g$joint_groups_selected),
      exchange_convergence_rate = mean(g$exchange_converged),
      exchange_numerical_stop_rate = mean(g$exchange_numerical_stop),
      exchange_count_mean = mean(g$exchange_count),
      exchange_elapsed_seconds_mean = mean(g$exchange_elapsed_seconds),
      legacy_unweighted_efficiency_mean = mean(g$legacy_unweighted_efficiency),
      legacy_unweighted_efficiency_gt1_rate = mean(g$legacy_unweighted_efficiency > 1),
      max_coupling_check_error = max(g$coupling_check_error),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, summary_rows)
  out[order(out$lambda), , drop = FALSE]
}

case <- get_arg("--case", "case3")
B <- as.integer(get_arg("--B", "300"))
R_groups <- as.integer(get_arg("--R", "20"))
N_total <- as.integer(get_arg("--N-total", "1000"))
n_sub <- as.integer(get_arg("--n-sub", "120"))
sigma_a2 <- as.numeric(get_arg("--sigma-a2", "0.5"))
sigma_e2 <- as.numeric(get_arg("--sigma-e2", "9"))
method <- get_arg("--method", "Dopt_0.9_GOSS")
jobs <- as.integer(get_arg("--jobs", Sys.getenv("NJOBS", "1")))
seed0 <- as.integer(get_arg("--seed", "20260603"))
max_exchanges <- as.integer(get_arg("--max-exchanges", "40"))
tolerance <- as.numeric(get_arg("--tolerance", "1e-9"))
resume <- as_flag(get_arg("--resume", "true"))
out_dir <- normalizePath(
  get_arg("--output-dir", file.path(script_dir, "results")),
  mustWork = FALSE
)

if (!is.finite(sigma_a2) || sigma_a2 < 0 || !is.finite(sigma_e2) || sigma_e2 <= 0) {
  stop("sigma-a2 must be non-negative and sigma-e2 must be positive.")
}
if (B < 1L || jobs < 1L || max_exchanges < 1L) stop("B, jobs, and max-exchanges must be positive.")
true_lambda <- sigma_a2 / sigma_e2
lambda_override <- trimws(get_arg("--lambda-grid", ""))
if (nzchar(lambda_override)) {
  lambda_grid <- sort(unique(csv_num(lambda_override)))
} else {
  multipliers <- csv_num(get_arg("--lambda-multipliers", "0,0.25,0.5,1,2,4"))
  lambda_grid <- sort(unique(true_lambda * multipliers))
}
if (any(!is.finite(lambda_grid)) || any(lambda_grid < 0)) {
  stop("All lambda values must be finite and non-negative.")
}
if (!any(abs(lambda_grid - true_lambda) < 1e-12)) {
  lambda_grid <- sort(unique(c(lambda_grid, true_lambda)))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(out_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

tasks <- lapply(seq_len(B), function(rep) list(
  case = case,
  rep = rep,
  seed = seed0 + rep * 1009L,
  R = R_groups,
  N_total = N_total,
  n = n_sub,
  sigma_a2 = sigma_a2,
  sigma_e2 = sigma_e2,
  true_lambda = true_lambda,
  method = method,
  lambda_grid = lambda_grid,
  max_exchanges = max_exchanges,
  tolerance = tolerance
))

run_or_load <- function(task) {
  checkpoint <- file.path(checkpoint_dir, sprintf("rep_%04d.csv", task$rep))
  if (resume && file.exists(checkpoint)) {
    cached <- read.csv(checkpoint, stringsAsFactors = FALSE)
    expected <- sort(task$lambda_grid)
    valid_cache <- nrow(cached) == length(expected) &&
      identical(as.character(cached$analysis_version[1]), analysis_version) &&
      identical(as.integer(cached$seed[1]), as.integer(task$seed)) &&
      identical(as.character(cached$case[1]), as.character(task$case)) &&
      identical(as.character(cached$method[1]), as.character(task$method)) &&
      isTRUE(all.equal(cached$sigma_a2_data[1], task$sigma_a2, tolerance = 1e-12)) &&
      isTRUE(all.equal(cached$sigma_e2_data[1], task$sigma_e2, tolerance = 1e-12)) &&
      identical(as.integer(cached$optimizer_max_exchanges[1]), as.integer(task$max_exchanges)) &&
      isTRUE(all.equal(cached$optimizer_tolerance[1], task$tolerance, tolerance = 1e-15)) &&
      isTRUE(all.equal(sort(cached$lambda), expected, tolerance = 1e-12))
    if (valid_cache) {
      return(cached)
    }
    unlink(checkpoint)
  }
  result <- run_replication(task)
  temporary <- paste0(checkpoint, ".", Sys.getpid(), ".tmp")
  write.csv(result, temporary, row.names = FALSE)
  if (!file.rename(temporary, checkpoint)) {
    unlink(temporary)
    stop("Could not atomically save checkpoint: ", checkpoint)
  }
  result
}

jobs <- min(jobs, length(tasks))
if (.Platform$OS.type == "windows") jobs <- 1L
message(sprintf(
  "Running %d replications over %d lambda values with NJOBS=%d",
  B, length(lambda_grid), jobs
))
raw_list <- parallel::mclapply(
  tasks,
  run_or_load,
  mc.cores = jobs,
  mc.preschedule = FALSE
)
raw <- do.call(rbind, raw_list)
raw$lambda <- vapply(
  raw$lambda,
  function(l) lambda_grid[which.min(abs(lambda_grid - l))],
  numeric(1)
)
raw$lambda_multiplier <- if (true_lambda == 0) NA_real_ else raw$lambda / true_lambda
raw$assumed_icc <- raw$lambda / (1 + raw$lambda)
raw <- raw[order(raw$rep, raw$lambda), , drop = FALSE]
summary <- summarize_results(raw, true_lambda)

prefix <- paste0(case, "_weighted_sequential_vs_joint")
raw_path <- file.path(out_dir, paste0(prefix, "_raw.csv"))
summary_path <- file.path(out_dir, paste0(prefix, "_summary.csv"))
config_path <- file.path(out_dir, paste0(prefix, "_config.csv"))
session_path <- file.path(out_dir, paste0(prefix, "_sessionInfo.txt"))

write.csv(raw, raw_path, row.names = FALSE)
write.csv(summary, summary_path, row.names = FALSE)
write.csv(data.frame(
  analysis_version = analysis_version,
  case = case,
  B = B,
  R = R_groups,
  N_total = N_total,
  n_sub = n_sub,
  sigma_a2 = sigma_a2,
  sigma_e2 = sigma_e2,
  true_lambda = true_lambda,
  lambda_grid = paste(format(lambda_grid, digits = 16), collapse = ","),
  method = method,
  jobs = jobs,
  seed = seed0,
  max_exchanges = max_exchanges,
  tolerance = tolerance,
  resume = resume,
  stringsAsFactors = FALSE
), config_path, row.names = FALSE)
writeLines(capture.output(sessionInfo()), session_path)

print(summary, row.names = FALSE)
message("Raw results: ", raw_path)
message("Summary: ", summary_path)
message("Config: ", config_path)
