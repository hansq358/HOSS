# Comparison with balanced subsampling within the HOSS-D stage, Case 3
# (Supplementary Table S17). Within-group rules: HOSS-I, an adaptation of
# the balanced subsampling of Wang (2026), and uniform sampling. The
# balanced algorithm itself is in balanced_subsampling_check.R.
# Run from this folder:  Rscript Balanced_Subsampling.R
source("common.R")
source("balanced_subsampling_check.R")

select_within <- function(pop, selected_groups, n, w_strat_name) {
  m <- ceiling(n / pop$R)
  unlist(sapply(selected_groups, function(g_idx) {
    rows <- (pop$C_cum[g_idx] + 1):pop$C_cum[g_idx + 1]
    size <- min(m, pop$C_i[g_idx])
    if (size == 0) return(NULL)
    if (w_strat_name == "GOSS") {
      return(rows[OAJ2_cpp(apply(pop$X_within_select[rows, , drop = FALSE], 2, scalex), size, tPow = 2)])
    }
    if (w_strat_name == "GUNIF") {
      return(rows[sample.int(pop$C_i[g_idx], size, replace = FALSE)])
    }
    if (w_strat_name == "GLEV") {
      return(rows[sample.int(pop$C_i[g_idx], size, replace = FALSE, prob = lev(pop$X_within_select_std[rows, , drop = FALSE]))])
    }
    if (w_strat_name == "WANGBAL") {
      if (ncol(pop$X_ind_cat_raw) == 0) return(rows[sample.int(pop$C_i[g_idx], size, replace = FALSE)])
      local_cat <- pop$X_ind_cat_raw[rows, , drop = FALSE]
      for (nm in names(local_cat)) {
        local_cat[[nm]] <- factor(local_cat[[nm]], levels = sort(unique(pop$X_ind_cat_raw[[nm]])))
      }
      local_idx <- wang_balanced_indices(local_cat, size, tie = "first")
      return(rows[local_idx])
    }
    rows[iboss(pop$X_within_select_std[rows, , drop = FALSE], size)]
  }))
}

run_hoss_wang_balanced <- function() {
  case <- get_arg("--case", "case3")
  B <- as.integer(get_arg("--B", "300"))
  R <- as.integer(get_arg("--R", "100"))
  N_total <- as.integer(get_arg("--N-total", "50000"))
  n <- as.integer(get_arg("--n-sub", "1000"))
  sigma_a2 <- as.numeric(get_arg("--sigma-a2", "0.5"))
  sigma_e2 <- as.numeric(get_arg("--sigma-e2", "9"))
  jobs <- as.integer(get_arg("--jobs", Sys.getenv("NJOBS", "128")))
  seed0 <- as.integer(get_arg("--seed", "20260603"))
  out_dir <- get_arg("--output-dir", file.path(root, "results"))
  methods <- methods_from_arg(
    get_arg("--methods", "Dopt_0.9_GOSS,Dopt_0.9_WANGBAL,Unif_0.9_GUNIF"),
    c(0.9, 0.95)
  )

  tasks <- lapply(seq_len(B), function(rep) {
    list(
      case = case,
      rep = rep,
      seed = seed0 + rep * 1009,
      R = R,
      N_total = N_total,
      n = n,
      sigma_a2 = sigma_a2,
      sigma_e2 = sigma_e2,
      methods = methods,
      truth_model = "full",
      fit_model = "full",
      scenario = "wang_balanced"
    )
  })

  raw <- collect(tasks, task_prediction, jobs)
  prefix <- paste0(case, "_hoss_wang_balanced")
  write_outputs(out_dir, prefix, raw, c("case", "method"))
}

run_hoss_wang_balanced()
