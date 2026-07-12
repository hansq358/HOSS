# Small utilities for checking Lin Wang's balanced subsampling algorithm.
# The paper's method is named "balanced subsampling". Algorithm 1 is the
# sequential selection algorithm used below.

as_level_matrix <- function(x_cat) {
  x_df <- as.data.frame(x_cat, stringsAsFactors = FALSE)
  if (ncol(x_df) == 0) stop("x_cat must contain at least one categorical column")

  factors <- lapply(x_df, function(col) {
    if (is.factor(col)) col else factor(col)
  })
  q_levels <- vapply(factors, nlevels, integer(1))
  x_mat <- do.call(cbind, lapply(factors, as.integer))
  colnames(x_mat) <- names(x_df)
  attr(x_mat, "q_levels") <- q_levels
  x_mat
}

sample_values <- function(x, size) {
  if (size == 0) return(x[integer(0)])
  x[sample.int(length(x), size)]
}

same_level_delta2 <- function(candidate_mat, selected_row, q_levels) {
  candidate_mat <- as.matrix(candidate_mat)
  selected_row <- as.integer(selected_row)
  score <- numeric(nrow(candidate_mat))
  for (j in seq_along(q_levels)) {
    score <- score + q_levels[j] * (candidate_mat[, j] == selected_row[j])
  }
  as.numeric(score^2)
}

wang_balanced_indices <- function(x_cat, n, start = NULL, seed = NULL,
                                  tie = c("first", "random")) {
  tie <- match.arg(tie)
  x_mat <- as_level_matrix(x_cat)
  q_levels <- attr(x_mat, "q_levels")
  n_total <- nrow(x_mat)
  if (n > n_total) stop("n cannot exceed the number of rows in x_cat")
  if (!is.null(seed)) set.seed(seed)

  if (is.null(start)) {
    selected <- sample.int(n_total, 1)
  } else {
    if (start < 1 || start > n_total) stop("start must be a valid row index")
    selected <- as.integer(start)
  }

  candidates <- setdiff(seq_len(n_total), selected)
  phi <- rep(Inf, n_total)
  if (length(candidates) > 0) {
    phi[candidates] <- same_level_delta2(
      x_mat[candidates, , drop = FALSE],
      x_mat[selected[1], ],
      q_levels
    )
  }

  while (length(selected) < n) {
    best_score <- min(phi[candidates])
    best <- candidates[phi[candidates] == best_score]
    next_row <- if (tie == "random") sample_values(best, 1) else best[1]
    selected <- c(selected, next_row)
    candidates <- candidates[candidates != next_row]

    if (length(candidates) > 0) {
      phi[candidates] <- phi[candidates] + same_level_delta2(
        x_mat[candidates, , drop = FALSE],
        x_mat[next_row, ],
        q_levels
      )
    }
  }

  selected
}

combination_quota_indices <- function(x_cat, n, seed = NULL) {
  x_mat <- as_level_matrix(x_cat)
  n_total <- nrow(x_mat)
  if (n > n_total) stop("n cannot exceed the number of rows in x_cat")
  if (!is.null(seed)) set.seed(seed)

  keys <- apply(x_mat, 1, paste, collapse = "_")
  buckets <- split(seq_len(n_total), keys)
  out <- integer(0)
  quota <- max(1, floor(n / length(buckets)))

  for (bucket in buckets[sample.int(length(buckets))]) {
    take <- min(length(bucket), quota, n - length(out))
    if (take > 0) out <- c(out, sample_values(bucket, take))
    if (length(out) >= n) break
  }
  if (length(out) < n) {
    remaining <- setdiff(seq_len(n_total), out)
    out <- c(out, sample_values(remaining, n - length(out)))
  }
  out
}

wang_balance_components <- function(x_cat, indices, ordered_pairs = TRUE) {
  x_mat <- as_level_matrix(x_cat)
  q_levels <- attr(x_mat, "q_levels")
  xs <- x_mat[indices, , drop = FALSE]
  n <- nrow(xs)
  p <- ncol(xs)

  main_sum <- 0
  for (j in seq_len(p)) {
    counts <- tabulate(xs[, j], nbins = q_levels[j])
    main_sum <- main_sum + sum(q_levels[j]^2 * (1 / q_levels[j] - counts / n)^2)
  }

  pair_sum <- 0
  if (p > 1) {
    pair_grid <- if (ordered_pairs) {
      which(row(diag(p)) != col(diag(p)), arr.ind = TRUE)
    } else {
      which(upper.tri(matrix(FALSE, p, p)), arr.ind = TRUE)
    }

    for (r in seq_len(nrow(pair_grid))) {
      j <- pair_grid[r, 1]
      k <- pair_grid[r, 2]
      counts <- table(
        factor(xs[, j], levels = seq_len(q_levels[j])),
        factor(xs[, k], levels = seq_len(q_levels[k]))
      )
      target <- 1 / (q_levels[j] * q_levels[k])
      pair_sum <- pair_sum + sum(q_levels[j] * q_levels[k] * (target - counts / n)^2)
    }
  }

  f2 <- main_sum + pair_sum
  data.frame(
    selected_n = n,
    f = sqrt(f2),
    f2 = f2,
    main_balance = main_sum,
    pair_balance = pair_sum
  )
}

dummy_design_diagnostics <- function(x_cat, indices = NULL) {
  x_df <- as.data.frame(x_cat, stringsAsFactors = FALSE)
  x_df <- data.frame(lapply(x_df, function(col) if (is.factor(col)) col else factor(col)))
  if (!is.null(indices)) x_df <- x_df[indices, , drop = FALSE]
  design <- model.matrix(~ ., data = x_df)
  rank <- qr(design)$rank
  cond <- tryCatch(kappa(design, exact = FALSE), error = function(e) NA_real_)
  data.frame(
    design_rows = nrow(design),
    design_cols = ncol(design),
    design_rank = rank,
    full_rank = rank == ncol(design),
    condition_number = cond
  )
}

compare_balanced_samplers <- function(x_cat, n, reps = 100, seed = 1,
                                      methods = c("WANG_BAL", "COMB_QUOTA", "UNI"),
                                      tie = "first") {
  rows <- list()
  pos <- 1
  n_total <- nrow(as.data.frame(x_cat))

  for (rep_id in seq_len(reps)) {
    for (method in methods) {
      rep_seed <- seed + rep_id * 1009 + nchar(method)
      indices <- switch(
        method,
        WANG_BAL = wang_balanced_indices(x_cat, n, seed = rep_seed, tie = tie),
        COMB_QUOTA = combination_quota_indices(x_cat, n, seed = rep_seed),
        UNI = {
          set.seed(rep_seed)
          sample.int(n_total, n)
        },
        stop("unknown method: ", method)
      )

      rows[[pos]] <- cbind(
        data.frame(method = method, rep = rep_id),
        wang_balance_components(x_cat, indices),
        dummy_design_diagnostics(x_cat, indices)
      )
      pos <- pos + 1
    }
  }

  do.call(rbind, rows)
}

summarize_balanced_comparison <- function(raw) {
  metric_names <- c(
    "selected_n", "f", "f2", "main_balance", "pair_balance",
    "design_rank", "full_rank", "condition_number"
  )
  groups <- split(raw, raw$method)
  rows <- lapply(groups, function(df) {
    out <- data.frame(method = df$method[1], reps = nrow(df))
    for (metric in metric_names) out[[metric]] <- mean(df[[metric]], na.rm = TRUE)
    out
  })
  do.call(rbind, rows)
}

make_full_factorial_categorical <- function(q_levels, repeats = 1,
                                            shuffle = TRUE, seed = NULL) {
  factors <- expand.grid(lapply(q_levels, seq_len))
  names(factors) <- paste0("x", seq_along(q_levels))
  out <- factors[rep(seq_len(nrow(factors)), each = repeats), , drop = FALSE]
  if (shuffle) {
    if (!is.null(seed)) set.seed(seed)
    out <- out[sample.int(nrow(out)), , drop = FALSE]
  }
  rownames(out) <- NULL
  for (j in seq_along(q_levels)) {
    out[[j]] <- factor(out[[j]], levels = seq_len(q_levels[j]))
  }
  out
}

make_wang_simulation_categorical <- function(case = "case1", p = 20,
                                             n_total = 5000, seed = 1) {
  case <- match.arg(case, c("case1", "case2", "case3"))
  if (!is.null(seed)) set.seed(seed)
  q_levels <- seq_len(p) + 1
  out <- vector("list", p)

  if (case == "case1") {
    for (j in seq_len(p)) {
      out[[j]] <- sample.int(q_levels[j], n_total, replace = TRUE)
    }
  } else if (case == "case2") {
    for (j in seq_len(p)) {
      probs <- seq_len(q_levels[j])
      out[[j]] <- sample.int(q_levels[j], n_total, replace = TRUE, prob = probs)
    }
  } else {
    sigma <- matrix(0.5, p, p)
    diag(sigma) <- 1
    normal_draws <- matrix(rnorm(n_total * p), n_total, p) %*% chol(sigma)
    for (j in seq_len(p)) {
      breaks <- seq(-3, 3, length.out = q_levels[j] + 1)
      out[[j]] <- pmin(q_levels[j], pmax(1, findInterval(normal_draws[, j], breaks)))
    }
  }

  out <- as.data.frame(out)
  names(out) <- paste0("x", seq_len(p))
  for (j in seq_len(p)) out[[j]] <- factor(out[[j]], levels = seq_len(q_levels[j]))
  out
}

parse_cli_args <- function(args) {
  get_arg <- function(name, default = NULL) {
    eq <- paste0(name, "=")
    hit <- grep(paste0("^", eq), args)
    if (length(hit) > 0) return(sub(eq, "", args[hit[1]], fixed = TRUE))
    pos <- match(name, args)
    if (!is.na(pos) && pos < length(args)) return(args[pos + 1])
    default
  }
  csv_chr <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  list(
    get = get_arg,
    csv_chr = csv_chr,
    csv_int = function(x) as.integer(csv_chr(x))
  )
}

write_balanced_outputs <- function(raw, out_dir, prefix) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(raw, file.path(out_dir, paste0(prefix, "_raw.csv")), row.names = FALSE)
  summary <- summarize_balanced_comparison(raw)
  write.csv(summary, file.path(out_dir, paste0(prefix, "_summary.csv")), row.names = FALSE)
  summary
}

run_wang_balanced_smoke <- function(seed = 1, reps = 50) {
  x_cat <- make_full_factorial_categorical(c(5, 5), repeats = 2, seed = seed)
  raw <- compare_balanced_samplers(x_cat, n = 25, reps = reps, seed = seed)
  summarize_balanced_comparison(raw)
}

run_wang_balanced_check <- function(args = commandArgs(trailingOnly = TRUE)) {
  cli <- parse_cli_args(args)
  mode <- cli$get("--mode", "paper_example")
  reps <- as.integer(cli$get("--reps", "1000"))
  seed <- as.integer(cli$get("--seed", "20260621"))
  n_sub <- as.integer(cli$get("--n-sub", ifelse(mode == "paper_example", "25", "500")))
  methods <- cli$csv_chr(cli$get("--methods", "WANG_BAL,COMB_QUOTA,UNI"))
  out_dir <- cli$get("--output-dir", ".")
  tie <- cli$get("--tie", "first")

  if (mode == "paper_example") {
    q_levels <- cli$csv_int(cli$get("--q-levels", "5,5"))
    repeats <- as.integer(cli$get("--factorial-repeats", "2"))
    x_cat <- make_full_factorial_categorical(q_levels, repeats = repeats, seed = seed)
    prefix <- paste0("wang_balanced_paper_example_n", n_sub)
  } else if (mode == "wang_case") {
    case <- cli$get("--case", "case1")
    p <- as.integer(cli$get("--p", "20"))
    n_total <- as.integer(cli$get("--N-total", "5000"))
    x_cat <- make_wang_simulation_categorical(case = case, p = p, n_total = n_total, seed = seed)
    prefix <- paste0("wang_balanced_", case, "_N", n_total, "_n", n_sub)
  } else {
    stop("unknown --mode: ", mode, ". Use paper_example or wang_case.")
  }

  raw <- compare_balanced_samplers(
    x_cat = x_cat,
    n = n_sub,
    reps = reps,
    seed = seed,
    methods = methods,
    tie = tie
  )
  summary <- write_balanced_outputs(raw, out_dir, prefix)
  print(summary)
  invisible(summary)
}
