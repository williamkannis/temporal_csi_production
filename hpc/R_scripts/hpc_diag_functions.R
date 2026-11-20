tweedie_diagnostics <- function(
    fit,
    y,
    yr,
    st,
    rg,
    data = NULL,
    ndraws = 1000,
    seed = 123,
    ppc_groups = NULL
) {
  

  # Packages  ------------------------------------------------------------------
  
  requireNamespace("ggplot2")
  requireNamespace("bayesplot")
  set.seed(seed)
  

  # Helpers  -------------------------------------------------------------------
  
  # extract posterior draws from rstan OR cmdstanr
  extract_draws <- function(fit) {
    
    # cmdstanr
    if (inherits(fit, "CmdStanMCMC")) {
      
      draws <- posterior::as_draws_matrix(
        fit$draws()
      )
      
      # rstan
    } else if (inherits(fit, "stanfit")) {
      
      draws <- as.matrix(fit)
      
    } else {
      
      stop(
        "fit must be either a rstan::stanfit or cmdstanr::CmdStanMCMC object."
      )
      
    }
    draws
  }
  
  # Find variables such as y_rep[1], residuals[1], etc.
  get_array <- function(draws, name) {
    
    pattern <- paste0("^", name, "\\[")
    
    cols <- grep(pattern, colnames(draws))
    
    if (length(cols) == 0) {
      return(NULL)
    }
    out <- draws[, cols, drop = FALSE]
    
    # Make sure variables are ordered numerically
    indices <- as.numeric(
      sub(
        paste0("^", name, "\\[([0-9]+)\\]$"),
        "\\1",
        colnames(out)
      )
    )
    
    out <- out[, order(indices), drop = FALSE]
    out
  }
  

  # Extract generated quantities -----------------------------------------------
  
  draws <- extract_draws(fit)
  y_rep <- get_array(draws, "y_rep")
  raw_resid <- get_array(draws, "residuals")
  log_lik <- get_array(draws, "log_lik")
  p_zero_draws <- get_array(draws, "p_zero")
  pearson_draws <- get_array(draws, "pearson_residual")
  
  
  if (is.null(y_rep)) {
    stop(
      "Could not find y_rep in the fitted model. ",
      "Make sure y_rep is generated in generated quantities."
    )
  }
  
  
  # Thin posterior draws if requested ------------------------------------------

  n_total <- nrow(y_rep)
  
  if (ndraws < n_total) {
    
    keep <- sample(
      seq_len(n_total),
      size = ndraws,
      replace = FALSE
    )
    
    y_rep <- y_rep[keep, , drop = FALSE]
    
    if (!is.null(raw_resid))
      raw_resid <- raw_resid[keep, , drop = FALSE]
    
    if (!is.null(log_lik))
      log_lik <- log_lik[keep, , drop = FALSE]
    
    if (!is.null(p_zero_draws))
      p_zero_draws <- p_zero_draws[keep, , drop = FALSE]
    
    if (!is.null(pearson_draws))
      pearson_draws <- pearson_draws[keep, , drop = FALSE]
    
    draws <- draws[keep, , drop = FALSE]
  }
  
  
  N <- length(y)
  if (ncol(y_rep) != N) {
    stop(
      "Length of y does not match number of observations in y_rep."
    )
  }
  

  # Posterior predictive quantile/rank residuals  ------------------------------

  # For each observation:
  #
  #   lower = P(Yrep < Y)
  #   upper = P(Yrep <= Y)
  #
  # Randomization handles the point mass at zero.
  # Transform to standard normal scale.

  lower <- numeric(N)
  upper <- numeric(N)
  
  for (n in seq_len(N)) {
    lower[n] <- mean(y_rep[, n] < y[n])
    upper[n] <- mean(y_rep[, n] <= y[n])
  }
  
  # Randomized probability
  u <- lower + runif(N) * (upper - lower)
  
  # Avoid +/- Inf from qnorm
  eps <- 0.5/nrow(y_rep)
  
  u <- pmin(
    pmax(u, eps),
    1 - eps
  )
  
  rq_residual <- qnorm(u)
  
  
  # Fitted values --------------------------------------------------------------

  mu_draws <- get_array(draws, "mu")
  
  if (is.null(mu_draws)) {
    stop(
      "Could not find mu in the fitted model. ",
      "Make sure mu is saved in generated quantities."
    )
  }
  
  fitted <- apply(
    mu_draws,
    2,
    median
  )
  
  
  # Raw residuals  -------------------------------------------------------------

  if (!is.null(raw_resid)) {
    
    raw_residual <- apply(
      raw_resid,
      2,
      median
    )
    
  } else {
    raw_residual <- y - fitted
  }
  
  
  # Pearson residuals  ---------------------------------------------------------
  
  if (!is.null(pearson_draws)) {
    
    pearson_residual <- apply(
      pearson_draws,
      2,
      median
    )
    
  } else {
    
    phi_cols <- grep(
      "^phi$",
      colnames(draws)
    )
    
    power_cols <- grep(
      "^power$",
      colnames(draws)
    )
    
    if (length(phi_cols) == 1 &&
        length(power_cols) == 1 &&
        !is.null(raw_resid)) {
      
      phi_draw <- draws[, phi_cols]
      power_draw <- draws[, power_cols]
      
      pearson_draws <- matrix(
        NA_real_,
        nrow = nrow(draws),
        ncol = N
      )
      
      for (i in seq_len(nrow(draws))) {
        
        mu_i <- y - raw_resid[i, ]
        
        pearson_draws[i, ] <-
          raw_resid[i, ] /
          sqrt(
            phi_draw[i] *
              mu_i^power_draw[i]
          )
      }
      
      pearson_residual <- apply(
        pearson_draws,
        2,
        median
      )
      
    } else {
      
      pearson_residual <- NULL
    }
  }
  
  
  # Probability of zero  -------------------------------------------------------
  
  if (!is.null(p_zero_draws)) {
    
    p_zero <- apply(
      p_zero_draws,
      2,
      median
    )
    
  } else {
    
    phi_cols <- grep("^phi$", colnames(draws))
    power_cols <- grep("^power$", colnames(draws))
    
    if (length(phi_cols) == 1 &&
        length(power_cols) == 1 &&
        !is.null(raw_resid)) {
      
      phi_draw <- draws[, phi_cols]
      power_draw <- draws[, power_cols]
      
      p_zero_draws <- matrix(
        NA_real_,
        nrow = nrow(draws),
        ncol = N
      )
      
      for (i in seq_len(nrow(draws))) {
        
        mu_i <- y - raw_resid[i, ]
        
        lambda_i <-
          mu_i^(2 - power_draw[i]) /
          (
            phi_draw[i] *
              (2 - power_draw[i])
          )
        
        p_zero_draws[i, ] <- exp(-lambda_i)
      }
      
      p_zero <- apply(
        p_zero_draws,
        2,
        median
      )
      
    } else {
      
      p_zero <- NULL
    }
  }
  
  
  # Observation-level diagnostic data frame  -----------------------------------

  diagnostic_data <- data.frame(
    observation = seq_len(N),
    y = y,
    fitted = fitted,
    raw_residual = raw_residual,
    rq_residual = rq_residual,
    observed_zero = y == 0
  )
  
  if (!is.null(pearson_residual)) {
    diagnostic_data$pearson_residual <-
      pearson_residual
  }
  
  if (!is.null(p_zero)) {
    diagnostic_data$p_zero <-
      p_zero
  }
  

  # Posterior predictive summary statistics

  observed_stats <- c(
    mean = mean(y),
    median = median(y),
    variance = var(y),
    sd = sd(y),
    zero_proportion = mean(y == 0),
    q90 = quantile(y, 0.90),
    q95 = quantile(y, 0.95),
    q99 = quantile(y, 0.99),
    maximum = max(y)
  )
  
  replicated_stats <- data.frame(
    mean = apply(y_rep, 1, mean),
    median = apply(y_rep, 1, median),
    variance = apply(y_rep, 1, var),
    sd = apply(y_rep, 1, sd),
    zero_proportion =
      apply(y_rep, 1, function(x)
        mean(x == 0)
      ),
    q90 = apply(y_rep, 1, quantile, probs = 0.90),
    q95 = apply(y_rep, 1, quantile, probs = 0.95),
    q99 = apply(y_rep, 1, quantile, probs = 0.99),
    maximum = apply(y_rep, 1, max)
  )
  

  # Bayesian posterior predictive p-values  ------------------------------------

  bayesian_p <- sapply(
    names(observed_stats),
    function(stat) {
      mean(
        replicated_stats[[stat]] >=
          observed_stats[[stat]]
      )
    }
  )
  
  
  # Also calculate two-sided-style discrepancy
  # relative to posterior predictive median
  bayesian_p_two_sided <- sapply(
    names(observed_stats),
    function(stat) {
      
      obs <- observed_stats[[stat]]
      rep <- replicated_stats[[stat]]
      
      median_rep <- median(rep)
      
      mean(
        abs(rep - median_rep) >=
          abs(obs - median_rep)
      )
    }
  )
  
  ppc_summary <- data.frame(
    statistic = names(observed_stats),
    observed = as.numeric(observed_stats),
    predictive_mean =
      sapply(
        replicated_stats,
        mean
      ),
    predictive_median =
      sapply(
        replicated_stats,
        median
      ),
    predictive_lower =
      sapply(
        replicated_stats,
        quantile,
        probs = 0.025
      ),
    predictive_upper =
      sapply(
        replicated_stats,
        quantile,
        probs = 0.975
      ),
    bayesian_p =
      as.numeric(bayesian_p),
    bayesian_p_two_sided =
      as.numeric(bayesian_p_two_sided),
    row.names = NULL
  )
  
  
  # Overall PPC plots  ---------------------------------------------------------

  # Sample posterior predictive draws for visualization
  n_plot <- min(100, nrow(y_rep))
  
  plot_draws <- y_rep[
    sample(
      seq_len(nrow(y_rep)),
      n_plot
    ),
    ,
    drop = FALSE
  ]
  
  
  # Density overlay
  ppc_density <- bayesplot::ppc_dens_overlay(
    y = y,
    yrep = plot_draws
  ) +
    ggplot2::labs(
      title = "Posterior predictive distribution",
      x = "Response",
      y = "Density"
    )
  
  
  # ECDF overlay
  ppc_ecdf <- bayesplot::ppc_ecdf_overlay(
    y = y,
    yrep = plot_draws
  ) +
    ggplot2::labs(
      title = "Posterior predictive ECDF",
      x = "Response",
      y = "Cumulative probability"
    )
  
  
  # PPC statistics plots  ------------------------------------------------------
  
  ppc_stat_plot <- bayesplot::ppc_stat(
    y = y,
    yrep = plot_draws,
    stat = "mean"
  ) +
    ggplot2::labs(
      title = "Posterior predictive check: mean",
      x = "Mean"
    )
  
  
  ppc_zero_plot <- bayesplot::ppc_stat(
    y = y,
    yrep = plot_draws,
    stat = function(x) mean(x == 0)
  ) +
    ggplot2::labs(
      title = "Posterior predictive check: zero proportion",
      x = "Proportion of zeros"
    )
  
  
  ppc_variance_plot <- bayesplot::ppc_stat(
    y = y,
    yrep = plot_draws,
    stat = var
  ) +
    ggplot2::labs(
      title = "Posterior predictive check: variance",
      x = "Variance"
    )
  
  
  ppc_max_plot <- bayesplot::ppc_stat(
    y = y,
    yrep = plot_draws,
    stat = max
  ) +
    ggplot2::labs(
      title = "Posterior predictive check: maximum",
      x = "Maximum"
    )
  

  # RQR QQ plot  ---------------------------------------------------------------

  qq_data <- data.frame(
    rq_residual = rq_residual
  )
  
  qq_rqr <- ggplot2::ggplot(
    qq_data,
    ggplot2::aes(sample = rq_residual)
  ) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line() +
    ggplot2::labs(
      title = "Posterior predictive quantile residuals",
      x = "Theoretical quantiles",
      y = "Posterior predictive quantile residuals"
    ) +
    ggplot2::theme_classic()
  

  # RQR vs fitted  -------------------------------------------------------------
  
  rqr_vs_fitted <- ggplot2::ggplot(
    diagnostic_data,
    ggplot2::aes(
      x = fitted,
      y = rq_residual
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      se = TRUE
    ) +
    ggplot2::labs(
      title = "Quantile residuals vs fitted values",
      x = "Fitted value",
      y = "Quantile residual"
    ) +
    ggplot2::theme_classic()
  
  
  # Pearson residual vs fitted  ------------------------------------------------

  pearson_vs_fitted <- NULL
  
  if (!is.null(pearson_residual)) {
    
    pearson_vs_fitted <- ggplot2::ggplot(
      diagnostic_data,
      ggplot2::aes(
        x = fitted,
        y = pearson_residual
      )
    ) +
      ggplot2::geom_point(
        alpha = 0.5
      ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed"
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        se = TRUE
      ) +
      ggplot2::labs(
        title = "Pearson residuals vs fitted values",
        x = "Fitted value",
        y = "Pearson residual"
      ) +
      ggplot2::theme_classic()
  }
  
  
  # Raw residual vs fitted  ----------------------------------------------------
  
  raw_vs_fitted <- ggplot2::ggplot(
    diagnostic_data,
    ggplot2::aes(
      x = fitted,
      y = raw_residual
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      se = TRUE
    ) +
    ggplot2::labs(
      title = "Raw residuals vs fitted values",
      x = "Fitted value",
      y = "Raw residual"
    ) +
    ggplot2::theme_classic()
  
  
  # Residuals vs observation order  --------------------------------------------
  
  residual_order <- ggplot2::ggplot(
    diagnostic_data,
    ggplot2::aes(
      x = observation,
      y = rq_residual
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      se = TRUE
    ) +
    ggplot2::labs(
      title = "Quantile residuals vs observation order",
      x = "Observation",
      y = "Quantile residual"
    ) +
    ggplot2::theme_classic()
  
  
  # Zero probability diagnostic  -----------------------------------------------
  
  zero_plot <- NULL
  
  if (!is.null(p_zero)) {
    
    zero_plot <- ggplot2::ggplot(
      diagnostic_data,
      ggplot2::aes(
        x = p_zero,
        y = as.numeric(observed_zero)
      )
    ) +
      ggplot2::geom_jitter(
        height = 0.05,
        width = 0,
        alpha = 0.4
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        se = TRUE
      ) +
      ggplot2::labs(
        title = "Observed zeros vs predicted zero probability",
        x = "Predicted probability of zero",
        y = "Observed zero"
      ) +
      ggplot2::theme_classic()
  }
  
  
  # Region × year-level posterior predictive diagnostics  ----------------------
  
  # Create region × year grouping
  region_year <- interaction(
    rg[st],
    yr,
    drop = TRUE,
    sep = "_"
  )
  
  group_levels <- levels(region_year)
  
  group_diagnostics <- lapply(group_levels, function(g) {
    
    idx <- which(region_year == g)
    
    # Observed values
    y_obs <- y[idx]
    
    # Fitted values from posterior mu
    mu_group <- rowMeans(mu_draws[, idx, drop = FALSE])
    
    # Posterior predictive group means
    yrep_group <- rowMeans(y_rep[, idx, drop = FALSE])
    
    # Posterior predictive zero proportions
    yrep_zero <- rowMeans(y_rep[, idx, drop = FALSE] == 0)
    
    # Observed statistics
    obs_mean <- mean(y_obs)
    obs_zero <- mean(y_obs == 0)
    
    # Posterior predictive summaries
    pred_mean <- median(yrep_group)
    pred_mean_lower <- quantile(yrep_group, 0.025)
    pred_mean_upper <- quantile(yrep_group, 0.975)
    
    pred_zero <- median(yrep_zero)
    pred_zero_lower <- quantile(yrep_zero, 0.025)
    pred_zero_upper <- quantile(yrep_zero, 0.975)
    
    # Posterior predictive p-values
    p_mean <- mean(yrep_group >= obs_mean)
    p_zero <- mean(yrep_zero >= obs_zero)
    
    data.frame(
      group = g,
      n = length(idx),
      observed_mean = obs_mean,
      fitted_mean = median(mu_group),
      predicted_mean = pred_mean,
      predicted_mean_lower = pred_mean_lower,
      predicted_mean_upper = pred_mean_upper,
      observed_zero = obs_zero,
      predicted_zero = pred_zero,
      predicted_zero_lower = pred_zero_lower,
      predicted_zero_upper = pred_zero_upper,
      p_mean = p_mean,
      p_zero = p_zero
    )
  })
  
  group_diagnostics <- do.call(rbind, group_diagnostics)
  
  
  # observed vs posterior predictive region × year means plot  -----------------
  
  group_mean_plot <- ggplot2::ggplot(
    group_diagnostics,
    ggplot2::aes(x = predicted_mean, y = observed_mean)
  ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = predicted_mean_lower,
        xmax = predicted_mean_upper
      ),
      height = 0
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::labs(
      x = "Posterior predictive mean",
      y = "Observed mean",
      title = "Region × year posterior predictive check"
    ) +
    ggplot2::theme_classic()
  
  
  # observed vs posterior predictive zero proportions plot  --------------------
  
  group_zero_plot <- ggplot2::ggplot(
    group_diagnostics,
    ggplot2::aes(x = predicted_zero, y = observed_zero)
  ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = predicted_zero_lower,
        xmax = predicted_zero_upper
      ),
      height = 0
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::labs(
      x = "Posterior predictive zero proportion",
      y = "Observed zero proportion",
      title = "Region × year zero-inflation check"
    ) +
    ggplot2::theme_classic()
  
  
  # Return a single diagnostic object  -----------------------------------------

  diagnostics <- list(
    
    # Observation-level diagnostics
    data = diagnostic_data,
    
    # Posterior predictive simulations
    y_rep = y_rep,
    
    # Posterior predictive summary statistics
    ppc = list(
      observed = observed_stats,
      replicated = replicated_stats,
      summary = ppc_summary
    ),
    
    # Group posterior predictive summary statistics
    group_ppc = list(
      summary = group_diagnostics
      ),
      
    # Bayesian p-values
    bayesian_p = bayesian_p,
    bayesian_p_two_sided =
      bayesian_p_two_sided,
    
    # Residuals
    residuals = list(
      raw = raw_residual,
      pearson = pearson_residual,
      quantile = rq_residual
    ),
    
    # Predicted zero probabilities
    p_zero = p_zero,
    
    # Plots
    plots = list(
      rqr_qq = qq_rqr,
      rqr_vs_fitted = rqr_vs_fitted,
      pearson_vs_fitted = pearson_vs_fitted,
      raw_vs_fitted = raw_vs_fitted,
      residual_order = residual_order,
      zero_probability = zero_plot,
      zero_probability_group = group_zero_plot,
      group_mean = group_mean_plot,
      ppc_density = ppc_density,
      ppc_ecdf = ppc_ecdf,
      ppc_mean = ppc_stat_plot,
      ppc_zero = ppc_zero_plot,
      ppc_variance = ppc_variance_plot,
      ppc_max = ppc_max_plot
    ),
    
    # Settings
    settings = list(
      ndraws = nrow(y_rep),
      seed = seed
    )
  )
  
  
  class(diagnostics) <- "tweedie_diagnostics"
  
  return(diagnostics)
}


hurdle_diagnostics <- function(
    fit,
    y,
    yr,
    st,
    rg,
    data = NULL,
    ndraws = 1000,
    seed = 123,
    ppc_groups = NULL
) {
  
  # Packages ------------------------------------------------------------------
  
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")
  
  if (!requireNamespace("bayesplot", quietly = TRUE))
    stop("Package 'bayesplot' is required.")
  
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Package 'posterior' is required.")
  
  set.seed(seed)
  
  
  # Input validation ----------------------------------------------------------
  
  if (!is.numeric(y) || anyNA(y) || any(!is.finite(y)))
    stop("y must be a numeric vector without missing or non-finite values.")
  
  if (any(y < 0))
    stop("y must contain non-negative values.")
  
  N <- length(y)
  
  if (length(yr) != N || length(st) != N)
    stop("yr and st must have the same length as y.")
  
  if (length(rg) < max(st))
    stop("rg must contain a region assignment for every site.")
  
  if (anyNA(yr) || anyNA(st) || anyNA(rg))
    stop("yr, st, and rg cannot contain missing values.")
  
  if (any(st < 1) || any(st != as.integer(st)))
    stop("st must contain positive integer site indices.")
  
  if (any(yr < 1) || any(yr != as.integer(yr)))
    stop("yr must contain positive integer year indices.")
  
  if (length(ndraws) != 1 || is.na(ndraws) ||
      ndraws < 1 || ndraws != as.integer(ndraws))
    stop("ndraws must be a positive integer.")
  
  
  # Helpers -------------------------------------------------------------------
  
  # Extract posterior draws from rstan OR cmdstanr
  
  extract_draws <- function(fit) {
    
    # cmdstanr
    if (inherits(fit, "CmdStanMCMC")) {
      
      draws <- posterior::as_draws_matrix(
        fit$draws()
      )
      
      # rstan
    } else if (inherits(fit, "stanfit")) {
      
      draws <- as.matrix(fit)
      
    } else {
      
      stop(
        "fit must be either a rstan::stanfit or ",
        "cmdstanr::CmdStanMCMC object."
      )
      
    }
    
    draws
  }
  
  
  # Extract array variables such as y_rep[1], mu[1], etc.
  
  get_array <- function(draws, name) {
    
    pattern <- paste0("^", name, "\\[")
    
    cols <- grep(pattern, colnames(draws))
    
    if (length(cols) == 0)
      return(NULL)
    
    out <- draws[, cols, drop = FALSE]
    
    indices <- as.numeric(
      sub(
        paste0("^", name, "\\[([0-9]+)\\]$"),
        "\\1",
        colnames(out)
      )
    )
    
    out <- out[, order(indices), drop = FALSE]
    
    out
  }
  
  
  # Extract scalar parameters
  
  get_scalar <- function(draws, name) {
    
    cols <- which(colnames(draws) == name)
    
    if (length(cols) != 1)
      return(NULL)
    
    draws[, cols]
  }
  
  
  # Posterior predictive summary statistics
  
  calculate_stats <- function(x) {
    
    c(
      mean = mean(x),
      median = median(x),
      variance = var(x),
      sd = sd(x),
      q90 = unname(quantile(x, 0.90)),
      q95 = unname(quantile(x, 0.95)),
      q99 = unname(quantile(x, 0.99)),
      maximum = max(x)
    )
  }
  
  
  # Summarize replicated statistics
  
  summarize_ppc <- function(observed_stats, replicated_stats) {
    
    bayesian_p <- sapply(
      names(observed_stats),
      function(stat) {
        mean(
          replicated_stats[[stat]] >=
            observed_stats[[stat]]
        )
      }
    )
    
    bayesian_p_two_sided <- sapply(
      names(observed_stats),
      function(stat) {
        
        obs <- observed_stats[[stat]]
        rep <- replicated_stats[[stat]]
        
        median_rep <- median(rep)
        
        mean(
          abs(rep - median_rep) >=
            abs(obs - median_rep)
        )
      }
    )
    
    summary <- data.frame(
      statistic = names(observed_stats),
      observed = as.numeric(observed_stats),
      
      predictive_mean = sapply(
        replicated_stats,
        mean
      ),
      
      predictive_median = sapply(
        replicated_stats,
        median
      ),
      
      predictive_lower = sapply(
        replicated_stats,
        quantile,
        probs = 0.025
      ),
      
      predictive_upper = sapply(
        replicated_stats,
        quantile,
        probs = 0.975
      ),
      
      bayesian_p = as.numeric(bayesian_p),
      
      bayesian_p_two_sided =
        as.numeric(bayesian_p_two_sided),
      
      row.names = NULL
    )
    
    list(
      summary = summary,
      bayesian_p = bayesian_p,
      bayesian_p_two_sided = bayesian_p_two_sided
    )
  }
  
  
  # Extract posterior draws ---------------------------------------------------
  
  draws <- extract_draws(fit)
  
  y_rep <- get_array(draws, "y_rep")
  
  log_lik <- get_array(draws, "log_lik")
  
  mu_draws <- get_array(draws, "mu")
  
  hu_draws <- get_array(draws, "hu")
  
  hurdle_resid_draws <- get_array(
    draws,
    "pearson_residual_hurdle"
  )
  
  gamma_resid_draws <- get_array(
    draws,
    "pearson_residual_gamma"
  )
  
  
  # Required variables
  
  if (is.null(y_rep))
    stop("Could not find y_rep in the fitted model.")
  
  if (is.null(mu_draws))
    stop("Could not find mu in the fitted model.")
  
  if (is.null(hu_draws))
    stop("Could not find hu in the fitted model.")
  
  if (is.null(get_scalar(draws, "shape")))
    stop("Could not find shape in the fitted model.")
  
  if (ncol(y_rep) != N ||
      ncol(mu_draws) != N ||
      ncol(hu_draws) != N) {
    stop(
      "The number of observations in y, y_rep, mu, and hu must match."
    )
  }
  
  if (nrow(y_rep) != nrow(draws))
    stop("Posterior draw dimensions do not match.")
  
  if (any(!is.finite(hu_draws)) ||
      any(hu_draws < 0 | hu_draws > 1))
    stop("Posterior draws of hu contain invalid probabilities.")
  
  if (any(!is.finite(mu_draws)) ||
      any(mu_draws <= 0))
    stop("Posterior draws of mu must be positive and finite.")
  
  
  # Thin posterior draws ------------------------------------------------------
  
  n_total <- nrow(y_rep)
  
  keep <- if (ndraws < n_total) {
    sample(
      seq_len(n_total),
      size = ndraws,
      replace = FALSE
    )
  } else {
    seq_len(n_total)
  }
  
  y_rep <- y_rep[keep, , drop = FALSE]
  mu_draws <- mu_draws[keep, , drop = FALSE]
  hu_draws <- hu_draws[keep, , drop = FALSE]
  
  if (!is.null(log_lik))
    log_lik <- log_lik[keep, , drop = FALSE]
  
  if (!is.null(hurdle_resid_draws))
    hurdle_resid_draws <-
    hurdle_resid_draws[keep, , drop = FALSE]
  
  if (!is.null(gamma_resid_draws))
    gamma_resid_draws <-
    gamma_resid_draws[keep, , drop = FALSE]
  
  draws <- draws[keep, , drop = FALSE]
  
  shape_draws <- get_scalar(draws, "shape")
  
  n_draws <- nrow(draws)
  
  
  # Posterior fitted values ---------------------------------------------------
  
  # Conditional Gamma mean
  
  fitted_gamma <- apply(
    mu_draws,
    2,
    median
  )
  
  # Probability of zero
  
  p_zero <- apply(
    hu_draws,
    2,
    median
  )
  
  # Unconditional expected response:
  # E[Y] = (1 - hu) * mu
  
  expected_draws <- (1 - hu_draws) * mu_draws
  
  fitted <- apply(
    expected_draws,
    2,
    median
  )
  
  
  # Generate conditional Gamma posterior predictive draws --------------------
  
  positive_idx <- which(y > 0)
  
  gamma_rep <- NULL
  
  if (length(positive_idx) > 0) {
    
    gamma_rep <- matrix(
      NA_real_,
      nrow = n_draws,
      ncol = length(positive_idx)
    )
    
    for (i in seq_len(n_draws)) {
      
      gamma_rep[i, ] <- rgamma(
        n = length(positive_idx),
        shape = shape_draws[i],
        rate = shape_draws[i] /
          mu_draws[i, positive_idx]
      )
      
    }
    
    colnames(gamma_rep) <- positive_idx
  }
  
  
  # Randomized quantile residuals ---------------------------------------------
  
  lower <- numeric(N)
  upper <- numeric(N)
  
  for (n in seq_len(N)) {
    
    lower[n] <- mean(y_rep[, n] < y[n])
    
    upper[n] <- mean(y_rep[, n] <= y[n])
    
  }
  
  u <- lower + runif(N) * (upper - lower)
  
  eps <- 0.5 / n_draws
  
  u <- pmin(
    pmax(u, eps),
    1 - eps
  )
  
  rq_residual <- qnorm(u)
  
  
  # Raw residuals -------------------------------------------------------------
  
  # Use unconditional fitted values rather than mu alone.
  
  raw_residual <- y - fitted
  
  
  # Hurdle Pearson residuals ---------------------------------------------------
  
  if (!is.null(hurdle_resid_draws)) {
    
    pearson_hurdle <- apply(
      hurdle_resid_draws,
      2,
      median,
      na.rm = TRUE
    )
    
  } else {
    
    pearson_hurdle_draws <-
      matrix(
        NA_real_,
        nrow = n_draws,
        ncol = N
      )
    
    for (i in seq_len(n_draws)) {
      
      pearson_hurdle_draws[i, ] <-
        ((y == 0) - hu_draws[i, ]) /
        sqrt(
          hu_draws[i, ] *
            (1 - hu_draws[i, ])
        )
    }
    
    pearson_hurdle <- apply(
      pearson_hurdle_draws,
      2,
      median,
      na.rm = TRUE
    )
  }
  
  
  # Gamma Pearson residuals ---------------------------------------------------
  
  if (!is.null(gamma_resid_draws)) {
    
    pearson_gamma <- apply(
      gamma_resid_draws,
      2,
      function(x) {
        
        if (all(!is.finite(x)))
          return(NA_real_)
        
        median(x[is.finite(x)])
      }
    )
    
  } else {
    
    pearson_gamma_draws <-
      matrix(
        NA_real_,
        nrow = n_draws,
        ncol = N
      )
    
    for (i in seq_len(n_draws)) {
      
      pearson_gamma_draws[i, positive_idx] <-
        (y[positive_idx] -
           mu_draws[i, positive_idx]) /
        (
          mu_draws[i, positive_idx] /
            sqrt(shape_draws[i])
        )
    }
    
    pearson_gamma <- apply(
      pearson_gamma_draws,
      2,
      function(x) {
        
        if (all(!is.finite(x)))
          return(NA_real_)
        
        median(x[is.finite(x)])
      }
    )
  }
  
  
  # Observation-level diagnostic data frame -------------------------------
  
  diagnostic_data <- data.frame(
    observation = seq_len(N),
    y = y,
    fitted = fitted,
    fitted_gamma = fitted_gamma,
    raw_residual = raw_residual,
    rq_residual = rq_residual,
    observed_zero = y == 0,
    p_zero = p_zero,
    pearson_hurdle = pearson_hurdle,
    pearson_gamma = pearson_gamma
  )
  
  
  # Overall posterior predictive checks ---------------------------------------
  
  observed_stats <- c(
    calculate_stats(y),
    zero_proportion = mean(y == 0)
  )
  
  replicated_stats <- data.frame(
    mean = apply(y_rep, 1, mean),
    median = apply(y_rep, 1, median),
    variance = apply(y_rep, 1, var),
    sd = apply(y_rep, 1, sd),
    
    zero_proportion = apply(
      y_rep,
      1,
      function(x) mean(x == 0)
    ),
    
    q90 = apply(
      y_rep,
      1,
      quantile,
      probs = 0.90
    ),
    
    q95 = apply(
      y_rep,
      1,
      quantile,
      probs = 0.95
    ),
    
    q99 = apply(
      y_rep,
      1,
      quantile,
      probs = 0.99
    ),
    
    maximum = apply(y_rep, 1, max)
  )
  
  overall_ppc <- summarize_ppc(
    observed_stats,
    replicated_stats
  )
  
  
  # Positive-response posterior predictive checks ----------------------------
  
  positive_ppc <- NULL
  
  if (length(positive_idx) > 0) {
    
    observed_positive_stats <- calculate_stats(
      y[positive_idx]
    )
    
    positive_replicated_stats <- as.data.frame(
      t(
        apply(
          gamma_rep,
          1,
          calculate_stats
        )
      )
    )
    
    positive_ppc <- summarize_ppc(
      observed_positive_stats,
      positive_replicated_stats
    )
    
  }
  
  
  # Hurdle component diagnostics ----------------------------------------------
  
  observed_zero <- as.numeric(y == 0)
  
  # Brier score
  
  hurdle_brier <- mean(
    (observed_zero - p_zero)^2
  )
  
  # Log score
  
  p_zero_safe <- pmin(
    pmax(p_zero, 1e-10),
    1 - 1e-10
  )
  
  hurdle_log_score <- -mean(
    observed_zero * log(p_zero_safe) +
      (1 - observed_zero) * log(1 - p_zero_safe)
  )
  
  hurdle_summary <- data.frame(
    observed_zero_proportion = mean(observed_zero),
    predicted_zero_probability = mean(p_zero),
    brier_score = hurdle_brier,
    log_score = hurdle_log_score
  )
  
  
  # PPC plots -----------------------------------------------------------------
  
  n_plot <- min(100, nrow(y_rep))
  
  plot_idx <- sample(
    seq_len(nrow(y_rep)),
    n_plot
  )
  
  plot_draws <- y_rep[
    plot_idx,
    ,
    drop = FALSE
  ]
  
  
  # Overall density overlay
  
  ppc_density <- bayesplot::ppc_dens_overlay(
    y = y,
    yrep = plot_draws
  ) +
    ggplot2::labs(
      title = "Posterior predictive distribution",
      x = "Response",
      y = "Density"
    )
  
  
  # Overall ECDF overlay
  
  ppc_ecdf <- bayesplot::ppc_ecdf_overlay(
    y = y,
    yrep = plot_draws
  ) +
    ggplot2::labs(
      title = "Posterior predictive ECDF",
      x = "Response",
      y = "Cumulative probability"
    )
  
  
  # Overall mean PPC
  
  ppc_stat_plot <- bayesplot::ppc_stat(
    y = y,
    yrep = plot_draws,
    stat = "mean"
  ) +
    ggplot2::labs(
      title = "Posterior predictive check: mean",
      x = "Mean"
    )
  
  
  # Overall zero proportion PPC
  
  ppc_zero_plot <- bayesplot::ppc_stat(
    y = y,
    yrep = plot_draws,
    stat = function(x) mean(x == 0)
  ) +
    ggplot2::labs(
      title = "Posterior predictive check: zero proportion",
      x = "Proportion of zeros"
    )
  
  
  # Overall variance PPC
  
  ppc_variance_plot <- bayesplot::ppc_stat(
    y = y,
    yrep = plot_draws,
    stat = var
  ) +
    ggplot2::labs(
      title = "Posterior predictive check: variance",
      x = "Variance"
    )
  
  
  # Overall maximum PPC
  
  ppc_max_plot <- bayesplot::ppc_stat(
    y = y,
    yrep = plot_draws,
    stat = max
  ) +
    ggplot2::labs(
      title = "Posterior predictive check: maximum",
      x = "Maximum"
    )
  
  
  # Positive-response density and ECDF ----------------------------------------
  
  ppc_positive_density <- NULL
  ppc_positive_ecdf <- NULL
  ppc_positive_mean <- NULL
  
  if (length(positive_idx) > 0) {
    
    positive_plot_draws <- gamma_rep[
      plot_idx,
      ,
      drop = FALSE
    ]
    
    ppc_positive_density <-
      bayesplot::ppc_dens_overlay(
        y = y[positive_idx],
        yrep = positive_plot_draws
      ) +
      ggplot2::labs(
        title = "Conditional Gamma posterior predictive distribution",
        x = "Positive response",
        y = "Density"
      )
    
    ppc_positive_ecdf <-
      bayesplot::ppc_ecdf_overlay(
        y = y[positive_idx],
        yrep = positive_plot_draws
      ) +
      ggplot2::labs(
        title = "Conditional Gamma posterior predictive ECDF",
        x = "Positive response",
        y = "Cumulative probability"
      )
    
    ppc_positive_mean <- bayesplot::ppc_stat(
      y = y[positive_idx],
      yrep = positive_plot_draws,
      stat = "mean"
    ) +
      ggplot2::labs(
        title = "Conditional Gamma posterior predictive check: mean",
        x = "Positive-response mean"
      )
  }
  
  
  # Randomized quantile residual QQ plot --------------------------------------
  
  qq_data <- data.frame(
    rq_residual = rq_residual
  )
  
  qq_rqr <- ggplot2::ggplot(
    qq_data,
    ggplot2::aes(sample = rq_residual)
  ) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line() +
    ggplot2::labs(
      title = "Posterior predictive quantile residuals",
      x = "Theoretical quantiles",
      y = "Posterior predictive quantile residuals"
    ) +
    ggplot2::theme_classic()
  
  
  # RQR vs fitted -------------------------------------------------------------
  
  rqr_vs_fitted <- ggplot2::ggplot(
    diagnostic_data,
    ggplot2::aes(
      x = fitted,
      y = rq_residual
    )
  ) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      se = TRUE
    ) +
    ggplot2::labs(
      title = "Quantile residuals vs fitted values",
      x = "Unconditional fitted value",
      y = "Quantile residual"
    ) +
    ggplot2::theme_classic()
  
  
  # Hurdle Pearson residual vs predicted zero probability --------------------
  
  hurdle_residual_plot <- ggplot2::ggplot(
    diagnostic_data,
    ggplot2::aes(
      x = p_zero,
      y = pearson_hurdle
    )
  ) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      se = TRUE
    ) +
    ggplot2::labs(
      title = "Hurdle Pearson residuals vs predicted zero probability",
      x = "Predicted probability of zero",
      y = "Hurdle Pearson residual"
    ) +
    ggplot2::theme_classic()
  
  
  # Gamma Pearson residual vs conditional fitted values ----------------------
  
  gamma_residual_plot <- NULL
  
  if (length(positive_idx) > 0) {
    
    gamma_residual_plot <- ggplot2::ggplot(
      diagnostic_data[
        positive_idx,
        ,
        drop = FALSE
      ],
      ggplot2::aes(
        x = fitted_gamma,
        y = pearson_gamma
      )
    ) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed"
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        se = TRUE
      ) +
      ggplot2::labs(
        title = "Gamma Pearson residuals vs conditional fitted values",
        x = "Conditional Gamma mean",
        y = "Gamma Pearson residual"
      ) +
      ggplot2::theme_classic()
  }
  
  
  # Raw residual vs fitted ----------------------------------------------------
  
  raw_vs_fitted <- ggplot2::ggplot(
    diagnostic_data,
    ggplot2::aes(
      x = fitted,
      y = raw_residual
    )
  ) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      se = TRUE
    ) +
    ggplot2::labs(
      title = "Raw residuals vs fitted values",
      x = "Unconditional fitted value",
      y = "Raw residual"
    ) +
    ggplot2::theme_classic()
  
  
  # Residuals vs observation order --------------------------------------------
  
  residual_order <- ggplot2::ggplot(
    diagnostic_data,
    ggplot2::aes(
      x = observation,
      y = rq_residual
    )
  ) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      se = TRUE
    ) +
    ggplot2::labs(
      title = "Quantile residuals vs observation order",
      x = "Observation",
      y = "Quantile residual"
    ) +
    ggplot2::theme_classic()
  
  
  # Zero probability diagnostic -----------------------------------------------
  
  zero_plot <- ggplot2::ggplot(
    diagnostic_data,
    ggplot2::aes(
      x = p_zero,
      y = as.numeric(observed_zero)
    )
  ) +
    ggplot2::geom_jitter(
      height = 0.05,
      width = 0,
      alpha = 0.4
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      se = TRUE
    ) +
    ggplot2::labs(
      title = "Observed zeros vs predicted zero probability",
      x = "Predicted probability of zero",
      y = "Observed zero"
    ) +
    ggplot2::theme_classic()
  
  
  # Region × year posterior predictive diagnostics ----------------------------
  
  region_year <- interaction(
    rg[st],
    yr,
    drop = TRUE,
    sep = "_"
  )
  
  group_levels <- levels(region_year)
  
  group_diagnostics <- lapply(
    group_levels,
    function(g) {
      
      idx <- which(region_year == g)
      
      positive_group_idx <- which(y[idx] > 0)
      
      # Observed values
      
      y_obs <- y[idx]
      
      obs_mean <- mean(y_obs)
      obs_zero <- mean(y_obs == 0)
      
      # Fitted values
      
      fitted_group <- rowMeans(
        expected_draws[, idx, drop = FALSE]
      )
      
      # Posterior predictive group statistics
      
      yrep_group <- rowMeans(
        y_rep[, idx, drop = FALSE]
      )
      
      yrep_zero <- rowMeans(
        y_rep[, idx, drop = FALSE] == 0
      )
      
      expected_zero <- rowMeans(
        hu_draws[, idx, drop = FALSE]
      )
      
      # Positive-response component
      
      if (length(positive_group_idx) > 0) {
        
        fitted_gamma_group <- rowMeans(
          mu_draws[, idx, drop = FALSE]
        )
        
        positive_obs <- y_obs[positive_group_idx]
        
        gamma_columns <- idx[positive_group_idx]
        
        gamma_group <- rowMeans(
          gamma_rep[
            ,
            match(gamma_columns, positive_idx),
            drop = FALSE
          ]
        )
        
        obs_positive_mean <- mean(positive_obs)
        
        pred_positive_mean <- median(gamma_group)
        
        pred_positive_lower <- quantile(
          gamma_group,
          0.025
        )
        
        pred_positive_upper <- quantile(
          gamma_group,
          0.975
        )
        
        p_positive_mean <- mean(
          gamma_group >= obs_positive_mean
        )
        
      } else {
        
        obs_positive_mean <- NA_real_
        pred_positive_mean <- NA_real_
        pred_positive_lower <- NA_real_
        pred_positive_upper <- NA_real_
        p_positive_mean <- NA_real_
      }
      
      
      # Overall posterior predictive summaries
      
      pred_mean <- median(yrep_group)
      
      pred_mean_lower <- quantile(
        yrep_group,
        0.025
      )
      
      pred_mean_upper <- quantile(
        yrep_group,
        0.975
      )
      
      pred_zero <- median(yrep_zero)
      
      pred_zero_lower <- quantile(
        yrep_zero,
        0.025
      )
      
      pred_zero_upper <- quantile(
        yrep_zero,
        0.975
      )
      
      # Posterior predictive p-values
      
      p_mean <- mean(yrep_group >= obs_mean)
      
      p_zero <- mean(yrep_zero >= obs_zero)
      
      data.frame(
        group = g,
        n = length(idx),
        
        observed_mean = obs_mean,
        fitted_mean = median(fitted_group),
        
        predicted_mean = pred_mean,
        predicted_mean_lower = pred_mean_lower,
        predicted_mean_upper = pred_mean_upper,
        
        observed_zero = obs_zero,
        predicted_zero = pred_zero,
        predicted_zero_lower = pred_zero_lower,
        predicted_zero_upper = pred_zero_upper,
        
        expected_zero = median(expected_zero),
        
        observed_positive_mean = obs_positive_mean,
        fitted_positive_mean = median(fitted_gamma_group),
        
        predicted_positive_mean = pred_positive_mean,
        predicted_positive_lower = pred_positive_lower,
        predicted_positive_upper = pred_positive_upper,
        
        p_mean = p_mean,
        p_zero = p_zero,
        p_positive_mean = p_positive_mean
      )
    }
  )
  
  group_diagnostics <- do.call(
    rbind,
    group_diagnostics
  )
  
  
  # Observed vs predicted region × year means -------------------------------
  
  group_mean_plot <- ggplot2::ggplot(
    group_diagnostics,
    ggplot2::aes(
      x = predicted_mean,
      y = observed_mean
    )
  ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = predicted_mean_lower,
        xmax = predicted_mean_upper
      ),
      height = 0
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::labs(
      x = "Posterior predictive mean",
      y = "Observed mean",
      title = "Region × year posterior predictive check"
    ) +
    ggplot2::theme_classic()
  
  
  # Observed vs predicted zero proportions ------------------------------------
  
  group_zero_plot <- ggplot2::ggplot(
    group_diagnostics,
    ggplot2::aes(
      x = predicted_zero,
      y = observed_zero
    )
  ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = predicted_zero_lower,
        xmax = predicted_zero_upper
      ),
      height = 0
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::labs(
      x = "Posterior predictive zero proportion",
      y = "Observed zero proportion",
      title = "Region × year hurdle component check"
    ) +
    ggplot2::theme_classic()
  
  
  # Observed vs predicted positive means ---------------------------------------
  
  group_positive_plot <- NULL
  
  if (length(positive_idx) > 0) {
    
    group_positive_plot <- ggplot2::ggplot(
      group_diagnostics,
      ggplot2::aes(
        x = predicted_positive_mean,
        y = observed_positive_mean
      )
    ) +
      ggplot2::geom_errorbarh(
        ggplot2::aes(
          xmin = predicted_positive_lower,
          xmax = predicted_positive_upper
        ),
        height = 0
      ) +
      ggplot2::geom_point(
        na.rm = TRUE
      ) +
      ggplot2::geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed"
      ) +
      ggplot2::labs(
        x = "Posterior predictive positive mean",
        y = "Observed positive mean",
        title = "Region × year conditional Gamma check"
      ) +
      ggplot2::theme_classic()
  }
  
  
  # Return diagnostic object --------------------------------------------------
  
  diagnostics <- list(
    
    # Observation-level diagnostics
    
    data = diagnostic_data,
    
    # Posterior predictive simulations
    
    y_rep = y_rep,
    
    # Conditional Gamma posterior predictive simulations
    
    gamma_rep = gamma_rep,
    
    # Posterior predictive summary statistics
    
    ppc = list(
      
      observed = observed_stats,
      
      replicated = replicated_stats,
      
      summary = overall_ppc$summary
    ),
    
    # Conditional positive-response PPC
    
    positive_ppc = if (!is.null(positive_ppc)) {
      list(
        observed = observed_positive_stats,
        replicated = positive_replicated_stats,
        summary = positive_ppc$summary
      )
    } else {
      NULL
    },
    
    # Hurdle component diagnostics
    
    hurdle = list(
      summary = hurdle_summary,
      brier_score = hurdle_brier,
      log_score = hurdle_log_score
    ),
    
    # Region × year posterior predictive summaries
    
    group_ppc = list(
      summary = group_diagnostics
    ),
    
    # Bayesian p-values
    
    bayesian_p = overall_ppc$bayesian_p,
    
    bayesian_p_two_sided =
      overall_ppc$bayesian_p_two_sided,
    
    positive_bayesian_p = if (!is.null(positive_ppc)) {
      positive_ppc$bayesian_p
    } else {
      NULL
    },
    
    positive_bayesian_p_two_sided =
      if (!is.null(positive_ppc)) {
        positive_ppc$bayesian_p_two_sided
      } else {
        NULL
      },
    
    # Residuals
    
    residuals = list(
      raw = raw_residual,
      hurdle = pearson_hurdle,
      gamma = pearson_gamma,
      quantile = rq_residual
    ),
    
    # Predicted values
    
    fitted = fitted,
    fitted_gamma = fitted_gamma,
    
    # Predicted probability of zero
    
    p_zero = p_zero,
    
    # Posterior draws
    
    posterior = list(
      mu = mu_draws,
      hu = hu_draws,
      shape = shape_draws,
      expected = expected_draws
    ),
    
    # Log likelihood
    
    log_lik = log_lik,
    
    # Plots
    
    plots = list(
      
      rqr_qq = qq_rqr,
      rqr_vs_fitted = rqr_vs_fitted,
      
      hurdle_residuals = hurdle_residual_plot,
      gamma_residuals = gamma_residual_plot,
      
      raw_vs_fitted = raw_vs_fitted,
      residual_order = residual_order,
      
      zero_probability = zero_plot,
      
      zero_probability_group = group_zero_plot,
      
      group_mean = group_mean_plot,
      group_positive_mean = group_positive_plot,
      
      ppc_density = ppc_density,
      ppc_ecdf = ppc_ecdf,
      
      ppc_mean = ppc_stat_plot,
      ppc_zero = ppc_zero_plot,
      ppc_variance = ppc_variance_plot,
      ppc_max = ppc_max_plot,
      
      ppc_positive_density = ppc_positive_density,
      ppc_positive_ecdf = ppc_positive_ecdf,
      ppc_positive_mean = ppc_positive_mean
    ),
    
    # Settings
    
    settings = list(
      ndraws = nrow(y_rep),
      seed = seed
    )
  )
  
  class(diagnostics) <- "hurdle_diagnostics"
  
  return(diagnostics)
}

