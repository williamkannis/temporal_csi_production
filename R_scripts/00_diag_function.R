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
    
    theta_cols <- grep(
      "^theta$",
      colnames(draws)
    )
    
    if (length(phi_cols) == 1 &&
        length(theta_cols) == 1 &&
        !is.null(raw_resid)) {
      
      phi_draw <- draws[, phi_cols]
      theta_draw <- draws[, theta_cols]
      
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
              mu_i^theta_draw[i]
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
    theta_cols <- grep("^theta$", colnames(draws))
    
    if (length(phi_cols) == 1 &&
        length(theta_cols) == 1 &&
        !is.null(raw_resid)) {
      
      phi_draw <- draws[, phi_cols]
      theta_draw <- draws[, theta_cols]
      
      p_zero_draws <- matrix(
        NA_real_,
        nrow = nrow(draws),
        ncol = N
      )
      
      for (i in seq_len(nrow(draws))) {
        
        mu_i <- y - raw_resid[i, ]
        
        lambda_i <-
          mu_i^(2 - theta_draw[i]) /
          (
            phi_draw[i] *
              (2 - theta_draw[i])
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
