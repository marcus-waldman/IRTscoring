# Stan Model Access and Compilation
#
# Functions for compiling and accessing the Stan model

#' @importFrom utils methods
NULL

#' Get Compiled Stan Model
#'
#' Returns a compiled Stan model object. Compiles the model on first call
#' and caches it for subsequent calls. Supports both rstan and cmdstanr backends.
#'
#' @param backend Character string specifying which backend to use: "cmdstanr"
#'   (default) or "rstan". If NULL, will try cmdstanr first, then rstan.
#' @param recompile Logical; if TRUE, force recompilation even if model is cached
#'
#' @return Compiled Stan model object (class depends on backend)
#'
#' @details
#' The function attempts to use cmdstanr by default as it typically provides
#' better performance and more features. If cmdstanr is not available, it falls
#' back to rstan.
#'
#' The compiled model is cached in the package environment to avoid
#' recompilation on subsequent calls.
#'
#' @examples
#' \dontrun{
#' # Get model using default backend (cmdstanr)
#' model <- get_stan_model()
#'
#' # Force use of rstan
#' model <- get_stan_model(backend = "rstan")
#'
#' # Force recompilation
#' model <- get_stan_model(recompile = TRUE)
#' }
#'
#' @export
get_stan_model <- function(backend = NULL, recompile = FALSE) {
  env <- .get_irtscoring_env()

  # Check if model is already compiled and cached
  if (!recompile && env$stan_model_compiled && !is.null(env$stan_model)) {
    return(env$stan_model)
  }

  # Determine which backend to use
  if (is.null(backend)) {
    # Auto-detect: try cmdstanr first, then rstan
    if (requireNamespace("cmdstanr", quietly = TRUE)) {
      backend <- "cmdstanr"
    } else if (requireNamespace("rstan", quietly = TRUE)) {
      backend <- "rstan"
    } else {
      stop("Neither cmdstanr nor rstan is available. Please install one of these packages.")
    }
  }

  # Validate backend choice
  backend <- match.arg(backend, choices = c("cmdstanr", "rstan"))

  # Get path to Stan model file
  model_path <- system.file("stan", "multidim_irt.stan", package = "IRTscoring")

  if (model_path == "" || !file.exists(model_path)) {
    stop("Stan model file not found. Expected at: inst/stan/multidim_irt.stan")
  }

  # Compile model based on backend
  message("Compiling Stan model using ", backend, "...")

  if (backend == "cmdstanr") {
    model <- compile_cmdstanr_model(model_path)
  } else if (backend == "rstan") {
    model <- compile_rstan_model(model_path)
  }

  # Cache the compiled model
  env$stan_model <- model
  env$stan_model_compiled <- TRUE
  env$backend <- backend

  message("Stan model compiled successfully")

  return(model)
}

#' Compile Stan Model with cmdstanr
#'
#' @param model_path Path to .stan file
#' @return Compiled cmdstanr model
#' @keywords internal
compile_cmdstanr_model <- function(model_path) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("cmdstanr package is required but not installed")
  }

  tryCatch({
    model <- cmdstanr::cmdstan_model(model_path)
    return(model)
  }, error = function(e) {
    stop("Failed to compile Stan model with cmdstanr: ", e$message)
  })
}

#' Compile Stan Model with rstan
#'
#' @param model_path Path to .stan file
#' @return Compiled rstan model
#' @keywords internal
compile_rstan_model <- function(model_path) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("rstan package is required but not installed")
  }

  tryCatch({
    model <- rstan::stan_model(file = model_path)
    return(model)
  }, error = function(e) {
    stop("Failed to compile Stan model with rstan: ", e$message)
  })
}

#' Check if Stan Model is Compiled
#'
#' @return Logical indicating whether the Stan model has been compiled
#' @export
is_stan_model_compiled <- function() {
  env <- .get_irtscoring_env()
  env$stan_model_compiled && !is.null(env$stan_model)
}

#' Get Current Stan Backend
#'
#' @return Character string: "cmdstanr", "rstan", or NULL if no model compiled
#' @export
get_stan_backend <- function() {
  env <- .get_irtscoring_env()
  if (env$stan_model_compiled) {
    return(env$backend)
  } else {
    return(NULL)
  }
}

# ============================================================================
# IRT Estimation Functions
# ============================================================================

#' IRT Model Fit Object
#'
#' S3 class for storing IRT model fit results
#'
#' @param theta Matrix of estimated ability parameters (N x K)
#' @param beta0 Vector of dimension intercepts (length K)
#' @param b Matrix of covariate effects (K x J)
#' @param tau Vector of dimension standard deviations (length K)
#' @param Omega Correlation matrix of dimensions (K x K)
#' @param log_lik Log-likelihood value
#' @param method Estimation method used
#' @param backend Stan backend used
#' @param stan_fit Original Stan fit object (optional)
#' @param ... Additional components
#'
#' @return Object of class "irt_fit"
#' @keywords internal
new_irt_fit <- function(theta, beta0, b, tau, Omega, log_lik,
                        method, backend, stan_fit = NULL, ...) {
  structure(
    list(
      theta = theta,
      beta0 = beta0,
      b = b,
      tau = tau,
      Omega = Omega,
      log_lik = log_lik,
      method = method,
      backend = backend,
      stan_fit = stan_fit,
      ...
    ),
    class = "irt_fit"
  )
}

#' Fit IRT Model using MAP Estimation
#'
#' Fits multidimensional IRT model using Maximum A Posteriori (MAP) estimation.
#' This is the fastest method, finding the mode of the posterior distribution.
#'
#' @param stan_data List of data prepared by \code{\link{prepare_irt_data}}
#' @param backend Stan backend to use: "cmdstanr" or "rstan"
#' @param ... Additional arguments passed to the optimizer
#'
#' @return Object of class "irt_fit" containing:
#'   \item{theta}{N x K matrix of estimated ability scores}
#'   \item{beta0}{Length K vector of dimension intercepts}
#'   \item{b}{K x J matrix of covariate effects}
#'   \item{tau}{Length K vector of dimension standard deviations}
#'   \item{Omega}{K x K correlation matrix}
#'   \item{log_lik}{Log-likelihood at MAP estimate}
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' stan_data <- prepare_irt_data(response_data, loadings, thresh_L, thresh_R)
#'
#' # Fit using MAP
#' fit <- fit_irt_map(stan_data)
#'
#' # Extract ability estimates
#' theta_hat <- fit$theta
#' }
#'
#' @export
fit_irt_map <- function(stan_data, backend = NULL, ...) {
  # Get compiled model
  model <- get_stan_model(backend = backend)
  current_backend <- get_stan_backend()

  # Run optimization based on backend
  if (current_backend == "cmdstanr") {
    fit <- model$optimize(data = stan_data, ...)

    # Extract MAP estimates
    draws <- fit$draws(format = "df")
    theta <- as.matrix(draws[, grepl("^theta\\[", names(draws))])
    theta <- matrix(theta, nrow = stan_data$N, ncol = stan_data$K)

    beta0 <- as.numeric(draws[, grepl("^beta0\\[", names(draws))])
    b <- as.matrix(draws[, grepl("^b\\[", names(draws))])
    b <- matrix(b, nrow = stan_data$K, ncol = stan_data$J)

    tau <- as.numeric(draws[, grepl("^tau\\[", names(draws))])
    Omega <- as.matrix(draws[, grepl("^Omega\\[", names(draws))])
    Omega <- matrix(Omega, nrow = stan_data$K, ncol = stan_data$K)

    log_lik <- as.numeric(draws[, "log_lik"])

  } else if (current_backend == "rstan") {
    fit <- rstan::optimizing(model, data = stan_data, ...)

    # Extract MAP estimates
    theta <- matrix(fit$par[grepl("^theta\\[", names(fit$par))],
                    nrow = stan_data$N, ncol = stan_data$K)

    beta0 <- fit$par[grepl("^beta0\\[", names(fit$par))]
    b <- matrix(fit$par[grepl("^b\\[", names(fit$par))],
                nrow = stan_data$K, ncol = stan_data$J)

    tau <- fit$par[grepl("^tau\\[", names(fit$par))]
    Omega <- matrix(fit$par[grepl("^Omega\\[", names(fit$par))],
                    nrow = stan_data$K, ncol = stan_data$K)

    log_lik <- fit$par["log_lik"]
  }

  # Return irt_fit object
  new_irt_fit(
    theta = theta,
    beta0 = beta0,
    b = b,
    tau = tau,
    Omega = Omega,
    log_lik = log_lik,
    method = "MAP",
    backend = current_backend,
    stan_fit = fit
  )
}

#' Fit IRT Model using Laplace Approximation
#'
#' Fits multidimensional IRT model using Laplace approximation. Finds the MAP
#' estimate then approximates the posterior with a multivariate normal.
#'
#' @param stan_data List of data prepared by \code{\link{prepare_irt_data}}
#' @param backend Stan backend to use: "cmdstanr" or "rstan"
#' @param draws Number of draws to sample from the Laplace approximation
#' @param ... Additional arguments passed to the optimizer
#'
#' @return Object of class "irt_fit" containing MAP estimates plus:
#'   \item{posterior_samples}{Matrix of draws from Laplace approximation}
#'   \item{hessian}{Hessian matrix at MAP estimate}
#'
#' @examples
#' \dontrun{
#' fit <- fit_irt_laplace(stan_data, draws = 1000)
#' }
#'
#' @export
fit_irt_laplace <- function(stan_data, backend = NULL, draws = 1000, ...) {
  # Get compiled model
  model <- get_stan_model(backend = backend)
  current_backend <- get_stan_backend()

  # Run optimization with Hessian
  if (current_backend == "cmdstanr") {
    # cmdstanr doesn't support Laplace directly, use optimization + manual sampling
    fit <- model$optimize(data = stan_data, jacobian = TRUE, ...)

    # Extract MAP and Hessian
    map_estimates <- fit$mle()

    # Parse parameter samples (this will sample from Laplace approximation)
    samples <- parse_parameter_samples(fit, draws)

    # Extract point estimates
    theta <- matrix(map_estimates[grepl("^theta\\[", names(map_estimates))],
                    nrow = stan_data$N, ncol = stan_data$K)
    beta0 <- map_estimates[grepl("^beta0\\[", names(map_estimates))]
    b <- matrix(map_estimates[grepl("^b\\[", names(map_estimates))],
                nrow = stan_data$K, ncol = stan_data$J)
    tau <- map_estimates[grepl("^tau\\[", names(map_estimates))]
    Omega <- matrix(map_estimates[grepl("^Omega\\[", names(map_estimates))],
                    nrow = stan_data$K, ncol = stan_data$K)
    log_lik <- map_estimates["log_lik"]

  } else if (current_backend == "rstan") {
    fit <- rstan::optimizing(model, data = stan_data, hessian = TRUE, ...)

    # Extract MAP estimates
    theta <- matrix(fit$par[grepl("^theta\\[", names(fit$par))],
                    nrow = stan_data$N, ncol = stan_data$K)
    beta0 <- fit$par[grepl("^beta0\\[", names(fit$par))]
    b <- matrix(fit$par[grepl("^b\\[", names(fit$par))],
                nrow = stan_data$K, ncol = stan_data$J)
    tau <- fit$par[grepl("^tau\\[", names(fit$par))]
    Omega <- matrix(fit$par[grepl("^Omega\\[", names(fit$par))],
                    nrow = stan_data$K, ncol = stan_data$K)
    log_lik <- fit$par["log_lik"]

    # Sample from Laplace approximation using Hessian
    hessian <- fit$hessian
    covariance <- solve(-hessian)
    samples <- MASS::mvrnorm(n = draws, mu = fit$par, Sigma = covariance)
  }

  # Return irt_fit object with samples
  new_irt_fit(
    theta = theta,
    beta0 = beta0,
    b = b,
    tau = tau,
    Omega = Omega,
    log_lik = log_lik,
    method = "Laplace",
    backend = current_backend,
    stan_fit = fit,
    posterior_samples = samples,
    hessian = if (current_backend == "rstan") fit$hessian else NULL
  )
}

#' Parse Parameter Samples from Optimization Fit
#'
#' Helper function to extract parameter samples for Laplace approximation
#'
#' @param fit Stan optimization fit object
#' @param draws Number of draws to sample
#' @return Matrix of parameter samples
#' @keywords internal
parse_parameter_samples <- function(fit, draws) {
  # This is a placeholder - actual implementation would compute Laplace samples
  # For cmdstanr, we'd need to manually compute the Hessian approximation
  warning("Full Laplace approximation not yet implemented for cmdstanr. ",
          "Use backend='rstan' for complete Laplace functionality.")
  return(NULL)
}

#' Fit IRT Model using Variational Bayes
#'
#' Fits multidimensional IRT model using Variational Bayes (VB) approximation.
#' Tries to use Pathfinder algorithm if available, otherwise falls back to ADVI.
#'
#' @param stan_data List of data prepared by \code{\link{prepare_irt_data}}
#' @param backend Stan backend to use: "cmdstanr" or "rstan"
#' @param algorithm VB algorithm: "pathfinder" (cmdstanr only) or "meanfield"
#' @param draws Number of draws to sample from variational distribution
#' @param ... Additional arguments passed to VB algorithm
#'
#' @return Object of class "irt_fit" containing variational estimates
#'
#' @examples
#' \dontrun{
#' # Use Pathfinder (fast, cmdstanr only)
#' fit <- fit_irt_vb(stan_data, algorithm = "pathfinder")
#'
#' # Use mean-field ADVI
#' fit <- fit_irt_vb(stan_data, algorithm = "meanfield")
#' }
#'
#' @export
fit_irt_vb <- function(stan_data, backend = NULL,
                       algorithm = c("pathfinder", "meanfield"),
                       draws = 1000, ...) {
  algorithm <- match.arg(algorithm)

  # Get compiled model
  model <- get_stan_model(backend = backend)
  current_backend <- get_stan_backend()

  # Run VB based on backend and algorithm
  if (current_backend == "cmdstanr") {
    if (algorithm == "pathfinder" && has_pathfinder()) {
      fit <- run_pathfinder(model, stan_data, draws, ...)
    } else {
      # Fall back to variational
      fit <- model$variational(data = stan_data, output_samples = draws, ...)
    }

    # Extract samples
    samples <- extract_vb_samples(fit, current_backend)

  } else if (current_backend == "rstan") {
    fit <- rstan::vb(model, data = stan_data,
                     algorithm = ifelse(algorithm == "pathfinder", "meanfield", algorithm),
                     output_samples = draws, ...)

    # Extract samples
    samples <- extract_vb_samples(fit, current_backend)
  }

  # Compute posterior means
  theta <- apply(samples$theta_samples, c(2, 3), mean)
  beta0 <- apply(samples$beta0_samples, 2, mean)
  b <- apply(samples$b_samples, c(2, 3), mean)
  tau <- apply(samples$tau_samples, 2, mean)
  Omega <- apply(samples$Omega_samples, c(2, 3), mean)
  log_lik <- mean(samples$log_lik_samples)

  # Return irt_fit object
  new_irt_fit(
    theta = theta,
    beta0 = beta0,
    b = b,
    tau = tau,
    Omega = Omega,
    log_lik = log_lik,
    method = paste0("VB_", algorithm),
    backend = current_backend,
    stan_fit = fit,
    posterior_samples = samples
  )
}

#' Fit IRT Model using MCMC
#'
#' Fits multidimensional IRT model using full Bayesian MCMC sampling.
#' This is the most computationally intensive but provides complete
#' posterior distributions.
#'
#' @param stan_data List of data prepared by \code{\link{prepare_irt_data}}
#' @param backend Stan backend to use: "cmdstanr" or "rstan"
#' @param chains Number of MCMC chains to run
#' @param iter Number of iterations per chain (rstan) or total (cmdstanr)
#' @param warmup Number of warmup iterations
#' @param ... Additional arguments passed to MCMC sampler
#'
#' @return Object of class "irt_fit" containing MCMC samples
#'
#' @examples
#' \dontrun{
#' # Run MCMC with 4 chains
#' fit <- fit_irt_mcmc(stan_data, chains = 4, iter = 2000)
#' }
#'
#' @export
fit_irt_mcmc <- function(stan_data, backend = NULL,
                         chains = 4, iter = 2000, warmup = 1000, ...) {
  # Get compiled model
  model <- get_stan_model(backend = backend)
  current_backend <- get_stan_backend()

  # Run MCMC based on backend
  if (current_backend == "cmdstanr") {
    fit <- model$sample(
      data = stan_data,
      chains = chains,
      iter_warmup = warmup,
      iter_sampling = iter - warmup,
      ...
    )

    # Extract posterior means
    draws <- fit$draws(format = "df")

    theta_cols <- grepl("^theta\\[", names(draws))
    theta_samples <- as.matrix(draws[, theta_cols])
    theta <- array(theta_samples,
                   dim = c(nrow(theta_samples), stan_data$N, stan_data$K))
    theta_mean <- apply(theta, c(2, 3), mean)

    beta0 <- apply(as.matrix(draws[, grepl("^beta0\\[", names(draws))]), 2, mean)

    b_cols <- grepl("^b\\[", names(draws))
    b_samples <- as.matrix(draws[, b_cols])
    b <- matrix(apply(b_samples, 2, mean), nrow = stan_data$K, ncol = stan_data$J)

    tau <- apply(as.matrix(draws[, grepl("^tau\\[", names(draws))]), 2, mean)

    Omega_cols <- grepl("^Omega\\[", names(draws))
    Omega_samples <- as.matrix(draws[, Omega_cols])
    Omega <- matrix(apply(Omega_samples, 2, mean),
                    nrow = stan_data$K, ncol = stan_data$K)

    log_lik <- mean(draws$log_lik)

  } else if (current_backend == "rstan") {
    fit <- rstan::sampling(
      model,
      data = stan_data,
      chains = chains,
      iter = iter,
      warmup = warmup,
      ...
    )

    # Extract posterior means
    theta_mean <- apply(rstan::extract(fit, "theta")[[1]], c(2, 3), mean)
    beta0 <- apply(rstan::extract(fit, "beta0")[[1]], 2, mean)
    b <- apply(rstan::extract(fit, "b")[[1]], c(2, 3), mean)
    tau <- apply(rstan::extract(fit, "tau")[[1]], 2, mean)
    Omega <- apply(rstan::extract(fit, "Omega")[[1]], c(2, 3), mean)
    log_lik <- mean(rstan::extract(fit, "log_lik")[[1]])
  }

  # Return irt_fit object
  new_irt_fit(
    theta = theta_mean,
    beta0 = beta0,
    b = b,
    tau = tau,
    Omega = Omega,
    log_lik = log_lik,
    method = "MCMC",
    backend = current_backend,
    stan_fit = fit
  )
}

# ============================================================================
# Helper Functions
# ============================================================================

#' Check if cmdstanr is Available
#'
#' @return Logical indicating if cmdstanr package is installed
#' @export
has_cmdstanr <- function() {
  requireNamespace("cmdstanr", quietly = TRUE)
}

#' Check if rstan is Available
#'
#' @return Logical indicating if rstan package is installed
#' @export
has_rstan <- function() {
  requireNamespace("rstan", quietly = TRUE)
}

#' Check if Pathfinder is Available
#'
#' @return Logical indicating if cmdstanr has Pathfinder support
#' @keywords internal
has_pathfinder <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    return(FALSE)
  }

  # Check if pathfinder method exists
  tryCatch({
    model <- get_stan_model(backend = "cmdstanr")
    return("pathfinder" %in% methods(class(model)))
  }, error = function(e) {
    return(FALSE)
  })
}

#' Run Pathfinder Algorithm
#'
#' @param model Compiled cmdstanr model
#' @param stan_data Data list
#' @param draws Number of draws
#' @param ... Additional arguments
#' @return Pathfinder fit object
#' @keywords internal
run_pathfinder <- function(model, stan_data, draws, ...) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("cmdstanr required for Pathfinder")
  }

  model$pathfinder(data = stan_data, draws = draws, ...)
}

#' Extract VB Samples
#'
#' Helper to extract samples from VB fit objects
#'
#' @param fit VB fit object
#' @param backend Backend used
#' @return List of parameter samples
#' @keywords internal
extract_vb_samples <- function(fit, backend) {
  if (backend == "cmdstanr") {
    draws <- fit$draws(format = "df")

    list(
      theta_samples = as.matrix(draws[, grepl("^theta\\[", names(draws))]),
      beta0_samples = as.matrix(draws[, grepl("^beta0\\[", names(draws))]),
      b_samples = as.matrix(draws[, grepl("^b\\[", names(draws))]),
      tau_samples = as.matrix(draws[, grepl("^tau\\[", names(draws))]),
      Omega_samples = as.matrix(draws[, grepl("^Omega\\[", names(draws))]),
      log_lik_samples = draws$log_lik
    )
  } else {
    list(
      theta_samples = rstan::extract(fit, "theta")[[1]],
      beta0_samples = rstan::extract(fit, "beta0")[[1]],
      b_samples = rstan::extract(fit, "b")[[1]],
      tau_samples = rstan::extract(fit, "tau")[[1]],
      Omega_samples = rstan::extract(fit, "Omega")[[1]],
      log_lik_samples = rstan::extract(fit, "log_lik")[[1]]
    )
  }
}
