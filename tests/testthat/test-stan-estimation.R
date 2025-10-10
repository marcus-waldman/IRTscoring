# Tests for Stan-based estimation functions

# Helper function to create simple test data
create_simple_test_data <- function() {
  # 3 persons, 4 responses, 1 dimension
  response_data <- data.frame(
    pid = c(1, 1, 2, 2, 3, 3, 3),
    iid = c(1, 2, 1, 2, 1, 2, 3)
  )

  loadings <- matrix(c(1, 1.5, 1, 1.5, 2, 1.2, 0.8), ncol = 1)
  thresh_L <- c(-Inf, -1, -Inf, -1, -Inf, -2, 0)
  thresh_R <- c(-1, Inf, -1, Inf, -2, 0, Inf)

  stan_data <- prepare_irt_data(response_data, loadings, thresh_L, thresh_R)
  return(stan_data)
}

# Helper for multidimensional test data
create_multidim_test_data <- function() {
  # 2 persons, 4 responses, 2 dimensions
  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )

  loadings <- matrix(c(1, 0.5, 1, 0.5,    # Dim 1
                       0.5, 1, 0.5, 1),    # Dim 2
                     ncol = 2)
  thresh_L <- c(-Inf, -1, -Inf, -1)
  thresh_R <- c(-1, Inf, -1, Inf)

  # Add covariates
  covariates <- matrix(c(1, 1,      # Intercept
                         25, 30),    # Age
                       ncol = 2)

  stan_data <- prepare_irt_data(response_data, loadings, thresh_L, thresh_R,
                                 person_covariates = covariates)
  return(stan_data)
}

# ============================================================================
# Tests for irt_fit S3 class
# ============================================================================

test_that("new_irt_fit creates valid irt_fit object", {
  theta <- matrix(c(0.5, -0.3), nrow = 2, ncol = 1)
  beta0 <- 0
  b <- matrix(c(0, 0.1), nrow = 1, ncol = 2)
  tau <- 1.2
  Omega <- matrix(1, nrow = 1, ncol = 1)
  log_lik <- -15.3

  fit <- new_irt_fit(
    theta = theta,
    beta0 = beta0,
    b = b,
    tau = tau,
    Omega = Omega,
    log_lik = log_lik,
    method = "MAP",
    backend = "cmdstanr"
  )

  expect_s3_class(fit, "irt_fit")
  expect_equal(fit$theta, theta)
  expect_equal(fit$beta0, beta0)
  expect_equal(fit$method, "MAP")
  expect_equal(fit$backend, "cmdstanr")
})

test_that("irt_fit can store additional components", {
  fit <- new_irt_fit(
    theta = matrix(0, 1, 1),
    beta0 = 0,
    b = matrix(0, 1, 1),
    tau = 1,
    Omega = matrix(1, 1, 1),
    log_lik = 0,
    method = "MAP",
    backend = "cmdstanr",
    custom_field = "custom_value"
  )

  expect_equal(fit$custom_field, "custom_value")
})

# ============================================================================
# Tests for MAP Estimation
# ============================================================================

test_that("fit_irt_map works with cmdstanr on simple data", {
  skip_if_not_installed("cmdstanr")

  stan_data <- create_simple_test_data()

  # Suppress Stan output
  expect_no_error({
    fit <- fit_irt_map(stan_data, backend = "cmdstanr")
  })

  expect_s3_class(fit, "irt_fit")
  expect_equal(fit$method, "MAP")
  expect_equal(fit$backend, "cmdstanr")

  # Check dimensions
  expect_equal(nrow(fit$theta), stan_data$N)
  expect_equal(ncol(fit$theta), stan_data$K)
  expect_equal(length(fit$beta0), stan_data$K)
  expect_equal(nrow(fit$b), stan_data$K)
  expect_equal(ncol(fit$b), stan_data$J)
  expect_equal(length(fit$tau), stan_data$K)
  expect_equal(dim(fit$Omega), c(stan_data$K, stan_data$K))
})

test_that("fit_irt_map works with rstan on simple data", {
  skip_if_not_installed("rstan")

  stan_data <- create_simple_test_data()

  expect_no_error({
    fit <- fit_irt_map(stan_data, backend = "rstan")
  })

  expect_s3_class(fit, "irt_fit")
  expect_equal(fit$method, "MAP")
  expect_equal(fit$backend, "rstan")

  # Check dimensions
  expect_equal(nrow(fit$theta), stan_data$N)
  expect_equal(ncol(fit$theta), stan_data$K)
})

test_that("fit_irt_map handles multidimensional model", {
  skip_if_not_installed("cmdstanr")

  stan_data <- create_multidim_test_data()

  expect_no_error({
    fit <- fit_irt_map(stan_data)
  })

  # Check multidimensional structure
  expect_equal(ncol(fit$theta), 2)  # K = 2
  expect_equal(length(fit$beta0), 2)
  expect_equal(dim(fit$b), c(2, 2))  # K x J
  expect_equal(dim(fit$Omega), c(2, 2))
})

test_that("MAP estimates are reasonable", {
  skip_if_not_installed("cmdstanr")

  stan_data <- create_simple_test_data()
  fit <- fit_irt_map(stan_data)

  # Theta should be numeric and finite
  expect_true(all(is.finite(fit$theta)))

  # Tau (SD) should be positive
  expect_true(all(fit$tau > 0))

  # Omega should be positive definite correlation matrix
  expect_true(all(diag(fit$Omega) == 1))  # Diagonal should be 1
  expect_true(all(fit$Omega >= -1 & fit$Omega <= 1))  # Correlations in [-1, 1]

  # Log-likelihood should be finite
  expect_true(is.finite(fit$log_lik))
})

# ============================================================================
# Tests for Laplace Approximation
# ============================================================================

test_that("fit_irt_laplace works with rstan", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("MASS")

  stan_data <- create_simple_test_data()

  expect_no_error({
    fit <- fit_irt_laplace(stan_data, backend = "rstan", draws = 100)
  })

  expect_s3_class(fit, "irt_fit")
  expect_equal(fit$method, "Laplace")

  # Should have posterior samples
  expect_true(!is.null(fit$posterior_samples))
  expect_equal(nrow(fit$posterior_samples), 100)

  # Should have Hessian
  expect_true(!is.null(fit$hessian))
})

test_that("fit_irt_laplace warns for cmdstanr", {
  skip_if_not_installed("cmdstanr")

  stan_data <- create_simple_test_data()

  # Should warn about limited functionality
  expect_warning(
    fit <- fit_irt_laplace(stan_data, backend = "cmdstanr", draws = 100),
    "not yet implemented for cmdstanr"
  )

  expect_s3_class(fit, "irt_fit")
})

# ============================================================================
# Tests for Variational Bayes
# ============================================================================

test_that("fit_irt_vb works with meanfield algorithm", {
  skip_if_not_installed("cmdstanr")

  stan_data <- create_simple_test_data()

  expect_no_error({
    fit <- fit_irt_vb(stan_data, algorithm = "meanfield", draws = 100)
  })

  expect_s3_class(fit, "irt_fit")
  expect_equal(fit$method, "VB_meanfield")

  # Check dimensions
  expect_equal(nrow(fit$theta), stan_data$N)
  expect_equal(ncol(fit$theta), stan_data$K)

  # Should have posterior samples
  expect_true(!is.null(fit$posterior_samples))
})

test_that("fit_irt_vb uses pathfinder when available", {
  skip_if_not_installed("cmdstanr")

  # Skip if pathfinder not available
  skip_if(!has_pathfinder(), "Pathfinder not available")

  stan_data <- create_simple_test_data()

  expect_no_error({
    fit <- fit_irt_vb(stan_data, algorithm = "pathfinder", draws = 100)
  })

  expect_equal(fit$method, "VB_pathfinder")
})

test_that("fit_irt_vb falls back gracefully when pathfinder unavailable", {
  skip_if_not_installed("cmdstanr")

  stan_data <- create_simple_test_data()

  # Request pathfinder - will fall back if not available
  expect_no_error({
    fit <- fit_irt_vb(stan_data, algorithm = "pathfinder", draws = 100)
  })

  expect_s3_class(fit, "irt_fit")
})

test_that("VB samples have correct structure", {
  skip_if_not_installed("cmdstanr")

  stan_data <- create_simple_test_data()
  fit <- fit_irt_vb(stan_data, draws = 100)

  samples <- fit$posterior_samples

  expect_true(!is.null(samples$theta_samples))
  expect_true(!is.null(samples$beta0_samples))
  expect_true(!is.null(samples$b_samples))
  expect_true(!is.null(samples$tau_samples))
  expect_true(!is.null(samples$Omega_samples))
  expect_true(!is.null(samples$log_lik_samples))
})

# ============================================================================
# Tests for MCMC
# ============================================================================

test_that("fit_irt_mcmc works with cmdstanr", {
  skip_if_not_installed("cmdstanr")

  stan_data <- create_simple_test_data()

  # Use very few iterations for speed
  expect_no_error({
    fit <- fit_irt_mcmc(stan_data, backend = "cmdstanr",
                        chains = 2, iter = 100, warmup = 50)
  })

  expect_s3_class(fit, "irt_fit")
  expect_equal(fit$method, "MCMC")
  expect_equal(fit$backend, "cmdstanr")

  # Check dimensions
  expect_equal(nrow(fit$theta), stan_data$N)
  expect_equal(ncol(fit$theta), stan_data$K)
})

test_that("fit_irt_mcmc works with rstan", {
  skip_if_not_installed("rstan")

  stan_data <- create_simple_test_data()

  # Use very few iterations for speed
  expect_no_error({
    fit <- fit_irt_mcmc(stan_data, backend = "rstan",
                        chains = 2, iter = 100, warmup = 50)
  })

  expect_s3_class(fit, "irt_fit")
  expect_equal(fit$method, "MCMC")
  expect_equal(fit$backend, "rstan")
})

test_that("MCMC handles multidimensional model", {
  skip_if_not_installed("cmdstanr")

  stan_data <- create_multidim_test_data()

  expect_no_error({
    fit <- fit_irt_mcmc(stan_data, chains = 2, iter = 100, warmup = 50)
  })

  # Check multidimensional structure
  expect_equal(ncol(fit$theta), 2)  # K = 2
  expect_equal(dim(fit$Omega), c(2, 2))
})

# ============================================================================
# Tests for helper functions
# ============================================================================

test_that("has_pathfinder detects availability correctly", {
  skip_if_not_installed("cmdstanr")

  result <- has_pathfinder()
  expect_type(result, "logical")
})

test_that("extract_vb_samples works with cmdstanr", {
  skip_if_not_installed("cmdstanr")

  stan_data <- create_simple_test_data()
  model <- get_stan_model(backend = "cmdstanr")
  fit <- model$variational(data = stan_data, output_samples = 50)

  samples <- extract_vb_samples(fit, "cmdstanr")

  expect_type(samples, "list")
  expect_true("theta_samples" %in% names(samples))
  expect_true("beta0_samples" %in% names(samples))
})

test_that("extract_vb_samples works with rstan", {
  skip_if_not_installed("rstan")

  stan_data <- create_simple_test_data()
  model <- get_stan_model(backend = "rstan")
  fit <- rstan::vb(model, data = stan_data, output_samples = 50)

  samples <- extract_vb_samples(fit, "rstan")

  expect_type(samples, "list")
  expect_true("theta_samples" %in% names(samples))
})

# ============================================================================
# Comparison tests
# ============================================================================

test_that("MAP and VB give similar results", {
  skip_if_not_installed("cmdstanr")

  stan_data <- create_simple_test_data()

  fit_map <- fit_irt_map(stan_data)
  fit_vb <- fit_irt_vb(stan_data, draws = 500)

  # Estimates should be reasonably close
  # (allowing for VB approximation error)
  theta_diff <- mean(abs(fit_map$theta - fit_vb$theta))
  expect_lt(theta_diff, 0.5)  # Within 0.5 units on average
})

test_that("Different backends give consistent results", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("rstan")

  stan_data <- create_simple_test_data()

  fit_cmdstanr <- fit_irt_map(stan_data, backend = "cmdstanr")
  fit_rstan <- fit_irt_map(stan_data, backend = "rstan")

  # MAP estimates should be very similar between backends
  theta_diff <- mean(abs(fit_cmdstanr$theta - fit_rstan$theta))
  expect_lt(theta_diff, 0.1)  # Very close agreement
})
