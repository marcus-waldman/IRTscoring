# Tests for BridgeStan integration

# Helper to create simple test data
create_simple_bs_test_data <- function() {
  response_data <- data.frame(
    pid = c(1, 1, 2, 2, 3, 3),
    iid = c(1, 2, 1, 2, 1, 2)
  )

  loadings <- matrix(c(1, 1.5, 1, 1.5, 2, 1.2), ncol = 1)
  thresh_L <- c(-Inf, -1, -Inf, -1, -Inf, -2)
  thresh_R <- c(-1, Inf, -1, Inf, -2, Inf)

  prepare_irt_data(response_data, loadings, thresh_L, thresh_R)
}

# ============================================================================
# Tests for BridgeStan availability
# ============================================================================

test_that("has_bridgestan() returns logical", {
  result <- has_bridgestan()
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("has_bridgestan() detects package correctly", {
  # Should return TRUE if installed, FALSE otherwise
  bs_available <- requireNamespace("bridgestan", quietly = TRUE)
  expect_equal(has_bridgestan(), bs_available)
})

# ============================================================================
# Tests for JSON creation
# ============================================================================

test_that("create_stan_json() handles scalars", {
  data <- list(N = 10, K = 2)
  json <- create_stan_json(data)

  expect_type(json, "character")
  expect_true(grepl('"N": 10', json, fixed = TRUE))
  expect_true(grepl('"K": 2', json, fixed = TRUE))
})

test_that("create_stan_json() handles vectors", {
  data <- list(
    vec = c(1, 2, 3),
    inf = 999
  )
  json <- create_stan_json(data)

  expect_true(grepl('"vec": \\[1, 2, 3\\]', json))
  expect_true(grepl('"inf": 999', json, fixed = TRUE))
})

test_that("create_stan_json() handles matrices", {
  data <- list(
    A = matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  )
  json <- create_stan_json(data)

  # Should create array of arrays
  expect_true(grepl('"A":', json, fixed = TRUE))
  expect_true(grepl("\\[\\[", json))
})

test_that("create_stan_json() handles full stan_data", {
  stan_data <- create_simple_bs_test_data()
  json <- create_stan_json(stan_data)

  # Check all expected fields
  expect_true(grepl('"N":', json, fixed = TRUE))
  expect_true(grepl('"K":', json, fixed = TRUE))
  expect_true(grepl('"L":', json, fixed = TRUE))
  expect_true(grepl('"J":', json, fixed = TRUE))
  expect_true(grepl('"pid":', json, fixed = TRUE))
  expect_true(grepl('"dL":', json, fixed = TRUE))
  expect_true(grepl('"dR":', json, fixed = TRUE))
})

# ============================================================================
# Tests for BridgeStan model creation
# ============================================================================

test_that("create_bridgestan_model() fails gracefully without bridgestan", {
  skip_if(has_bridgestan(), "BridgeStan is installed")

  stan_data <- create_simple_bs_test_data()

  expect_error(
    create_bridgestan_model(stan_data),
    "bridgestan package is required"
  )
})

test_that("create_bridgestan_model() creates valid model", {
  skip_if_not_installed("bridgestan")

  stan_data <- create_simple_bs_test_data()

  expect_no_error({
    bs_model <- create_bridgestan_model(stan_data)
  })

  bs_model <- create_bridgestan_model(stan_data)
  expect_true(!is.null(bs_model))

  # Check that model has expected methods
  expect_true(is.function(bs_model$log_density))
  expect_true(is.function(bs_model$log_density_gradient))
  expect_true(is.function(bs_model$param_constrain))
})

test_that("create_bridgestan_model() validates model path", {
  skip_if_not_installed("bridgestan")

  stan_data <- create_simple_bs_test_data()

  expect_error(
    create_bridgestan_model(stan_data, model_path = "nonexistent.stan"),
    "not found"
  )
})

# ============================================================================
# Tests for helper functions
# ============================================================================

test_that("get_initial_values() returns correct length vector", {
  skip_if_not_installed("bridgestan")

  stan_data <- create_simple_bs_test_data()
  bs_model <- create_bridgestan_model(stan_data)

  init <- get_initial_values(stan_data, bs_model)

  expect_type(init, "double")
  expect_equal(length(init), bs_model$param_num())
})

test_that("get_initial_values() generates reasonable values", {
  skip_if_not_installed("bridgestan")

  stan_data <- create_simple_bs_test_data()
  bs_model <- create_bridgestan_model(stan_data)

  init <- get_initial_values(stan_data, bs_model)

  # Should be finite
  expect_true(all(is.finite(init)))

  # Tau parameters (SD) should be positive
  # These are typically at the end
  K <- stan_data$K
  tau_idx <- (length(init) - K + 1):length(init)
  expect_true(all(init[tau_idx] > 0))
})

test_that("objective_with_gradient() returns scalar", {
  skip_if_not_installed("bridgestan")

  stan_data <- create_simple_bs_test_data()
  bs_model <- create_bridgestan_model(stan_data)
  params <- get_initial_values(stan_data, bs_model)

  obj <- objective_with_gradient(params, bs_model, return_gradient = FALSE)

  expect_type(obj, "double")
  expect_length(obj, 1)
  expect_true(is.finite(obj))
})

test_that("objective_with_gradient() includes gradient when requested", {
  skip_if_not_installed("bridgestan")

  stan_data <- create_simple_bs_test_data()
  bs_model <- create_bridgestan_model(stan_data)
  params <- get_initial_values(stan_data, bs_model)

  obj <- objective_with_gradient(params, bs_model, return_gradient = TRUE)

  expect_true(!is.null(attr(obj, "gradient")))
  expect_equal(length(attr(obj, "gradient")), length(params))
  expect_true(all(is.finite(attr(obj, "gradient"))))
})

# ============================================================================
# Tests for fit_irt_optim()
# ============================================================================

test_that("fit_irt_optim() fails gracefully without bridgestan", {
  skip_if(has_bridgestan(), "BridgeStan is installed")

  stan_data <- create_simple_bs_test_data()

  expect_error(
    fit_irt_optim(stan_data),
    "bridgestan package is required"
  )
})

test_that("fit_irt_optim() runs with L-BFGS-B", {
  skip_if_not_installed("bridgestan")

  stan_data <- create_simple_bs_test_data()

  fit <- fit_irt_optim(
    stan_data,
    method = "L-BFGS-B",
    control = list(maxit = 50)
  )

  expect_s3_class(fit, "irt_fit")
  expect_equal(fit$method, "optim_L-BFGS-B")
  expect_equal(fit$backend, "bridgestan")

  # Check dimensions
  expect_equal(nrow(fit$theta), stan_data$N)
  expect_equal(ncol(fit$theta), stan_data$K)
  expect_equal(length(fit$beta0), stan_data$K)
  expect_equal(dim(fit$b), c(stan_data$K, stan_data$J))
})

test_that("fit_irt_optim() accepts custom initial values", {
  skip_if_not_installed("bridgestan")

  stan_data <- create_simple_bs_test_data()
  bs_model <- create_bridgestan_model(stan_data)
  custom_init <- get_initial_values(stan_data, bs_model) * 2

  expect_no_error({
    fit <- fit_irt_optim(stan_data, init = custom_init, control = list(maxit = 10))
  })
})

test_that("fit_irt_optim() works with BFGS method", {
  skip_if_not_installed("bridgestan")

  stan_data <- create_simple_bs_test_data()

  fit <- fit_irt_optim(
    stan_data,
    method = "BFGS",
    control = list(maxit = 50)
  )

  expect_s3_class(fit, "irt_fit")
  expect_equal(fit$method, "optim_BFGS")
})

test_that("fit_irt_optim() returns finite estimates", {
  skip_if_not_installed("bridgestan")

  stan_data <- create_simple_bs_test_data()

  fit <- fit_irt_optim(stan_data, control = list(maxit = 100))

  expect_true(all(is.finite(fit$theta)))
  expect_true(all(is.finite(fit$beta0)))
  expect_true(all(is.finite(fit$b)))
  expect_true(all(is.finite(fit$tau)))
  expect_true(all(fit$tau > 0))  # SDs should be positive
})

# ============================================================================
# Tests for fit_irt_nloptr()
# ============================================================================

test_that("fit_irt_nloptr() fails without nloptr", {
  skip_if_not_installed("bridgestan")
  skip_if(requireNamespace("nloptr", quietly = TRUE), "nloptr is installed")

  stan_data <- create_simple_bs_test_data()

  expect_error(
    fit_irt_nloptr(stan_data),
    "nloptr package is required"
  )
})

test_that("fit_irt_nloptr() runs with LBFGS", {
  skip_if_not_installed("bridgestan")
  skip_if_not_installed("nloptr")

  stan_data <- create_simple_bs_test_data()

  fit <- fit_irt_nloptr(
    stan_data,
    algorithm = "NLOPT_LD_LBFGS",
    opts = list(maxeval = 50)
  )

  expect_s3_class(fit, "irt_fit")
  expect_equal(fit$method, "nloptr_NLOPT_LD_LBFGS")
  expect_equal(fit$backend, "bridgestan")

  # Check structure
  expect_equal(nrow(fit$theta), stan_data$N)
  expect_equal(ncol(fit$theta), stan_data$K)
})

test_that("fit_irt_nloptr() accepts custom options", {
  skip_if_not_installed("bridgestan")
  skip_if_not_installed("nloptr")

  stan_data <- create_simple_bs_test_data()

  expect_no_error({
    fit <- fit_irt_nloptr(
      stan_data,
      opts = list(
        maxeval = 20,
        xtol_rel = 1e-4,
        ftol_rel = 1e-4
      )
    )
  })
})

test_that("fit_irt_nloptr() returns finite estimates", {
  skip_if_not_installed("bridgestan")
  skip_if_not_installed("nloptr")

  stan_data <- create_simple_bs_test_data()

  fit <- fit_irt_nloptr(stan_data, opts = list(maxeval = 100))

  expect_true(all(is.finite(fit$theta)))
  expect_true(all(is.finite(fit$tau)))
  expect_true(all(fit$tau > 0))
})

# ============================================================================
# Tests for parameter extraction
# ============================================================================

test_that("Parameter extraction functions work correctly", {
  skip_if_not_installed("bridgestan")

  stan_data <- create_simple_bs_test_data()
  bs_model <- create_bridgestan_model(stan_data)
  params <- get_initial_values(stan_data, bs_model)

  # Test each extraction function
  theta <- extract_theta_from_params(params, stan_data, bs_model)
  expect_equal(dim(theta), c(stan_data$N, stan_data$K))

  beta0 <- extract_beta0_from_params(params, stan_data, bs_model)
  expect_length(beta0, stan_data$K)

  b <- extract_b_from_params(params, stan_data, bs_model)
  expect_equal(dim(b), c(stan_data$K, stan_data$J))

  tau <- extract_tau_from_params(params, stan_data, bs_model)
  expect_length(tau, stan_data$K)

  Omega <- extract_omega_from_params(params, stan_data, bs_model)
  expect_equal(dim(Omega), c(stan_data$K, stan_data$K))
})

# ============================================================================
# Integration tests
# ============================================================================

test_that("BridgeStan and Stan MAP give similar results", {
  skip_if_not_installed("bridgestan")
  skip_if_not_installed("cmdstanr")

  stan_data <- create_simple_bs_test_data()

  # Fit with both methods
  fit_stan <- fit_irt_map(stan_data, backend = "cmdstanr")
  fit_bs <- fit_irt_optim(stan_data, control = list(maxit = 200))

  # Estimates should be reasonably close
  theta_diff <- mean(abs(fit_stan$theta - fit_bs$theta))
  expect_lt(theta_diff, 0.5)  # Within 0.5 units on average
})

test_that("optim() and nloptr give consistent results", {
  skip_if_not_installed("bridgestan")
  skip_if_not_installed("nloptr")

  stan_data <- create_simple_bs_test_data()

  fit_optim <- fit_irt_optim(stan_data, control = list(maxit = 200))
  fit_nloptr <- fit_irt_nloptr(stan_data, opts = list(maxeval = 200))

  # Should be similar
  theta_diff <- mean(abs(fit_optim$theta - fit_nloptr$theta))
  expect_lt(theta_diff, 0.3)
})
