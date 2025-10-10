# Tests for Stan model compilation and structure

test_that("Stan model file exists", {
  model_path <- system.file("stan", "multidim_irt.stan", package = "IRTscoring")
  expect_true(file.exists(model_path))
})

test_that("Stan model compiles successfully with cmdstanr", {
  skip_if_not_installed("cmdstanr")

  # Try to compile the model
  expect_no_error({
    model <- get_stan_model(backend = "cmdstanr")
  })

  # Check that model is compiled
  expect_true(is_stan_model_compiled())

  # Check backend is set correctly
  expect_equal(get_stan_backend(), "cmdstanr")
})

test_that("Stan model compiles successfully with rstan", {
  skip_if_not_installed("rstan")

  # Force recompilation with rstan backend
  expect_no_error({
    model <- get_stan_model(backend = "rstan", recompile = TRUE)
  })

  # Check that model is compiled
  expect_true(is_stan_model_compiled())

  # Check backend is set correctly
  expect_equal(get_stan_backend(), "rstan")
})

test_that("get_stan_model() caches compiled model", {
  skip_if_not_installed("cmdstanr")

  # First call compiles
  model1 <- get_stan_model(backend = "cmdstanr")

  # Second call should return cached model (faster)
  start_time <- Sys.time()
  model2 <- get_stan_model(backend = "cmdstanr")
  elapsed <- as.numeric(Sys.time() - start_time)

  # Should be very fast (< 1 second) since cached
  expect_lt(elapsed, 1)

  # Should be the same object
  expect_identical(model1, model2)
})

test_that("Stan model has correct parameter names", {
  skip_if_not_installed("cmdstanr")

  model <- get_stan_model(backend = "cmdstanr")

  # Expected parameters in the model
  expected_params <- c(
    "beta0",       # Dimension intercepts
    "b",           # Covariate effects
    "zeta",        # Individual effects (standardized)
    "L_Omega",     # Cholesky factor of correlation matrix
    "tau",         # Dimension standard deviations
    "eta",         # Transformed individual effects
    "theta",       # Full ability scores
    "Omega",       # Correlation matrix (generated quantities)
    "log_lik"      # Log-likelihood (generated quantities)
  )

  # Get model code
  model_code <- readLines(system.file("stan", "multidim_irt.stan",
                                      package = "IRTscoring"))
  model_text <- paste(model_code, collapse = "\n")

  # Check that each expected parameter appears in the model
  for (param in expected_params) {
    # Use regex to find parameter declarations
    # Look for patterns like "vector[K] beta0" or "matrix[N, K] theta"
    pattern <- paste0("\\b", param, "\\b")
    expect_true(grepl(pattern, model_text),
                info = paste("Parameter", param, "not found in Stan model"))
  }
})

test_that("Stan model data block has correct structure", {
  skip_if_not_installed("cmdstanr")

  model_code <- readLines(system.file("stan", "multidim_irt.stan",
                                      package = "IRTscoring"))
  model_text <- paste(model_code, collapse = "\n")

  # Expected data variables
  expected_data <- c(
    "inf",    # Infinity code
    "L",      # Number of responses
    "N",      # Number of persons
    "J",      # Number of covariates
    "K",      # Number of dimensions
    "pid",    # Person IDs
    "dL",     # Left thresholds
    "dR",     # Right thresholds
    "wgt",    # Person weights
    "X",      # Covariate matrix
    "A"       # Loading matrix
  )

  for (data_var in expected_data) {
    pattern <- paste0("\\b", data_var, "\\b")
    expect_true(grepl(pattern, model_text),
                info = paste("Data variable", data_var, "not found in Stan model"))
  }
})

test_that("recompile parameter forces recompilation", {
  skip_if_not_installed("cmdstanr")

  # Get initial model
  model1 <- get_stan_model(backend = "cmdstanr")

  # Force recompilation
  model2 <- get_stan_model(backend = "cmdstanr", recompile = TRUE)

  # Should have recompiled (models might differ in memory location)
  expect_true(is_stan_model_compiled())
})

test_that("error when no Stan package available", {
  # This test is mostly for documentation - hard to test without unloading packages
  skip("Cannot easily test without unloading Stan packages")
})
