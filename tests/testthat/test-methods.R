# Tests for S3 Methods (print, summary, coef, predict)
#
# Comprehensive tests for IRT fit object methods

library(testthat)
library(IRTscoring)

# ============================================================================
# Helper Function to Create Test Fit Objects
# ============================================================================

create_simple_fit <- function() {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  response_data <- data.frame(
    pid = c(1, 1, 2, 2, 3, 3),
    iid = c(1, 2, 1, 2, 1, 2)
  )

  item_loadings <- matrix(c(1, 1.5, 1, 1.5, 2, 1.2), ncol = 1)
  threshold_left <- c(-Inf, -1, -Inf, -1, -Inf, -2)
  threshold_right <- c(-1, Inf, -1, Inf, -2, Inf)

  fit_ability(
    response_data = response_data,
    item_loadings = item_loadings,
    threshold_left = threshold_left,
    threshold_right = threshold_right,
    method = "map"
  )
}

create_multidim_fit <- function() {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  response_data <- data.frame(
    pid = c(1, 1, 2, 2, 3, 3),
    iid = c(1, 2, 1, 2, 1, 2)
  )

  item_loadings <- matrix(c(
    1, 0.5,
    1.5, 0.3,
    1, 0.5,
    1.5, 0.3,
    2, 0.8,
    1.2, 0.4
  ), ncol = 2, byrow = TRUE)

  threshold_left <- c(-Inf, -1, -Inf, -1, -Inf, -2)
  threshold_right <- c(-1, Inf, -1, Inf, -2, Inf)

  fit_ability(
    response_data = response_data,
    item_loadings = item_loadings,
    threshold_left = threshold_left,
    threshold_right = threshold_right,
    method = "map"
  )
}

# ============================================================================
# Tests for print.irt_fit()
# ============================================================================

test_that("print.irt_fit() produces output", {
  fit <- create_simple_fit()

  expect_output(print(fit), "IRT Model Fit")
  expect_output(print(fit), "Method:")
  expect_output(print(fit), "Backend:")
  expect_output(print(fit), "Dimensions:")
  expect_output(print(fit), "Persons:")
  expect_output(print(fit), "Log-likelihood:")
})

test_that("print.irt_fit() returns object invisibly", {
  fit <- create_simple_fit()

  result <- withVisible(print(fit))
  expect_false(result$visible)
  expect_s3_class(result$value, "irt_fit")
})

test_that("print.irt_fit() works with multidimensional fits", {
  fit <- create_multidim_fit()

  expect_output(print(fit), "Dimensions:.*2")
})

# ============================================================================
# Tests for summary.irt_fit()
# ============================================================================

test_that("summary.irt_fit() produces detailed output", {
  fit <- create_simple_fit()

  expect_output(summary(fit), "IRT Model Summary")
  expect_output(summary(fit), "Estimation Method:")
  expect_output(summary(fit), "Dimension Intercepts")
  expect_output(summary(fit), "Dimension Standard Deviations")
  expect_output(summary(fit), "Ability Estimates")
})

test_that("summary.irt_fit() shows correlations for multidimensional models", {
  fit <- create_multidim_fit()

  expect_output(summary(fit), "Dimension Correlations")
})

test_that("summary.irt_fit() returns object invisibly", {
  fit <- create_simple_fit()

  result <- withVisible(summary(fit))
  expect_false(result$visible)
  expect_s3_class(result$value, "irt_fit")
})

# ============================================================================
# Tests for coef.irt_fit()
# ============================================================================

test_that("coef.irt_fit() extracts beta0", {
  fit <- create_simple_fit()

  beta0 <- coef(fit, pars = "beta0")

  expect_type(beta0, "double")
  expect_length(beta0, 1)  # 1 dimension
})

test_that("coef.irt_fit() extracts b matrix", {
  fit <- create_simple_fit()

  b <- coef(fit, pars = "b")

  expect_true(is.matrix(b))
  expect_equal(nrow(b), 1)  # 1 dimension
})

test_that("coef.irt_fit() extracts tau", {
  fit <- create_simple_fit()

  tau <- coef(fit, pars = "tau")

  expect_type(tau, "double")
  expect_length(tau, 1)  # 1 dimension
  expect_true(all(tau > 0))  # Standard deviations must be positive
})

test_that("coef.irt_fit() extracts Omega", {
  fit <- create_multidim_fit()

  Omega <- coef(fit, pars = "Omega")

  expect_true(is.matrix(Omega))
  expect_equal(nrow(Omega), 2)
  expect_equal(ncol(Omega), 2)
  expect_equal(diag(Omega), c(1, 1))  # Diagonal should be 1
})

test_that("coef.irt_fit() extracts theta", {
  fit <- create_simple_fit()

  theta <- coef(fit, pars = "theta")

  expect_true(is.matrix(theta))
  expect_equal(nrow(theta), 3)  # 3 persons
  expect_equal(ncol(theta), 1)  # 1 dimension
})

test_that("coef.irt_fit() extracts multiple parameters", {
  fit <- create_simple_fit()

  params <- coef(fit, pars = c("beta0", "tau"))

  expect_type(params, "list")
  expect_equal(names(params), c("beta0", "tau"))
  expect_length(params$beta0, 1)
  expect_length(params$tau, 1)
})

test_that("coef.irt_fit() with default pars", {
  fit <- create_simple_fit()

  params <- coef(fit)

  expect_type(params, "list")
  expect_true(all(c("beta0", "b", "tau", "Omega") %in% names(params)))
})

# ============================================================================
# Tests for predict.irt_fit()
# ============================================================================

test_that("predict.irt_fit() returns abilities by default", {
  fit <- create_simple_fit()

  abilities <- predict(fit)

  expect_true(is.matrix(abilities))
  expect_equal(nrow(abilities), 3)  # 3 persons
  expect_equal(ncol(abilities), 1)  # 1 dimension
})

test_that("predict.irt_fit() with type='abilities'", {
  fit <- create_simple_fit()

  abilities <- predict(fit, type = "abilities")

  expect_true(is.matrix(abilities))
  expect_equal(dim(abilities), c(3, 1))
})

test_that("predict.irt_fit() with type='person_scores' (alias)", {
  fit <- create_simple_fit()

  scores1 <- predict(fit, type = "abilities")
  scores2 <- predict(fit, type = "person_scores")

  expect_equal(scores1, scores2)
})

test_that("predict.irt_fit() works with multidimensional models", {
  fit <- create_multidim_fit()

  abilities <- predict(fit)

  expect_equal(nrow(abilities), 3)  # 3 persons
  expect_equal(ncol(abilities), 2)  # 2 dimensions
})

# ============================================================================
# Tests for Utility Functions
# ============================================================================

test_that("get_dimension_correlations() returns correlation matrix", {
  fit <- create_multidim_fit()

  Omega <- get_dimension_correlations(fit)

  expect_true(is.matrix(Omega))
  expect_equal(dim(Omega), c(2, 2))
  expect_equal(diag(Omega), c(1, 1))
  expect_true(all(Omega >= -1 & Omega <= 1))
})

test_that("get_dimension_correlations() returns 1x1 for unidimensional", {
  fit <- create_simple_fit()

  Omega <- get_dimension_correlations(fit)

  expect_equal(Omega, matrix(1, 1, 1))
})

test_that("get_individual_effects() extracts deviations", {
  fit <- create_simple_fit()

  ind_effects <- get_individual_effects(fit)

  expect_true(is.matrix(ind_effects))
  expect_equal(nrow(ind_effects), 3)  # 3 persons
  expect_equal(ncol(ind_effects), 1)  # 1 dimension

  # Individual effects should be centered (roughly, given small sample)
  # For intercept-only model: ind_effects = theta - beta0
  beta0 <- coef(fit, pars = "beta0")
  theta <- coef(fit, pars = "theta")
  expected <- theta - beta0
  expect_equal(ind_effects, expected)
})

test_that("get_covariate_effects() returns intercept effects", {
  fit <- create_simple_fit()

  cov_effects <- get_covariate_effects(fit)

  expect_true(is.matrix(cov_effects))
  expect_equal(nrow(cov_effects), 3)  # 3 persons
  expect_equal(ncol(cov_effects), 1)  # 1 dimension

  # For intercept-only model, all persons should have same covariate effect
  beta0 <- coef(fit, pars = "beta0")
  expect_equal(cov_effects[, 1], rep(beta0, 3))
})

# ============================================================================
# Tests for Score Decomposition
# ============================================================================

test_that("abilities equal sum of covariate and individual effects", {
  fit <- create_simple_fit()

  abilities <- predict(fit, type = "abilities")
  cov_effects <- get_covariate_effects(fit)
  ind_effects <- get_individual_effects(fit)

  # theta = beta0 + individual_effects
  # Which means: theta = cov_effects + ind_effects (for intercept-only)
  reconstructed <- cov_effects + ind_effects

  expect_equal(abilities, reconstructed, tolerance = 1e-10)
})

test_that("abilities decomposition works for multidimensional models", {
  fit <- create_multidim_fit()

  abilities <- predict(fit)
  cov_effects <- get_covariate_effects(fit)
  ind_effects <- get_individual_effects(fit)

  reconstructed <- cov_effects + ind_effects

  expect_equal(abilities, reconstructed, tolerance = 1e-10)
})

# ============================================================================
# Tests for Different Backends
# ============================================================================

test_that("print method works with cmdstanr backend", {
  skip_if_not(has_cmdstanr(), "cmdstanr not available")

  response_data <- data.frame(pid = c(1, 1), iid = c(1, 2))
  loadings <- matrix(c(1, 1.5), ncol = 1)

  fit <- fit_ability(
    response_data = response_data,
    item_loadings = loadings,
    threshold_left = c(-Inf, -1),
    threshold_right = c(-1, Inf),
    method = "map",
    backend = "cmdstanr"
  )

  expect_output(print(fit), "Backend:.*cmdstanr")
})

test_that("coef works regardless of backend", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend")

  fit <- create_simple_fit()

  params <- coef(fit)

  expect_type(params, "list")
  expect_true("beta0" %in% names(params))
})

# ============================================================================
# Tests for Edge Cases
# ============================================================================

test_that("methods work with single person", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  response_data <- data.frame(pid = c(1, 1), iid = c(1, 2))
  loadings <- matrix(c(1, 1.5), ncol = 1)

  fit <- fit_ability(
    response_data = response_data,
    item_loadings = loadings,
    threshold_left = c(-Inf, -1),
    threshold_right = c(-1, Inf),
    method = "map"
  )

  expect_output(print(fit), "Persons:.*1")
  expect_output(summary(fit), "Number of Persons:.*1")

  theta <- predict(fit)
  expect_equal(nrow(theta), 1)
})

test_that("methods handle parameter names correctly", {
  fit <- create_simple_fit()

  # Test that we can extract parameters by different combinations
  expect_no_error(coef(fit, pars = "beta0"))
  expect_no_error(coef(fit, pars = c("beta0", "tau")))
  expect_no_error(coef(fit, pars = c("beta0", "b", "tau", "Omega")))
})
