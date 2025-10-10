# Tests for Main User Interface Functions
#
# Comprehensive tests for fit_ability(), fit_unidim_ability(), and fscores()

library(testthat)
library(IRTscoring)

# ============================================================================
# Test Data Setup
# ============================================================================

# Helper function to create basic test data
create_test_data <- function() {
  list(
    response_data = data.frame(
      pid = c(1, 1, 2, 2, 3, 3),
      iid = c(1, 2, 1, 2, 1, 2)
    ),
    item_loadings = matrix(c(1, 1.5, 1, 1.5, 2, 1.2), ncol = 1),
    threshold_left = c(-Inf, -1, -Inf, -1, -Inf, -2),
    threshold_right = c(-1, Inf, -1, Inf, -2, Inf)
  )
}

# Helper function to create multidimensional test data
create_multidim_data <- function() {
  list(
    response_data = data.frame(
      pid = c(1, 1, 2, 2, 3, 3),
      iid = c(1, 2, 1, 2, 1, 2)
    ),
    item_loadings = matrix(c(
      1, 0.5,
      1.5, 0.3,
      1, 0.5,
      1.5, 0.3,
      2, 0.8,
      1.2, 0.4
    ), ncol = 2, byrow = TRUE),
    threshold_left = c(-Inf, -1, -Inf, -1, -Inf, -2),
    threshold_right = c(-1, Inf, -1, Inf, -2, Inf)
  )
}

# ============================================================================
# Tests for fit_ability()
# ============================================================================

test_that("fit_ability() works with MAP estimation", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  data <- create_test_data()

  fit <- fit_ability(
    response_data = data$response_data,
    item_loadings = data$item_loadings,
    threshold_left = data$threshold_left,
    threshold_right = data$threshold_right,
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
  expect_equal(fit$method, "MAP")
  expect_true(nrow(fit$theta) == 3)  # 3 persons
  expect_true(ncol(fit$theta) == 1)  # 1 dimension
  expect_true(length(fit$beta0) == 1)
  expect_true(is.numeric(fit$log_lik))
})

test_that("fit_ability() works with multidimensional models", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  data <- create_multidim_data()

  fit <- fit_ability(
    response_data = data$response_data,
    item_loadings = data$item_loadings,
    threshold_left = data$threshold_left,
    threshold_right = data$threshold_right,
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
  expect_true(nrow(fit$theta) == 3)  # 3 persons
  expect_true(ncol(fit$theta) == 2)  # 2 dimensions
  expect_true(length(fit$beta0) == 2)
  expect_true(nrow(fit$Omega) == 2 && ncol(fit$Omega) == 2)
})

test_that("fit_ability() accepts person covariates", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  data <- create_test_data()

  # Add covariates (with intercept)
  covariates <- matrix(c(
    1, 25,
    1, 30,
    1, 35
  ), ncol = 2, byrow = TRUE)

  fit <- fit_ability(
    response_data = data$response_data,
    item_loadings = data$item_loadings,
    threshold_left = data$threshold_left,
    threshold_right = data$threshold_right,
    person_covariates = covariates,
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
  expect_true(ncol(fit$b) == 2)  # 2 covariates (intercept + age)
})

test_that("fit_ability() accepts weights", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  data <- create_test_data()
  weights <- c(1, 2, 1)

  fit <- fit_ability(
    response_data = data$response_data,
    item_loadings = data$item_loadings,
    threshold_left = data$threshold_left,
    threshold_right = data$threshold_right,
    weights = weights,
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
})

test_that("fit_ability() routes to correct estimation method", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  data <- create_test_data()

  # Test MAP routing
  fit_map <- fit_ability(
    response_data = data$response_data,
    item_loadings = data$item_loadings,
    threshold_left = data$threshold_left,
    threshold_right = data$threshold_right,
    method = "map"
  )
  expect_equal(fit_map$method, "MAP")

  # Test Laplace routing (if rstan available)
  if (has_rstan()) {
    fit_laplace <- fit_ability(
      response_data = data$response_data,
      item_loadings = data$item_loadings,
      threshold_left = data$threshold_left,
      threshold_right = data$threshold_right,
      method = "laplace"
    )
    expect_equal(fit_laplace$method, "Laplace")
  }
})

test_that("fit_ability() validates method argument", {
  data <- create_test_data()

  expect_error(
    fit_ability(
      response_data = data$response_data,
      item_loadings = data$item_loadings,
      threshold_left = data$threshold_left,
      threshold_right = data$threshold_right,
      method = "invalid_method"
    ),
    "should be one of"
  )
})

# ============================================================================
# Tests for fit_unidim_ability()
# ============================================================================

test_that("fit_unidim_ability() converts discrimination vector to loadings matrix", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )

  fit <- fit_unidim_ability(
    response_data = response_data,
    item_discrimination = c(1, 1.5, 1, 1.5),
    threshold_left = c(-Inf, -1, -Inf, -1),
    threshold_right = c(-1, Inf, -1, Inf),
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
  expect_true(ncol(fit$theta) == 1)  # Unidimensional
  expect_equal(fit$method, "MAP")
})

test_that("fit_unidim_ability() passes additional arguments", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )

  weights <- c(1, 2)

  fit <- fit_unidim_ability(
    response_data = response_data,
    item_discrimination = c(1, 1.5, 1, 1.5),
    threshold_left = c(-Inf, -1, -Inf, -1),
    threshold_right = c(-1, Inf, -1, Inf),
    weights = weights,
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
})

# ============================================================================
# Tests for fscores() - Core Functionality
# ============================================================================

test_that("fscores() works with wide format data (0-indexed responses)", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  # Create wide format data with 0-indexed responses
  wide_data <- data.frame(
    person_id = 1:3,
    item1 = c(0, 0, 1),  # 0-indexed!
    item2 = c(1, 1, 0)   # 0-indexed!
  )

  # Define item parameters
  item_params <- list(
    item1 = list(
      loadings = c(1.2),
      thresholds = c(0.5)  # 1 threshold for 2 categories (0, 1)
    ),
    item2 = list(
      loadings = c(0.8),
      thresholds = c(-0.5)
    )
  )

  fit <- fscores(
    data = wide_data,
    item_params = item_params,
    person_id_col = "person_id",
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
  expect_true(nrow(fit$theta) == 3)  # 3 persons
  expect_true(ncol(fit$theta) == 1)  # 1 dimension
})

test_that("fscores() works with multidimensional items", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  wide_data <- data.frame(
    person_id = 1:3,
    item1 = c(0, 1, 0),
    item2 = c(1, 0, 1)
  )

  # 2-dimensional items
  item_params <- list(
    item1 = list(
      loadings = c(1.2, 0.1),
      thresholds = c(0.5)
    ),
    item2 = list(
      loadings = c(0.8, 1.5),
      thresholds = c(-0.5)
    )
  )

  fit <- fscores(
    data = wide_data,
    item_params = item_params,
    person_id_col = "person_id",
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
  expect_true(ncol(fit$theta) == 2)  # 2 dimensions
})

test_that("fscores() works with 3-category items (0-indexed)", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  wide_data <- data.frame(
    person_id = 1:4,
    item1 = c(0, 1, 2, 1)  # 0-indexed: categories 0, 1, 2
  )

  item_params <- list(
    item1 = list(
      loadings = c(1.0),
      thresholds = c(-0.5, 0.5)  # 2 thresholds for 3 categories
    )
  )

  fit <- fscores(
    data = wide_data,
    item_params = item_params,
    person_id_col = "person_id",
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
  expect_true(nrow(fit$theta) == 4)
})

test_that("fscores() handles missing responses (NA)", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  wide_data <- data.frame(
    person_id = 1:3,
    item1 = c(0, NA, 1),
    item2 = c(1, 1, NA)
  )

  item_params <- list(
    item1 = list(loadings = c(1.2), thresholds = c(0.5)),
    item2 = list(loadings = c(0.8), thresholds = c(-0.5))
  )

  fit <- fscores(
    data = wide_data,
    item_params = item_params,
    person_id_col = "person_id",
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
  expect_true(nrow(fit$theta) == 3)  # All 3 persons should be included
})

test_that("fscores() works with covariates", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  wide_data <- data.frame(
    person_id = 1:3,
    age = c(25, 30, 35),
    item1 = c(0, 0, 1),
    item2 = c(1, 1, 0)
  )

  item_params <- list(
    item1 = list(loadings = c(1.2), thresholds = c(0.5)),
    item2 = list(loadings = c(0.8), thresholds = c(-0.5))
  )

  fit <- fscores(
    data = wide_data,
    item_params = item_params,
    person_id_col = "person_id",
    covariate_cols = "age",
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
  # Should have intercept + age covariate
  expect_true(ncol(fit$b) == 2)
})

# ============================================================================
# Tests for fscores() - Validation (0-indexed requirement)
# ============================================================================

test_that("fscores() rejects 1-indexed responses", {
  # Data with 1-indexed responses (incorrect)
  # Using categories 1, 2, 3 instead of 0, 1, 2
  wide_data <- data.frame(
    person_id = 1:3,
    item1 = c(1, 2, 3)  # WRONG: should be 0, 1, 2 for 3 categories
  )

  item_params <- list(
    item1 = list(
      loadings = c(1.0),
      thresholds = c(-0.5, 0.5)  # 2 thresholds for what should be 3 categories
    )
  )

  # This should fail because response=3 (0-indexed category 3) needs at least 3 thresholds
  # The error message will mention needing thresholds, which is the symptom of using 1-indexed data
  expect_error(
    fscores(
      data = wide_data,
      item_params = item_params,
      person_id_col = "person_id",
      method = "map"
    ),
    "requires at least 3 threshold"
  )
})

test_that("fscores() explicitly detects when user likely provided 1-indexed data", {
  # This test documents the actual error users will see
  # When they use 1-indexed responses (1, 2) instead of 0-indexed (0, 1)
  wide_data <- data.frame(
    person_id = 1:2,
    item1 = c(1, 2)  # If meant as binary, should be c(0, 1)
  )

  item_params <- list(
    item1 = list(
      loadings = c(1.0),
      thresholds = c(0.5)  # 1 threshold for 2 categories (0, 1)
    )
  )

  # Will fail because response=2 needs at least 2 thresholds
  expect_error(
    fscores(
      data = wide_data,
      item_params = item_params,
      person_id_col = "person_id",
      method = "map"
    ),
    "requires at least 2 threshold"
  )
})

test_that("fscores() rejects negative response values", {
  wide_data <- data.frame(
    person_id = 1:3,
    item1 = c(0, -1, 1)  # Invalid: negative value
  )

  item_params <- list(
    item1 = list(loadings = c(1.0), thresholds = c(0.5))
  )

  expect_error(
    fscores(
      data = wide_data,
      item_params = item_params,
      person_id_col = "person_id",
      method = "map"
    ),
    "Response values must be >= 0 \\(0-indexed\\)"
  )
})

# ============================================================================
# Tests for fscores() - Threshold Validation
# ============================================================================

test_that("fscores() validates sufficient thresholds for observed responses", {
  wide_data <- data.frame(
    person_id = 1:2,
    item1 = c(0, 2)  # Response 2 requires at least 2 thresholds
  )

  # Only 1 threshold provided (insufficient for category 2)
  item_params <- list(
    item1 = list(
      loadings = c(1.0),
      thresholds = c(0.5)  # Only 1 threshold, but we have response=2
    )
  )

  expect_error(
    fscores(
      data = wide_data,
      item_params = item_params,
      person_id_col = "person_id",
      method = "map"
    ),
    "has observed response 2.*requires at least 2 threshold"
  )
})

test_that("fscores() accepts sufficient thresholds for all observed responses", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  wide_data <- data.frame(
    person_id = 1:3,
    item1 = c(0, 1, 2)  # Categories 0, 1, 2
  )

  # 2 thresholds for 3 categories (correct)
  item_params <- list(
    item1 = list(
      loadings = c(1.0),
      thresholds = c(-0.5, 0.5)  # 2 thresholds
    )
  )

  # Should not error
  fit <- fscores(
    data = wide_data,
    item_params = item_params,
    person_id_col = "person_id",
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
})

# ============================================================================
# Tests for fscores() - Input Validation
# ============================================================================

test_that("fscores() validates person_id_col exists", {
  wide_data <- data.frame(
    pid = 1:3,  # Named 'pid', not 'person_id'
    item1 = c(0, 1, 0)
  )

  item_params <- list(
    item1 = list(loadings = c(1.0), thresholds = c(0.5))
  )

  expect_error(
    fscores(
      data = wide_data,
      item_params = item_params,
      person_id_col = "person_id",  # Doesn't exist
      method = "map"
    ),
    "person_id_col 'person_id' not found in data columns"
  )
})

test_that("fscores() validates item parameters are provided for all items", {
  wide_data <- data.frame(
    person_id = 1:3,
    item1 = c(0, 1, 0),
    item2 = c(1, 0, 1)
  )

  # Missing parameters for item2
  item_params <- list(
    item1 = list(loadings = c(1.0), thresholds = c(0.5))
  )

  expect_error(
    fscores(
      data = wide_data,
      item_params = item_params,
      person_id_col = "person_id",
      method = "map"
    ),
    "Missing item parameters for columns: item2"
  )
})

test_that("fscores() validates all items have same number of dimensions", {
  wide_data <- data.frame(
    person_id = 1:3,
    item1 = c(0, 1, 0),
    item2 = c(1, 0, 1)
  )

  # item1 has 2 dimensions, item2 has 1 dimension (mismatch)
  item_params <- list(
    item1 = list(
      loadings = c(1.0, 0.5),  # 2 dimensions
      thresholds = c(0.5)
    ),
    item2 = list(
      loadings = c(0.8),  # 1 dimension
      thresholds = c(-0.5)
    )
  )

  expect_error(
    fscores(
      data = wide_data,
      item_params = item_params,
      person_id_col = "person_id",
      method = "map"
    ),
    "Item 'item2' has 1 loadings but expected 2"
  )
})

test_that("fscores() errors when no item columns found", {
  wide_data <- data.frame(
    person_id = 1:3,
    age = c(25, 30, 35)
  )

  item_params <- list()

  expect_error(
    fscores(
      data = wide_data,
      item_params = item_params,
      person_id_col = "person_id",
      covariate_cols = "age",
      method = "map"
    ),
    "No item columns found"
  )
})

test_that("fscores() errors when no valid responses found", {
  wide_data <- data.frame(
    person_id = 1:3,
    item1 = c(NA, NA, NA)  # All missing
  )

  item_params <- list(
    item1 = list(loadings = c(1.0), thresholds = c(0.5))
  )

  expect_error(
    fscores(
      data = wide_data,
      item_params = item_params,
      person_id_col = "person_id",
      method = "map"
    ),
    "No valid responses found in data"
  )
})

# ============================================================================
# Tests for fscores() - Edge Cases
# ============================================================================

test_that("fscores() works with single person", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  wide_data <- data.frame(
    person_id = 1,
    item1 = 0,
    item2 = 1
  )

  item_params <- list(
    item1 = list(loadings = c(1.2), thresholds = c(0.5)),
    item2 = list(loadings = c(0.8), thresholds = c(-0.5))
  )

  fit <- fscores(
    data = wide_data,
    item_params = item_params,
    person_id_col = "person_id",
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
  expect_true(nrow(fit$theta) == 1)
})

test_that("fscores() works with single item", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  wide_data <- data.frame(
    person_id = 1:3,
    item1 = c(0, 1, 0)
  )

  item_params <- list(
    item1 = list(loadings = c(1.2), thresholds = c(0.5))
  )

  fit <- fscores(
    data = wide_data,
    item_params = item_params,
    person_id_col = "person_id",
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
  expect_true(nrow(fit$theta) == 3)
})

test_that("fscores() handles custom person_id_col name", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  wide_data <- data.frame(
    student_id = 1:3,  # Custom column name
    item1 = c(0, 1, 0)
  )

  item_params <- list(
    item1 = list(loadings = c(1.2), thresholds = c(0.5))
  )

  fit <- fscores(
    data = wide_data,
    item_params = item_params,
    person_id_col = "student_id",
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
})

test_that("fscores() works with factor covariates", {
  skip_if_not(has_cmdstanr() || has_rstan(), "No Stan backend available")

  wide_data <- data.frame(
    person_id = 1:4,
    group = factor(c("A", "B", "A", "B")),
    item1 = c(0, 1, 0, 1)
  )

  item_params <- list(
    item1 = list(loadings = c(1.2), thresholds = c(0.5))
  )

  fit <- fscores(
    data = wide_data,
    item_params = item_params,
    person_id_col = "person_id",
    covariate_cols = "group",
    method = "map"
  )

  expect_s3_class(fit, "irt_fit")
  # model.matrix() should create intercept + dummy variable
  expect_true(ncol(fit$b) == 2)
})
