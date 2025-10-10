# Tests for data preparation and validation functions

test_that("prepare_irt_data handles simple unidimensional case", {
  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )
  loadings <- matrix(c(1, 1.5, 1, 1.5), ncol = 1)  # L x 1 matrix
  thresh_L <- c(-Inf, -1, -Inf, -1)  # Length L
  thresh_R <- c(-1, Inf, -1, Inf)    # Length L

  result <- prepare_irt_data(response_data, loadings, thresh_L, thresh_R)

  expect_equal(result$K, 1)
  expect_equal(result$L, 4)
  expect_equal(result$N, 2)
  expect_equal(result$inf, -999)
  # Check that infinities were replaced
  expect_equal(result$dL[1], -999)  # Was -Inf
  expect_equal(result$dR[2], -999)  # Was Inf
})

test_that("prepare_irt_data handles multidimensional case", {
  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )
  # L x K loading matrix (4 responses x 2 dimensions)
  loadings <- matrix(c(1, 0.5, 1, 0.5,    # Dim 1
                       0.5, 1, 0.5, 1),    # Dim 2
                     ncol = 2)
  thresh_L <- c(-Inf, -1, -Inf, -1)
  thresh_R <- c(-1, Inf, -1, Inf)

  result <- prepare_irt_data(response_data, loadings, thresh_L, thresh_R)

  expect_equal(result$K, 2)
  expect_equal(ncol(result$A), 2)
  expect_equal(nrow(result$A), 4)
})

test_that("prepare_irt_data handles person covariates", {
  response_data <- data.frame(
    pid = c(1, 1, 2, 2, 3, 3),
    iid = c(1, 2, 1, 2, 1, 2)
  )
  loadings <- matrix(rep(1, 6), ncol = 1)
  thresh_L <- rep(-Inf, 6)
  thresh_R <- rep(Inf, 6)

  # With covariates
  covariates <- matrix(c(1, 1, 1,    # Intercept
                         25, 30, 35), # Age
                       ncol = 2)

  result <- prepare_irt_data(
    response_data, loadings, thresh_L, thresh_R,
    person_covariates = covariates
  )

  expect_equal(result$J, 2)
  expect_equal(nrow(result$X), 3)
  expect_equal(ncol(result$X), 2)
})

test_that("prepare_irt_data creates intercept-only covariate matrix by default", {
  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )
  loadings <- matrix(rep(1, 4), ncol = 1)
  thresh_L <- rep(-Inf, 4)
  thresh_R <- rep(Inf, 4)

  result <- prepare_irt_data(response_data, loadings, thresh_L, thresh_R)

  expect_equal(result$J, 1)
  expect_true(all(result$X == 1))
})

test_that("prepare_irt_data handles custom weights", {
  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )
  loadings <- matrix(rep(1, 4), ncol = 1)
  thresh_L <- rep(-Inf, 4)
  thresh_R <- rep(Inf, 4)
  weights <- c(1.5, 0.5)

  result <- prepare_irt_data(
    response_data, loadings, thresh_L, thresh_R,
    weights = weights
  )

  expect_equal(result$wgt, weights)
})

test_that("prepare_irt_data creates default weights", {
  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )
  loadings <- matrix(rep(1, 4), ncol = 1)
  thresh_L <- rep(-Inf, 4)
  thresh_R <- rep(Inf, 4)

  result <- prepare_irt_data(response_data, loadings, thresh_L, thresh_R)

  expect_equal(result$wgt, c(1, 1))
})

test_that("validate_response_data rejects missing required columns", {
  bad_data <- data.frame(person = 1:10)  # Missing pid and iid

  expect_error(
    validate_response_data(bad_data),
    "must.include"
  )
})

test_that("validate_response_data rejects non-numeric person IDs", {
  bad_data <- data.frame(
    pid = c("A", "B"),
    iid = c(1, 2)
  )

  expect_error(
    validate_response_data(bad_data),
    "Person IDs.*must be numeric"
  )
})

test_that("validate_response_data rejects negative person IDs", {
  bad_data <- data.frame(
    pid = c(-1, 1),
    iid = c(1, 2)
  )

  expect_error(
    validate_response_data(bad_data),
    "Person IDs.*must be >= 1"
  )
})

test_that("validate_response_data rejects non-integer person IDs", {
  bad_data <- data.frame(
    pid = c(1.5, 2.5),
    iid = c(1, 2)
  )

  expect_error(
    validate_response_data(bad_data),
    "Person IDs.*must be integers"
  )
})

test_that("validate_response_data warns about non-consecutive person IDs", {
  data_with_gaps <- data.frame(
    pid = c(1, 1, 5, 5),  # Gap between 1 and 5
    iid = c(1, 2, 1, 2)
  )

  expect_warning(
    validate_response_data(data_with_gaps),
    "not consecutive"
  )
})

test_that("validate_item_loadings rejects wrong number of columns", {
  loadings <- matrix(1:6, ncol = 3)

  expect_error(
    validate_item_loadings(loadings, K = 2),
    "must have K columns"
  )
})

test_that("validate_item_loadings rejects missing values", {
  loadings <- matrix(c(1, NA, 3, 4), ncol = 2)

  expect_error(
    validate_item_loadings(loadings, K = 2),
    "cannot contain missing"
  )
})

test_that("validate_thresholds rejects wrong length", {
  left <- c(-Inf, -1)
  right <- c(0, Inf)

  expect_error(
    validate_thresholds(left, right, L = 5),
    "len"
  )
})

test_that("validate_thresholds rejects left > right", {
  left <- c(-1, 1)    # Second element: left=1 > right=0 (violation)
  right <- c(0, 0)

  expect_error(
    validate_thresholds(left, right, L = 2),
    "Left threshold must be <= right"
  )
})

test_that("validate_thresholds warns when both thresholds are infinite", {
  left <- c(-Inf, -Inf)
  right <- c(Inf, Inf)

  expect_warning(
    validate_thresholds(left, right, L = 2),
    "both.*infinite"
  )
})

test_that("prepare_irt_data rejects wrong weight length", {
  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )
  loadings <- matrix(rep(1, 4), ncol = 1)
  thresh_L <- rep(-Inf, 4)
  thresh_R <- rep(Inf, 4)
  bad_weights <- c(1, 0.5, 0.3)  # Should be length 2

  expect_error(
    prepare_irt_data(response_data, loadings, thresh_L, thresh_R,
                     weights = bad_weights),
    "weights must have length N"
  )
})

test_that("prepare_irt_data rejects negative weights", {
  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )
  loadings <- matrix(rep(1, 4), ncol = 1)
  thresh_L <- rep(-Inf, 4)
  thresh_R <- rep(Inf, 4)
  bad_weights <- c(-1, 1)

  expect_error(
    prepare_irt_data(response_data, loadings, thresh_L, thresh_R,
                     weights = bad_weights),
    "weights must be non-negative"
  )
})

test_that("prepare_irt_data converts covariate vector to matrix", {
  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )
  loadings <- matrix(rep(1, 4), ncol = 1)
  thresh_L <- rep(-Inf, 4)
  thresh_R <- rep(Inf, 4)
  covariates <- c(25, 30)  # Vector instead of matrix

  result <- prepare_irt_data(
    response_data, loadings, thresh_L, thresh_R,
    person_covariates = covariates
  )

  expect_true(is.matrix(result$X))
  expect_equal(ncol(result$X), 1)
})

test_that("prepare_irt_data rejects wrong covariate matrix size", {
  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )
  loadings <- matrix(rep(1, 4), ncol = 1)
  thresh_L <- rep(-Inf, 4)
  thresh_R <- rep(Inf, 4)
  bad_covariates <- matrix(c(1, 1, 1), ncol = 1)  # Should have 2 rows

  expect_error(
    prepare_irt_data(response_data, loadings, thresh_L, thresh_R,
                     person_covariates = bad_covariates),
    "person_covariates must have N rows"
  )
})

test_that("prepare_irt_data produces complete Stan data structure", {
  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )
  loadings <- matrix(c(1, 1.5, 1, 1.5), ncol = 1)  # L x 1 matrix
  thresh_L <- c(-Inf, -1, -Inf, -1)
  thresh_R <- c(-1, Inf, -1, Inf)

  result <- prepare_irt_data(response_data, loadings, thresh_L, thresh_R)

  # Check all required Stan data components are present
  expected_names <- c("inf", "L", "N", "J", "K", "pid", "dL", "dR", "wgt", "X", "A")
  expect_true(all(expected_names %in% names(result)))

  # Check dimensions are consistent
  expect_equal(length(result$pid), result$L)
  expect_equal(length(result$dL), result$L)
  expect_equal(length(result$dR), result$L)
  expect_equal(length(result$wgt), result$N)
  expect_equal(nrow(result$X), result$N)
  expect_equal(ncol(result$X), result$J)
  expect_equal(nrow(result$A), result$L)
  expect_equal(ncol(result$A), result$K)
})
