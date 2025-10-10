# Data Preparation and Validation Functions
#
# Functions for preparing and validating IRT data before Stan estimation

#' Prepare IRT Data for Estimation
#'
#' Validates and formats IRT data for Stan model estimation. Converts response
#' data and parameters into the format required by the Stan model.
#'
#' @param response_data Data frame with columns: pid (person ID), iid (item ID).
#'   Person IDs should be consecutive integers from 1 to N.
#' @param item_loadings Matrix of factor loadings (L x K), where L is the number
#'   of response observations and K is the number of dimensions.
#' @param threshold_left Numeric vector of left thresholds for each response (length L).
#'   Use Inf for -infinity boundaries.
#' @param threshold_right Numeric vector of right thresholds for each response (length L).
#'   Use Inf for +infinity boundaries.
#' @param person_covariates Optional matrix of person-level covariates (N x J).
#'   If NULL, an intercept-only model is used.
#' @param weights Optional vector of person weights (length N). If NULL, all
#'   persons receive equal weight of 1.
#' @param K Optional number of dimensions. If NULL, auto-detected from item_loadings.
#'
#' @return List formatted for Stan model with components:
#'   \item{inf}{Code for infinity (-999)}
#'   \item{L}{Number of response observations}
#'   \item{N}{Number of persons}
#'   \item{J}{Number of covariates}
#'   \item{K}{Number of dimensions}
#'   \item{pid}{Person ID vector}
#'   \item{dL}{Left thresholds (infinities replaced with -999)}
#'   \item{dR}{Right thresholds (infinities replaced with -999)}
#'   \item{wgt}{Person weights}
#'   \item{X}{Covariate matrix}
#'   \item{A}{Loading matrix}
#'
#' @details
#' The function performs comprehensive validation:
#' - Checks that response_data has required columns (pid, iid)
#' - Validates that person IDs are consecutive integers
#' - Ensures item_loadings is a matrix with correct dimensions
#' - Validates threshold vectors have correct length
#' - Checks that left thresholds <= right thresholds
#' - Replaces infinite thresholds with code -999 for Stan
#'
#' @examples
#' \dontrun{
#' # Simple unidimensional example
#' response_data <- data.frame(
#'   pid = c(1, 1, 2, 2),
#'   iid = c(1, 2, 1, 2)
#' )
#' loadings <- matrix(c(1, 1.5), ncol = 1)
#' thresh_L <- c(-Inf, -1)
#' thresh_R <- c(-1, Inf)
#'
#' stan_data <- prepare_irt_data(
#'   response_data = response_data,
#'   item_loadings = loadings,
#'   threshold_left = thresh_L,
#'   threshold_right = thresh_R
#' )
#' }
#'
#' @export
prepare_irt_data <- function(response_data,
                              item_loadings,
                              threshold_left,
                              threshold_right,
                              person_covariates = NULL,
                              weights = NULL,
                              K = NULL) {

  # Validate inputs
  validate_response_data(response_data)

  # Auto-detect K if not provided
  if (is.null(K)) {
    K <- ncol(item_loadings)
  }

  validate_item_loadings(item_loadings, K)
  validate_thresholds(threshold_left, threshold_right, nrow(response_data))

  # Get dimensions
  L <- nrow(response_data)
  N <- length(unique(response_data$pid))

  # Handle infinite thresholds
  # Map -Inf to -999 and +Inf to +999 to preserve probability bounds
  inf_code <- 999
  threshold_left[threshold_left == -Inf] <- -inf_code
  threshold_left[threshold_left == Inf] <- inf_code
  threshold_right[threshold_right == -Inf] <- -inf_code
  threshold_right[threshold_right == Inf] <- inf_code

  # Create person covariate matrix
  if (is.null(person_covariates)) {
    # Intercept-only model
    person_covariates <- matrix(1, nrow = N, ncol = 1)
  } else {
    # Validate covariate matrix
    if (!is.matrix(person_covariates)) {
      person_covariates <- as.matrix(person_covariates)
    }
    if (nrow(person_covariates) != N) {
      stop("person_covariates must have N rows (one per person)")
    }
  }
  J <- ncol(person_covariates)

  # Create weights
  if (is.null(weights)) {
    weights <- rep(1, N)
  } else {
    if (length(weights) != N) {
      stop("weights must have length N (one per person)")
    }
    if (any(weights < 0)) {
      stop("weights must be non-negative")
    }
  }

  # Format for Stan
  stan_data <- list(
    inf = inf_code,  # Not actually used by Stan, but kept for reference
    L = L,
    N = N,
    J = J,
    K = K,
    pid = response_data$pid,
    dL = threshold_left,
    dR = threshold_right,
    wgt = weights,
    X = person_covariates,
    A = item_loadings
  )

  return(stan_data)
}

#' Validate Response Data
#'
#' Checks that response data has the required structure
#'
#' @param data Data frame to validate
#' @keywords internal
validate_response_data <- function(data) {
  checkmate::assert_data_frame(data)
  checkmate::assert_names(colnames(data), must.include = c("pid", "iid"))

  # Check person IDs are valid
  if (!is.numeric(data$pid)) {
    stop("Person IDs (pid) must be numeric")
  }

  if (any(data$pid < 1)) {
    stop("Person IDs (pid) must be >= 1")
  }

  # Check that person IDs are integers
  if (!all(data$pid == floor(data$pid))) {
    stop("Person IDs (pid) must be integers")
  }

  # Check item IDs are valid
  if (!is.numeric(data$iid)) {
    stop("Item IDs (iid) must be numeric")
  }

  if (any(data$iid < 1)) {
    stop("Item IDs (iid) must be >= 1")
  }

  # Warn if person IDs are not consecutive
  unique_pids <- sort(unique(data$pid))
  expected_pids <- seq_len(length(unique_pids))
  if (!identical(unique_pids, expected_pids)) {
    warning("Person IDs are not consecutive integers from 1 to N. ",
            "This may cause issues with covariate alignment.")
  }

  invisible(TRUE)
}

#' Validate Item Loadings
#'
#' Checks that item loadings matrix has the correct structure
#'
#' @param loadings Loading matrix to validate
#' @param K Number of dimensions
#' @keywords internal
validate_item_loadings <- function(loadings, K) {
  checkmate::assert_matrix(loadings, mode = "numeric")

  if (!is.null(K)) {
    if (ncol(loadings) != K) {
      stop("item_loadings must have K columns (one per dimension). ",
           "Expected ", K, " columns but got ", ncol(loadings))
    }
  }

  if (any(is.na(loadings))) {
    stop("item_loadings cannot contain missing values (NA)")
  }

  invisible(TRUE)
}

#' Validate Thresholds
#'
#' Checks that threshold vectors have correct structure and ordering
#'
#' @param left Left thresholds
#' @param right Right thresholds
#' @param L Number of responses
#' @keywords internal
validate_thresholds <- function(left, right, L) {
  checkmate::assert_numeric(left, len = L, any.missing = FALSE)
  checkmate::assert_numeric(right, len = L, any.missing = FALSE)

  # Check that left <= right (except for infinities)
  finite_both <- is.finite(left) & is.finite(right)
  if (any(left[finite_both] > right[finite_both])) {
    bad_idx <- which(finite_both & (left > right))[1]
    stop("Left threshold must be <= right threshold. ",
         "Violation at index ", bad_idx, ": ",
         "left=", left[bad_idx], ", right=", right[bad_idx])
  }

  # Warn if both thresholds are infinite
  both_inf <- is.infinite(left) & is.infinite(right)
  if (any(both_inf)) {
    warning("Some responses have both left and right thresholds as infinite. ",
            "This creates a degenerate probability of 1.")
  }

  invisible(TRUE)
}
