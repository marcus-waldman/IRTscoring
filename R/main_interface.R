# Main User Interface Functions
#
# High-level functions for fitting IRT models with user-friendly interfaces

#' @importFrom stats model.matrix as.formula
NULL

#' Fit IRT Model for Ability Scoring
#'
#' Main interface for fitting unidimensional or multidimensional IRT models.
#' This function prepares data and routes to the appropriate estimation method.
#'
#' @param response_data Data frame with columns \code{pid} (person ID) and
#'   \code{iid} (item ID). Person IDs should be consecutive integers from 1 to N.
#' @param item_loadings Matrix of factor loadings (L x K), where L is the number
#'   of response observations and K is the number of dimensions.
#' @param threshold_left Numeric vector of left thresholds for each response (length L).
#'   Use \code{-Inf} for left-unbounded categories.
#' @param threshold_right Numeric vector of right thresholds for each response (length L).
#'   Use \code{Inf} for right-unbounded categories.
#' @param person_covariates Optional matrix of person-level covariates (N x J).
#'   If NULL, an intercept-only model is used.
#' @param weights Optional vector of person weights (length N). If NULL, all
#'   persons receive equal weight of 1.
#' @param method Estimation method to use:
#'   \describe{
#'     \item{"map"}{Maximum A Posteriori (fastest)}
#'     \item{"laplace"}{Laplace approximation (MAP + uncertainty)}
#'     \item{"vb"}{Variational Bayes}
#'     \item{"mcmc"}{Full Bayesian MCMC (slowest, most accurate)}
#'     \item{"optim"}{Base R optim() via BridgeStan (requires setup)}
#'     \item{"nloptr"}{nloptr algorithms via BridgeStan (requires setup)}
#'   }
#' @param ... Additional arguments passed to the specific estimation function.
#'   See \code{\link{fit_irt_map}}, \code{\link{fit_irt_laplace}},
#'   \code{\link{fit_irt_vb}}, \code{\link{fit_irt_mcmc}},
#'   \code{\link{fit_irt_optim}}, or \code{\link{fit_irt_nloptr}} for details.
#'
#' @return Object of class "irt_fit" containing:
#'   \item{theta}{N x K matrix of estimated ability scores}
#'   \item{beta0}{Length K vector of dimension intercepts}
#'   \item{b}{K x J matrix of covariate effects}
#'   \item{tau}{Length K vector of dimension standard deviations}
#'   \item{Omega}{K x K correlation matrix}
#'   \item{log_lik}{Log-likelihood value}
#'   \item{method}{Estimation method used}
#'   \item{backend}{Backend used (cmdstanr, rstan, or bridgestan)}
#'   \item{stan_fit}{Original fit object}
#'
#' @details
#' This is the main user-facing function for fitting IRT models. It handles
#' data preparation and routes to the appropriate estimation method based on
#' the \code{method} argument.
#'
#' **Recommended methods:**
#' - \code{"map"}: Fast point estimates, good for most applications
#' - \code{"laplace"}: MAP plus uncertainty quantification
#' - \code{"vb"}: Fast approximate Bayesian inference
#' - \code{"mcmc"}: Complete posterior distributions, gold standard
#'
#' **Advanced methods** (require BridgeStan setup):
#' - \code{"optim"}: Access to optim() algorithms
#' - \code{"nloptr"}: Access to NLopt algorithms
#'
#' @examples
#' \dontrun{
#' # Unidimensional IRT with MAP estimation
#' response_data <- data.frame(
#'   pid = c(1, 1, 2, 2, 3, 3),
#'   iid = c(1, 2, 1, 2, 1, 2)
#' )
#' loadings <- matrix(c(1, 1.5, 1, 1.5, 2, 1.2), ncol = 1)
#' thresh_L <- c(-Inf, -1, -Inf, -1, -Inf, -2)
#' thresh_R <- c(-1, Inf, -1, Inf, -2, Inf)
#'
#' fit <- fit_ability(
#'   response_data = response_data,
#'   item_loadings = loadings,
#'   threshold_left = thresh_L,
#'   threshold_right = thresh_R,
#'   method = "map"
#' )
#'
#' # Multidimensional IRT with MCMC
#' loadings_2d <- matrix(c(1, 0.5, 1, 0.5, 0.5, 1), ncol = 2)
#' fit_mcmc <- fit_ability(
#'   response_data = response_data,
#'   item_loadings = loadings_2d,
#'   threshold_left = thresh_L,
#'   threshold_right = thresh_R,
#'   method = "mcmc",
#'   chains = 4,
#'   iter = 2000
#' )
#' }
#'
#' @export
fit_ability <- function(response_data,
                        item_loadings,
                        threshold_left,
                        threshold_right,
                        person_covariates = NULL,
                        weights = NULL,
                        method = c("map", "laplace", "vb", "mcmc", "optim", "nloptr"),
                        ...) {

  # Validate and match method
  method <- match.arg(method)

  # Prepare data for Stan
  stan_data <- prepare_irt_data(
    response_data = response_data,
    item_loadings = item_loadings,
    threshold_left = threshold_left,
    threshold_right = threshold_right,
    person_covariates = person_covariates,
    weights = weights
  )

  # Route to appropriate estimation method
  fit <- switch(
    method,
    map = fit_irt_map(stan_data, ...),
    laplace = fit_irt_laplace(stan_data, ...),
    vb = fit_irt_vb(stan_data, ...),
    mcmc = fit_irt_mcmc(stan_data, ...),
    optim = fit_irt_optim(stan_data, ...),
    nloptr = fit_irt_nloptr(stan_data, ...),
    stop("Unknown method: ", method, ". This should not happen.")
  )

  return(fit)
}

#' Fit Unidimensional IRT Model
#'
#' Convenience wrapper for unidimensional IRT models. Accepts discrimination
#' parameters as a vector and converts to the L x 1 loading matrix format.
#'
#' @param response_data Data frame with columns \code{pid} and \code{iid}
#' @param item_discrimination Numeric vector of discrimination parameters (length L)
#' @param threshold_left Numeric vector of left thresholds (length L)
#' @param threshold_right Numeric vector of right thresholds (length L)
#' @param ... Additional arguments passed to \code{\link{fit_ability}}
#'
#' @return Object of class "irt_fit"
#'
#' @examples
#' \dontrun{
#' response_data <- data.frame(
#'   pid = c(1, 1, 2, 2),
#'   iid = c(1, 2, 1, 2)
#' )
#'
#' fit <- fit_unidim_ability(
#'   response_data = response_data,
#'   item_discrimination = c(1, 1.5, 1, 1.5),
#'   threshold_left = c(-Inf, -1, -Inf, -1),
#'   threshold_right = c(-1, Inf, -1, Inf),
#'   method = "map"
#' )
#' }
#'
#' @export
fit_unidim_ability <- function(response_data,
                                item_discrimination,
                                threshold_left,
                                threshold_right,
                                ...) {

  # Convert discrimination vector to L x 1 loading matrix
  item_loadings <- matrix(item_discrimination, ncol = 1)

  # Call main function
  fit_ability(
    response_data = response_data,
    item_loadings = item_loadings,
    threshold_left = threshold_left,
    threshold_right = threshold_right,
    ...
  )
}

#' Calculate Factor Scores Using IRT Model
#'
#' User-friendly function that accepts data in wide format (one row per person)
#' and returns IRT ability scores. Provides an interface similar to \code{mirt::fscores()}.
#'
#' @param data Data frame in wide format with one row per person. Columns should
#'   include a person ID column, optional covariate columns, and item response columns.
#' @param item_params Named list of item parameters. Each element should be a list
#'   containing:
#'   \describe{
#'     \item{loadings}{Numeric vector of length K (number of dimensions)}
#'     \item{thresholds}{Numeric vector of thresholds (length = n_categories - 1)}
#'   }
#' @param person_id_col Name of the column containing person IDs (default: "person_id")
#' @param covariate_cols Optional character vector of column names to use as covariates
#' @param method Estimation method (see \code{\link{fit_ability}} for options)
#' @param weights Optional vector of person weights
#' @param ... Additional arguments passed to \code{\link{fit_ability}}
#'
#' @return Object of class "irt_fit"
#'
#' @details
#' **CRITICAL: Response values must be 0-indexed!**
#'
#' Item responses must start at 0, not 1. For an item with 3 categories,
#' valid responses are 0, 1, 2 (not 1, 2, 3).
#'
#' **Threshold specification:**
#'
#' Thresholds define boundaries between categories:
#' - Category 0: (-Inf, threshold\[1\]]
#' - Category 1: (threshold\[1\], threshold\[2\]]
#' - Category 2: (threshold\[2\], Inf)
#'
#' For k+1 categories (0 through k), you need k thresholds.
#'
#' **Validation:**
#'
#' The function validates:
#' - Response values are >= 0 (0-indexed)
#' - Sufficient thresholds provided for observed responses
#' - All items have parameter specifications
#'
#' @examples
#' \dontrun{
#' # Create wide format data (responses are 0-indexed!)
#' wide_data <- data.frame(
#'   person_id = 1:100,
#'   age = rnorm(100, 30, 10),
#'   item1 = sample(0:2, 100, replace = TRUE),  # 0-indexed!
#'   item2 = sample(0:2, 100, replace = TRUE)   # 0-indexed!
#' )
#'
#' # Define item parameters
#' item_params <- list(
#'   item1 = list(
#'     loadings = c(1.2, 0.1),
#'     thresholds = c(-1, 0.5)  # 2 thresholds for 3 categories
#'   ),
#'   item2 = list(
#'     loadings = c(0.8, 1.5),
#'     thresholds = c(-0.5, 1.0)
#'   )
#' )
#'
#' # Fit model
#' fit <- fscores(
#'   data = wide_data,
#'   item_params = item_params,
#'   person_id_col = "person_id",
#'   covariate_cols = "age",
#'   method = "map"
#' )
#'
#' # Extract ability scores
#' theta <- fit$theta
#' }
#'
#' @export
fscores <- function(data,
                    item_params,
                    person_id_col = "person_id",
                    covariate_cols = NULL,
                    method = "map",
                    weights = NULL,
                    ...) {

  # Validate inputs
  if (!person_id_col %in% colnames(data)) {
    stop("person_id_col '", person_id_col, "' not found in data columns: ",
         paste(colnames(data), collapse = ", "))
  }

  # Identify item columns (exclude person ID and covariates)
  all_cols <- colnames(data)
  item_cols <- setdiff(all_cols, c(person_id_col, covariate_cols))

  if (length(item_cols) == 0) {
    stop("No item columns found. All columns are either person_id_col or covariate_cols.")
  }

  # Check that all item columns have parameters
  missing_params <- setdiff(item_cols, names(item_params))
  if (length(missing_params) > 0) {
    stop("Missing item parameters for columns: ", paste(missing_params, collapse = ", "))
  }

  # Get number of dimensions from first item
  K <- length(item_params[[item_cols[1]]]$loadings)

  # Validate all items have same number of dimensions
  for (item_name in item_cols) {
    if (length(item_params[[item_name]]$loadings) != K) {
      stop("Item '", item_name, "' has ", length(item_params[[item_name]]$loadings),
           " loadings but expected ", K, " (from first item)")
    }
  }

  # Convert wide to long format
  response_list <- list()
  loading_list <- list()
  thresh_left_list <- list()
  thresh_right_list <- list()

  for (i in 1:nrow(data)) {
    person_id <- data[[person_id_col]][i]

    for (item_name in item_cols) {
      response_value <- data[[item_name]][i]

      if (!is.na(response_value)) {
        # Validate response is 0-indexed
        if (response_value < 0) {
          stop("Response values must be >= 0 (0-indexed). ",
               "Found response ", response_value, " for person ", person_id,
               ", item '", item_name, "'. ",
               "Categories should be coded as 0, 1, 2, ... (not 1, 2, 3, ...)")
        }

        # Get item parameters
        item_info <- item_params[[item_name]]

        # Validate sufficient thresholds for observed response
        # For response value k (0-indexed), we need at least k thresholds
        n_thresholds_needed <- response_value

        if (length(item_info$thresholds) < n_thresholds_needed) {
          stop("Item '", item_name, "' has observed response ", response_value,
               " (category ", response_value, "), which requires at least ",
               n_thresholds_needed, " threshold(s), but only ",
               length(item_info$thresholds), " provided. ",
               "For categories 0 through ", response_value, ", you need ",
               response_value, " threshold(s).")
        }

        # Create threshold boundaries: (-Inf, t1, t2, ..., tn, Inf)
        thresholds <- c(-Inf, item_info$thresholds, Inf)

        # For 0-indexed response k, category spans (thresholds[k+1], thresholds[k+2]]
        # (adding 1 because R vectors are 1-indexed)
        response_list <- c(response_list, list(data.frame(
          pid = person_id,
          iid = match(item_name, item_cols)
        )))

        loading_list <- c(loading_list, list(item_info$loadings))
        thresh_left_list <- c(thresh_left_list, thresholds[response_value + 1])
        thresh_right_list <- c(thresh_right_list, thresholds[response_value + 2])
      }
    }
  }

  # Check if any responses were found
  if (length(response_list) == 0) {
    stop("No valid responses found in data. Check for missing values or incorrect column names.")
  }

  # Combine into matrices
  response_data <- do.call(rbind, response_list)
  item_loadings <- do.call(rbind, loading_list)
  threshold_left <- unlist(thresh_left_list)
  threshold_right <- unlist(thresh_right_list)

  # Prepare covariates
  person_covariates <- NULL
  if (!is.null(covariate_cols)) {
    # Get unique persons
    unique_pids <- sort(unique(response_data$pid))
    N <- length(unique_pids)

    # Create design matrix with intercept
    covariate_data <- data[data[[person_id_col]] %in% unique_pids,
                           c(person_id_col, covariate_cols), drop = FALSE]

    # Order by person ID to match response_data
    covariate_data <- covariate_data[order(covariate_data[[person_id_col]]), ]

    # Create model matrix (handles factors automatically, includes intercept)
    formula_str <- paste("~", paste(covariate_cols, collapse = " + "))
    person_covariates <- model.matrix(as.formula(formula_str), data = covariate_data)
  }

  # Prepare weights
  if (is.null(weights)) {
    weights <- rep(1, length(unique(response_data$pid)))
  }

  # Call the main function
  fit_ability(
    response_data = response_data,
    item_loadings = item_loadings,
    threshold_left = threshold_left,
    threshold_right = threshold_right,
    person_covariates = person_covariates,
    weights = weights,
    method = method,
    ...
  )
}
