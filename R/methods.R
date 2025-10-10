# S3 Methods for irt_fit Objects
#
# Print, summary, coef, and predict methods for IRT model fits

#' Print IRT Model Fit
#'
#' @param x An object of class "irt_fit"
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the input object
#' @export
print.irt_fit <- function(x, ...) {
  cat("IRT Model Fit\n")
  cat("=============\n\n")
  cat("Method:          ", x$method, "\n", sep = "")
  cat("Backend:         ", x$backend, "\n", sep = "")
  cat("Dimensions:      ", ncol(x$theta), "\n", sep = "")
  cat("Persons:         ", nrow(x$theta), "\n", sep = "")
  cat("Log-likelihood:  ", round(x$log_lik, 2), "\n", sep = "")

  cat("\nAbility Score Summary:\n")
  print(summary(as.vector(x$theta)))

  invisible(x)
}

#' Summary of IRT Model Fit
#'
#' @param object An object of class "irt_fit"
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the input object
#' @export
summary.irt_fit <- function(object, ...) {
  cat("IRT Model Summary\n")
  cat("=================\n\n")

  cat("Estimation Method: ", object$method, "\n", sep = "")
  cat("Backend:           ", object$backend, "\n", sep = "")
  cat("Number of Persons: ", nrow(object$theta), "\n", sep = "")
  cat("Dimensions:        ", ncol(object$theta), "\n", sep = "")
  cat("Log-likelihood:    ", round(object$log_lik, 2), "\n\n", sep = "")

  cat("Dimension Intercepts (beta0):\n")
  print(round(object$beta0, 3))
  cat("\n")

  if (ncol(object$b) > 1 || any(object$b[, 1] != 0)) {
    cat("Covariate Effects (b):\n")
    print(round(object$b, 3))
    cat("\n")
  }

  cat("Dimension Standard Deviations (tau):\n")
  print(round(object$tau, 3))
  cat("\n")

  if (ncol(object$theta) > 1) {
    cat("Dimension Correlations (Omega):\n")
    print(round(object$Omega, 3))
    cat("\n")
  }

  cat("Ability Estimates (theta):\n")
  theta_summary <- apply(object$theta, 2, function(x) {
    c(Mean = mean(x), SD = sd(x), Min = min(x), Max = max(x))
  })
  colnames(theta_summary) <- paste0("Dim", 1:ncol(object$theta))
  print(round(theta_summary, 3))

  invisible(object)
}

#' Extract Coefficients from IRT Model Fit
#'
#' @param object An object of class "irt_fit"
#' @param pars Character vector of parameter names to extract.
#'   Options: "beta0", "b", "tau", "Omega", "theta"
#' @param ... Additional arguments (currently unused)
#'
#' @return Named list of extracted parameters
#' @export
coef.irt_fit <- function(object,
                         pars = c("beta0", "b", "tau", "Omega"),
                         ...) {

  pars <- match.arg(pars, several.ok = TRUE,
                    choices = c("beta0", "b", "tau", "Omega", "theta"))

  result <- list()

  if ("beta0" %in% pars) {
    result$beta0 <- object$beta0
  }

  if ("b" %in% pars) {
    result$b <- object$b
  }

  if ("tau" %in% pars) {
    result$tau <- object$tau
  }

  if ("Omega" %in% pars) {
    result$Omega <- object$Omega
  }

  if ("theta" %in% pars) {
    result$theta <- object$theta
  }

  if (length(result) == 1) {
    return(result[[1]])
  }

  return(result)
}

#' Predict Method for IRT Model Fits
#'
#' Extract ability estimates and other predictions from fitted IRT models
#'
#' @param object An object of class "irt_fit"
#' @param type Type of prediction to return:
#'   \describe{
#'     \item{"abilities"}{Full ability scores (theta = beta0 + X*b + zeta)}
#'     \item{"person_scores"}{Same as "abilities" (alias)}
#'   }
#' @param ... Additional arguments (currently unused)
#'
#' @return Matrix of predictions (N x K)
#'
#' @details
#' The "abilities" or "person_scores" type returns the complete ability
#' estimates including both the fixed effects (intercepts and covariate effects)
#' and the random effects (individual deviations).
#'
#' @examples
#' \dontrun{
#' fit <- fit_ability(response_data, item_loadings,
#'                    threshold_left, threshold_right, method = "map")
#'
#' # Get ability estimates
#' abilities <- predict(fit)
#' abilities <- predict(fit, type = "abilities")
#' }
#'
#' @export
predict.irt_fit <- function(object,
                            type = c("abilities", "person_scores"),
                            ...) {

  type <- match.arg(type)

  if (type %in% c("abilities", "person_scores")) {
    return(object$theta)
  }
}

# ============================================================================
# Utility Functions for Parameter Extraction
# ============================================================================

#' Get Dimension Correlations from IRT Fit
#'
#' Extract the correlation matrix between dimensions (Omega)
#'
#' @param fit An object of class "irt_fit"
#'
#' @return K x K correlation matrix
#' @keywords internal
get_dimension_correlations <- function(fit) {
  K <- ncol(fit$theta)

  if (K == 1) {
    return(matrix(1, 1, 1))
  }

  return(fit$Omega)
}

#' Get Covariate Effects from IRT Fit
#'
#' Calculate the contribution of covariates to ability scores (beta0 + X*b)
#'
#' @param fit An object of class "irt_fit"
#' @param stan_data Optional stan_data list. If NULL, attempts to reconstruct
#'   from fit object.
#'
#' @return N x K matrix of covariate effects
#' @keywords internal
get_covariate_effects <- function(fit, stan_data = NULL) {
  # This is a simplified version - full implementation would need
  # access to the original X matrix which we should store in irt_fit

  # For now, return beta0 replicated for each person
  N <- nrow(fit$theta)
  K <- ncol(fit$theta)

  covariate_effects <- matrix(0, nrow = N, ncol = K)
  for (k in 1:K) {
    covariate_effects[, k] <- fit$beta0[k]
  }

  # If we have covariate effects in b matrix beyond intercept
  if (ncol(fit$b) > 1) {
    warning("Covariate effects beyond intercept cannot be computed without X matrix. ",
            "Returning intercept-only effects.")
  }

  return(covariate_effects)
}

#' Get Individual Effects from IRT Fit
#'
#' Extract the individual-specific deviations (zeta/eta component)
#'
#' @param fit An object of class "irt_fit"
#'
#' @return N x K matrix of individual effects
#' @keywords internal
get_individual_effects <- function(fit) {
  # Individual effects are: theta - (beta0 + X*b)
  # For intercept-only models (default), this is: theta - beta0

  N <- nrow(fit$theta)
  K <- ncol(fit$theta)

  individual_effects <- matrix(0, nrow = N, ncol = K)
  for (k in 1:K) {
    individual_effects[, k] <- fit$theta[, k] - fit$beta0[k]
  }

  return(individual_effects)
}
