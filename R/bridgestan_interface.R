# BridgeStan Interface for Alternative Optimization
#
# Functions for using BridgeStan to interface with external optimizers

#' Check if BridgeStan is Available
#'
#' @return Logical indicating if bridgestan package is installed
#' @export
has_bridgestan <- function() {
  requireNamespace("bridgestan", quietly = TRUE)
}

#' Create BridgeStan Model Object
#'
#' Creates a BridgeStan model for use with external optimizers
#'
#' @param stan_data List of data prepared by \code{\link{prepare_irt_data}}
#' @param model_path Optional path to Stan model file. If NULL, uses package model.
#'
#' @return BridgeStan model object
#'
#' @details
#' BridgeStan provides a lightweight interface to Stan models, allowing use
#' of external optimization algorithms like optim(), optimx, or nloptr.
#'
#' @examples
#' \dontrun{
#' stan_data <- prepare_irt_data(response_data, loadings, thresh_L, thresh_R)
#' bs_model <- create_bridgestan_model(stan_data)
#' }
#'
#' @export
create_bridgestan_model <- function(stan_data, model_path = NULL) {
  if (!has_bridgestan()) {
    stop("bridgestan package is required for this functionality. ",
         "Install with: install.packages('bridgestan')")
  }

  # Get model path
  if (is.null(model_path)) {
    model_path <- system.file("stan", "multidim_irt.stan", package = "IRTscoring")
  }

  if (!file.exists(model_path)) {
    stop("Stan model file not found at: ", model_path)
  }

  # Create temporary JSON file for data
  data_json <- create_stan_json(stan_data)
  data_file <- tempfile(fileext = ".json")
  writeLines(data_json, data_file)

  # Initialize BridgeStan model
  # Note: BridgeStan API may vary - try common parameter names
  tryCatch({
    bs_model <- bridgestan::StanModel$new(
      model = model_path,
      data = data_file
    )
    return(bs_model)
  }, error = function(e) {
    stop("Failed to create BridgeStan model: ", e$message)
  }, finally = {
    # Clean up temporary file
    if (file.exists(data_file)) {
      unlink(data_file)
    }
  })
}

#' Create Stan JSON from Data List
#'
#' Converts R data list to Stan JSON format
#'
#' @param stan_data Data list from prepare_irt_data()
#' @return JSON string
#' @keywords internal
create_stan_json <- function(stan_data) {
  # Convert list to JSON
  # Stan expects specific formatting for arrays and matrices
  json_parts <- character()

  for (name in names(stan_data)) {
    value <- stan_data[[name]]

    if (is.matrix(value)) {
      # Matrix: need to specify as array of arrays (row-major)
      rows <- lapply(1:nrow(value), function(i) {
        paste0("[", paste(value[i, ], collapse = ", "), "]")
      })
      json_value <- paste0("[", paste(rows, collapse = ", "), "]")
    } else if (is.vector(value) && length(value) > 1) {
      # Vector
      json_value <- paste0("[", paste(value, collapse = ", "), "]")
    } else {
      # Scalar
      json_value <- as.character(value)
    }

    json_parts <- c(json_parts, paste0('"', name, '": ', json_value))
  }

  paste0("{", paste(json_parts, collapse = ", "), "}")
}

#' Objective Function with Gradient for Optimization
#'
#' Returns negative log posterior and gradient for use with optimizers
#'
#' @param params Parameter vector
#' @param bs_model BridgeStan model object
#' @param return_gradient Logical, whether to compute and return gradient
#'
#' @return Negative log posterior (scalar) with gradient as attribute if requested
#' @keywords internal
objective_with_gradient <- function(params, bs_model, return_gradient = TRUE) {
  # Evaluate log posterior
  lp <- bs_model$log_density(params)

  # Return negative for minimization
  nlp <- -lp

  if (return_gradient) {
    # Evaluate gradient
    grad <- bs_model$log_density_gradient(params)
    attr(nlp, "gradient") <- -grad
  }

  return(nlp)
}

#' Get Initial Parameter Values
#'
#' Generates reasonable initial values for optimization
#'
#' @param stan_data Data list from prepare_irt_data()
#' @param bs_model BridgeStan model object
#'
#' @return Numeric vector of initial parameter values
#' @keywords internal
get_initial_values <- function(stan_data, bs_model) {
  # Get parameter dimensions from data
  N <- stan_data$N
  K <- stan_data$K
  J <- stan_data$J

  # Calculate total number of parameters
  # beta0[K] + b[K,J] + zeta[N,K] + L_Omega[K*(K-1)/2] + tau[K]
  n_params <- bs_model$param_num()

  # Generate initial values
  # Use small random values around zero for most parameters
  init <- rnorm(n_params, mean = 0, sd = 0.1)

  # Adjust tau (SD parameters) to be positive
  # Tau parameters are typically near the end
  tau_idx <- (n_params - K + 1):n_params
  init[tau_idx] <- abs(init[tau_idx]) + 0.5

  return(init)
}

#' Fit IRT Model Using optim() via BridgeStan
#'
#' Uses base R optim() function with BridgeStan interface for optimization
#'
#' @param stan_data List of data prepared by \code{\link{prepare_irt_data}}
#' @param method Optimization method for optim() (default: "L-BFGS-B")
#' @param init Optional initial parameter values
#' @param control List of control parameters passed to optim()
#' @param ... Additional arguments passed to optim()
#'
#' @return Object of class "irt_fit"
#'
#' @details
#' This function provides access to R's built-in optimization algorithms
#' via the BridgeStan interface. Common methods include:
#' - "L-BFGS-B": Limited-memory BFGS with box constraints
#' - "BFGS": Quasi-Newton method
#' - "CG": Conjugate gradient
#'
#' Requires the bridgestan package to be installed.
#'
#' @examples
#' \dontrun{
#' stan_data <- prepare_irt_data(response_data, loadings, thresh_L, thresh_R)
#' fit <- fit_irt_optim(stan_data, method = "L-BFGS-B")
#' }
#'
#' @export
fit_irt_optim <- function(stan_data,
                          method = "L-BFGS-B",
                          init = NULL,
                          control = list(),
                          ...) {
  if (!has_bridgestan()) {
    stop("bridgestan package is required for fit_irt_optim(). ",
         "Use fit_irt_map() instead for Stan-based optimization.")
  }

  # Create BridgeStan model
  bs_model <- create_bridgestan_model(stan_data)

  # Get initial values
  if (is.null(init)) {
    init <- get_initial_values(stan_data, bs_model)
  }

  # Set default control parameters
  default_control <- list(
    maxit = 1000,
    trace = 0,
    REPORT = 10
  )
  control <- c(control, default_control[!names(default_control) %in% names(control)])

  # Define objective function (returns scalar)
  obj_fn <- function(params) {
    objective_with_gradient(params, bs_model, return_gradient = FALSE)
  }

  # Define gradient function
  grad_fn <- function(params) {
    attr(objective_with_gradient(params, bs_model, return_gradient = TRUE), "gradient")
  }

  # Run optimization
  fit <- optim(
    par = init,
    fn = obj_fn,
    gr = grad_fn,
    method = method,
    control = control,
    ...
  )

  # Extract parameter estimates
  theta <- extract_theta_from_params(fit$par, stan_data, bs_model)
  beta0 <- extract_beta0_from_params(fit$par, stan_data, bs_model)
  b <- extract_b_from_params(fit$par, stan_data, bs_model)
  tau <- extract_tau_from_params(fit$par, stan_data, bs_model)
  Omega <- extract_omega_from_params(fit$par, stan_data, bs_model)

  # Return irt_fit object
  new_irt_fit(
    theta = theta,
    beta0 = beta0,
    b = b,
    tau = tau,
    Omega = Omega,
    log_lik = -fit$value,  # Convert back to log posterior
    method = paste0("optim_", method),
    backend = "bridgestan",
    stan_fit = fit,
    convergence = fit$convergence
  )
}

#' Fit IRT Model Using nloptr via BridgeStan
#'
#' Uses nloptr package algorithms with BridgeStan interface
#'
#' @param stan_data List of data prepared by \code{\link{prepare_irt_data}}
#' @param algorithm Algorithm for nloptr (default: "NLOPT_LD_LBFGS")
#' @param init Optional initial parameter values
#' @param opts List of options passed to nloptr
#' @param ... Additional arguments
#'
#' @return Object of class "irt_fit"
#'
#' @details
#' Provides access to advanced optimization algorithms from the NLopt library
#' via nloptr package. Common algorithms:
#' - "NLOPT_LD_LBFGS": Low-storage BFGS
#' - "NLOPT_LD_VAR1": Shifted limited-memory variable-metric
#' - "NLOPT_LD_TNEWTON_PRECOND": Truncated Newton with preconditioning
#'
#' Requires both bridgestan and nloptr packages.
#'
#' @examples
#' \dontrun{
#' fit <- fit_irt_nloptr(stan_data, algorithm = "NLOPT_LD_LBFGS")
#' }
#'
#' @export
fit_irt_nloptr <- function(stan_data,
                           algorithm = "NLOPT_LD_LBFGS",
                           init = NULL,
                           opts = list(),
                           ...) {
  if (!has_bridgestan()) {
    stop("bridgestan package is required for fit_irt_nloptr()")
  }
  if (!requireNamespace("nloptr", quietly = TRUE)) {
    stop("nloptr package is required. Install with: install.packages('nloptr')")
  }

  # Create BridgeStan model
  bs_model <- create_bridgestan_model(stan_data)

  # Get initial values
  if (is.null(init)) {
    init <- get_initial_values(stan_data, bs_model)
  }

  # Set default options
  default_opts <- list(
    algorithm = algorithm,
    maxeval = 1000,
    xtol_rel = 1e-6,
    ftol_rel = 1e-6
  )
  opts <- c(opts, default_opts[!names(default_opts) %in% names(opts)])

  # Define objective and gradient
  eval_f <- function(params) {
    objective_with_gradient(params, bs_model, return_gradient = FALSE)
  }

  eval_grad_f <- function(params) {
    attr(objective_with_gradient(params, bs_model, return_gradient = TRUE), "gradient")
  }

  # Run optimization
  fit <- nloptr::nloptr(
    x0 = init,
    eval_f = eval_f,
    eval_grad_f = eval_grad_f,
    opts = opts
  )

  # Extract parameters
  theta <- extract_theta_from_params(fit$solution, stan_data, bs_model)
  beta0 <- extract_beta0_from_params(fit$solution, stan_data, bs_model)
  b <- extract_b_from_params(fit$solution, stan_data, bs_model)
  tau <- extract_tau_from_params(fit$solution, stan_data, bs_model)
  Omega <- extract_omega_from_params(fit$solution, stan_data, bs_model)

  # Return irt_fit object
  new_irt_fit(
    theta = theta,
    beta0 = beta0,
    b = b,
    tau = tau,
    Omega = Omega,
    log_lik = -fit$objective,
    method = paste0("nloptr_", algorithm),
    backend = "bridgestan",
    stan_fit = fit,
    status = fit$status
  )
}

# ============================================================================
# Parameter Extraction Helpers
# ============================================================================

#' Extract Theta from Parameter Vector
#' @keywords internal
extract_theta_from_params <- function(params, stan_data, bs_model) {
  # Use BridgeStan to extract constrained parameters
  constrained <- bs_model$param_constrain(params)
  param_names <- bs_model$param_names()

  # Find theta parameters
  theta_idx <- grep("^theta\\[", param_names)
  theta_vals <- constrained[theta_idx]

  # Reshape to N x K matrix
  matrix(theta_vals, nrow = stan_data$N, ncol = stan_data$K, byrow = TRUE)
}

#' Extract Beta0 from Parameter Vector
#' @keywords internal
extract_beta0_from_params <- function(params, stan_data, bs_model) {
  constrained <- bs_model$param_constrain(params)
  param_names <- bs_model$param_names()

  beta0_idx <- grep("^beta0\\[", param_names)
  constrained[beta0_idx]
}

#' Extract b from Parameter Vector
#' @keywords internal
extract_b_from_params <- function(params, stan_data, bs_model) {
  constrained <- bs_model$param_constrain(params)
  param_names <- bs_model$param_names()

  b_idx <- grep("^b\\[", param_names)
  b_vals <- constrained[b_idx]

  # Reshape to K x J matrix
  matrix(b_vals, nrow = stan_data$K, ncol = stan_data$J, byrow = TRUE)
}

#' Extract Tau from Parameter Vector
#' @keywords internal
extract_tau_from_params <- function(params, stan_data, bs_model) {
  constrained <- bs_model$param_constrain(params)
  param_names <- bs_model$param_names()

  tau_idx <- grep("^tau\\[", param_names)
  constrained[tau_idx]
}

#' Extract Omega from Parameter Vector
#' @keywords internal
extract_omega_from_params <- function(params, stan_data, bs_model) {
  constrained <- bs_model$param_constrain(params)
  param_names <- bs_model$param_names()

  Omega_idx <- grep("^Omega\\[", param_names)
  Omega_vals <- constrained[Omega_idx]

  # Reshape to K x K matrix
  matrix(Omega_vals, nrow = stan_data$K, ncol = stan_data$K, byrow = TRUE)
}
