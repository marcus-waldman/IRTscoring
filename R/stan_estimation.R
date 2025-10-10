# Stan Model Access and Compilation
#
# Functions for compiling and accessing the Stan model

#' Get Compiled Stan Model
#'
#' Returns a compiled Stan model object. Compiles the model on first call
#' and caches it for subsequent calls. Supports both rstan and cmdstanr backends.
#'
#' @param backend Character string specifying which backend to use: "cmdstanr"
#'   (default) or "rstan". If NULL, will try cmdstanr first, then rstan.
#' @param recompile Logical; if TRUE, force recompilation even if model is cached
#'
#' @return Compiled Stan model object (class depends on backend)
#'
#' @details
#' The function attempts to use cmdstanr by default as it typically provides
#' better performance and more features. If cmdstanr is not available, it falls
#' back to rstan.
#'
#' The compiled model is cached in the package environment to avoid
#' recompilation on subsequent calls.
#'
#' @examples
#' \dontrun{
#' # Get model using default backend (cmdstanr)
#' model <- get_stan_model()
#'
#' # Force use of rstan
#' model <- get_stan_model(backend = "rstan")
#'
#' # Force recompilation
#' model <- get_stan_model(recompile = TRUE)
#' }
#'
#' @export
get_stan_model <- function(backend = NULL, recompile = FALSE) {
  env <- .get_irtscoring_env()

  # Check if model is already compiled and cached
  if (!recompile && env$stan_model_compiled && !is.null(env$stan_model)) {
    return(env$stan_model)
  }

  # Determine which backend to use
  if (is.null(backend)) {
    # Auto-detect: try cmdstanr first, then rstan
    if (requireNamespace("cmdstanr", quietly = TRUE)) {
      backend <- "cmdstanr"
    } else if (requireNamespace("rstan", quietly = TRUE)) {
      backend <- "rstan"
    } else {
      stop("Neither cmdstanr nor rstan is available. Please install one of these packages.")
    }
  }

  # Validate backend choice
  backend <- match.arg(backend, choices = c("cmdstanr", "rstan"))

  # Get path to Stan model file
  model_path <- system.file("stan", "multidim_irt.stan", package = "IRTscoring")

  if (model_path == "" || !file.exists(model_path)) {
    stop("Stan model file not found. Expected at: inst/stan/multidim_irt.stan")
  }

  # Compile model based on backend
  message("Compiling Stan model using ", backend, "...")

  if (backend == "cmdstanr") {
    model <- compile_cmdstanr_model(model_path)
  } else if (backend == "rstan") {
    model <- compile_rstan_model(model_path)
  }

  # Cache the compiled model
  env$stan_model <- model
  env$stan_model_compiled <- TRUE
  env$backend <- backend

  message("Stan model compiled successfully")

  return(model)
}

#' Compile Stan Model with cmdstanr
#'
#' @param model_path Path to .stan file
#' @return Compiled cmdstanr model
#' @keywords internal
compile_cmdstanr_model <- function(model_path) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("cmdstanr package is required but not installed")
  }

  tryCatch({
    model <- cmdstanr::cmdstan_model(model_path)
    return(model)
  }, error = function(e) {
    stop("Failed to compile Stan model with cmdstanr: ", e$message)
  })
}

#' Compile Stan Model with rstan
#'
#' @param model_path Path to .stan file
#' @return Compiled rstan model
#' @keywords internal
compile_rstan_model <- function(model_path) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("rstan package is required but not installed")
  }

  tryCatch({
    model <- rstan::stan_model(file = model_path)
    return(model)
  }, error = function(e) {
    stop("Failed to compile Stan model with rstan: ", e$message)
  })
}

#' Check if Stan Model is Compiled
#'
#' @return Logical indicating whether the Stan model has been compiled
#' @export
is_stan_model_compiled <- function() {
  env <- .get_irtscoring_env()
  env$stan_model_compiled && !is.null(env$stan_model)
}

#' Get Current Stan Backend
#'
#' @return Character string: "cmdstanr", "rstan", or NULL if no model compiled
#' @export
get_stan_backend <- function() {
  env <- .get_irtscoring_env()
  if (env$stan_model_compiled) {
    return(env$backend)
  } else {
    return(NULL)
  }
}
