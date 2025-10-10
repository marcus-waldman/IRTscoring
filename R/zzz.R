# Package initialization and Stan model compilation
#
# This file handles package loading and Stan model compilation

# Package environment to store compiled Stan model
.irtscoring_env <- new.env(parent = emptyenv())

#' Package Load Hook
#'
#' Called when the package is loaded. Compiles the Stan model if needed.
#'
#' @param libname Library name
#' @param pkgname Package name
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Initialize package environment variables
  .irtscoring_env$stan_model_compiled <- FALSE
  .irtscoring_env$stan_model <- NULL
  .irtscoring_env$backend <- NULL

  # Note: Stan model compilation is deferred until first use
  # to avoid long package load times
}

#' Package Attach Hook
#'
#' Called when the package is attached. Provides informational messages.
#'
#' @param libname Library name
#' @param pkgname Package name
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("IRTscoring: Flexible Item Response Theory Scoring")
  packageStartupMessage("Stan models will be compiled on first use")
}

#' Get Package Environment
#'
#' Returns the package environment containing compiled Stan models
#'
#' @return Environment
#' @keywords internal
.get_irtscoring_env <- function() {
  .irtscoring_env
}
