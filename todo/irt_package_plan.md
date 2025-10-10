# R Package Implementation Plan: IRTscoring

## Package Overview

### Package Name
**IRTscoring** - Flexible Item Response Theory Scoring in R

### Purpose
Create an R package for unidimensional and multidimensional Item Response Theory (IRT) scoring with flexible estimation methods.

### Core Functionality
- **Model**: Generalized multidimensional IRT with full factor loadings, correlated dimensions, and person covariates
- **Dimensions**: Supports K = 1 (unidimensional) through K = any positive integer (multidimensional)
- **Estimation Methods**:
  - Stan-based: MAP, Laplace Approximation, Variational Bayes (including Pathfinder), Full Bayesian MCMC
  - Alternative optimizers (via BridgeStan): `optim`, `trustOptim`, `nloptr`, etc.

### Key Features
- Response-level factor loadings (L × K matrix)
- Infinite threshold handling for boundary categories
- Person-level covariates
- Correlation structure between dimensions
- Efficient computation with vectorization
- **Wide format wrapper `fscores()`** for user-friendly data input
- **0-indexed responses** (responses must start at 0, not 1)
- **Automatic threshold validation** ensures sufficient thresholds for observed responses

---

## Quick Start: fscores() Function

The package provides a user-friendly function `fscores()` that accepts data in the typical wide format (one row per person, one column per item).

### Basic Example with 0-Indexed Responses

```r
library(IRTscoring)

# ============================================
# WIDE FORMAT DATA
# ============================================
# IMPORTANT: Item responses must be 0-indexed (start at 0, not 1)

wide_data <- data.frame(
  person_id = c(1, 2, 3, 4, 5),
  age = c(25, 30, 22, 28, 35),
  gender = c("F", "M", "F", "M", "F"),
  item_Q1 = c(0, 1, 0, 2, 1),   # 3 categories: 0, 1, 2
  item_Q2 = c(1, 2, 1, 1, 2),   # 3 categories: 0, 1, 2
  item_Q3 = c(0, 0, 1, 0, 1),   # 2 categories: 0, 1
  item_Q4 = c(2, 1, 2, 1, 0),   # 3 categories: 0, 1, 2
  item_Q5 = c(1, 1, 0, 2, 1)    # 3 categories: 0, 1, 2
)

# ============================================
# ITEM PARAMETERS (named list)
# ============================================
# For K=2 dimensions

item_params <- list(
  # Item Q1: 3 response categories (0, 1, 2)
  # Category 0: (-Inf, -1.0]
  # Category 1: (-1.0, 0.5]  
  # Category 2: (0.5, Inf)
  item_Q1 = list(
    loadings = c(1.2, 0.1),      # Vector of length K=2
    thresholds = c(-1.0, 0.5)    # 2 thresholds for 3 categories
  ),
  
  # Item Q2: 3 response categories (0, 1, 2)
  item_Q2 = list(
    loadings = c(1.5, 0.2),
    thresholds = c(-0.5, 1.0)
  ),
  
  # Item Q3: 2 response categories (0, 1)
  # Category 0: (-Inf, 1.0]
  # Category 1: (1.0, Inf)
  item_Q3 = list(
    loadings = c(0.1, 1.3),      # Loads primarily on dimension 2
    thresholds = c(1.0)          # 1 threshold for 2 categories (scalar)
  ),
  
  # Item Q4: 3 response categories (0, 1, 2)
  item_Q4 = list(
    loadings = c(0.2, 1.6),
    thresholds = c(-0.8, 0.8)
  ),
  
  # Item Q5: 3 response categories (0, 1, 2)
  item_Q5 = list(
    loadings = c(1.0, 1.0),      # Loads equally on both dimensions
    thresholds = c(-0.3, 0.3)
  )
)

# ============================================
# FIT THE MODEL
# ============================================

fit <- fscores(
  data = wide_data,
  item_params = item_params,
  person_id_col = "person_id",
  covariate_cols = c("age", "gender"),
  method = "map"
)

# Get full ability scores (N x K matrix)
abilities <- predict(fit, type = "abilities")
print(abilities)

# Decompose scores
individual_effects <- predict(fit, type = "individual_effects")
covariate_effects <- predict(fit, type = "covariate_effects")

# Get dimension correlation
omega <- get_dimension_correlations(fit)
print(omega)
```

### Unidimensional Example (K=1)

```r
# Unidimensional case
item_params_1d <- list(
  item_Q1 = list(
    loadings = c(1.2),           # Just one dimension (scalar wrapped in vector)
    thresholds = c(-1.0, 0.5)    # 2 thresholds for 3 categories (0, 1, 2)
  ),
  item_Q2 = list(
    loadings = c(1.5),
    thresholds = c(-0.5, 1.0)
  ),
  item_Q3 = list(
    loadings = c(1.3),
    thresholds = c(1.0)          # Binary item: categories 0 and 1
  )
)

# Responses must be 0-indexed
wide_data_1d <- data.frame(
  person_id = 1:5,
  item_Q1 = c(0, 1, 0, 2, 1),
  item_Q2 = c(1, 2, 1, 1, 2),
  item_Q3 = c(0, 0, 1, 0, 1)
)

fit_1d <- fscores(
  data = wide_data_1d,
  item_params = item_params_1d,
  person_id_col = "person_id",
  method = "map"
)

abilities_1d <- predict(fit_1d, type = "abilities")
print(abilities_1d)  # N x 1 matrix
```

---

## Project Architecture

### Package Structure
```
IRTscoring/
├── .gitignore
├── .Rbuildignore
├── DESCRIPTION
├── NAMESPACE
├── README.md
├── NEWS.md
├── LICENSE
├── R/
│   ├── fscores.R              # Wide format user-facing function
│   ├── fit_ability.R          # Main long format function
│   ├── prepare_data.R         # Data validation and formatting
│   ├── stan_estimation.R      # Stan MAP/VB/MCMC methods
│   ├── bridgestan_estimation.R # Alternative optimizer methods
│   ├── methods.R              # S3 methods: print, summary, coef, predict
│   ├── utilities.R            # Helper functions
│   └── zzz.R                  # Package initialization
├── inst/
│   └── stan/
│       └── multidim_irt.stan  # Stan model file
├── src/
│   └── (compiled Stan code)
├── tests/
│   └── testthat/
│       ├── test-data-preparation.R
│       ├── test-stan-estimation.R
│       ├── test-bridgestan.R
│       ├── test-mirt-comparison.R
│       └── test-methods.R
├── vignettes/
│   ├── quickstart.Rmd
│   ├── unidimensional-irt.Rmd
│   ├── multidimensional-irt.Rmd
│   └── simulation-study.Rmd
└── man/                       # Documentation
```

---

## Implementation Phases

### Phase 0: Project Setup (Foundation)

#### Tasks
1. **Initialize R package structure**
   ```r
   usethis::create_package("IRTscoring")
   usethis::use_git()
   usethis::use_mit_license()
   usethis::use_roxygen_md()
   ```

2. **Create .gitignore file**
   ```
   # .gitignore
   
   # History files
   .Rhistory
   .Rapp.history
   
   # Session Data files
   .RData
   .RDataTmp
   
   # User-specific files
   .Ruserdata
   
   # RStudio files
   .Rproj.user/
   *.Rproj
   
   # OAuth2 token
   .httr-oauth
   
   # knitr and R markdown default cache directories
   *_cache/
   /cache/
   
   # Temporary files created by R markdown
   *.utf8.md
   *.knit.md
   
   # R Environment Variables
   .Renviron
   
   # pkgdown site
   docs/
   
   # translation temp files
   po/*~
   
   # Compiled Stan models
   src/*.o
   src/*.so
   src/*.dll
   
   # Mac files
   .DS_Store
   
   # Test outputs
   tests/testthat/*.rds
   tests/testthat/*.png
   
   # Simulation results
   simulation_results.rds
   ```

3. **Add package dependencies to DESCRIPTION**
   ```
   Package: IRTscoring
   Title: Flexible Item Response Theory Scoring
   Version: 0.1.0
   Authors@R: 
       person("First", "Last", email = "first.last@example.com", 
              role = c("aut", "cre"))
   Description: Provides flexible estimation methods for unidimensional and 
       multidimensional Item Response Theory (IRT) models. Supports MAP, 
       Variational Bayes (including Pathfinder), and full Bayesian MCMC 
       estimation via Stan, as well as alternative optimizers via BridgeStan.
   License: MIT + file LICENSE
   Encoding: UTF-8
   Roxygen: list(markdown = TRUE)
   RoxygenNote: 7.0.0
   Imports:
     rstan (>= 2.26.0),
     cmdstanr,
     checkmate,
     Matrix,
     posterior,
     methods
   Suggests:
     bridgestan,
     optimx,
     trustOptim,
     nloptr,
     mirt,
     MASS,
     testthat (>= 3.0.0),
     knitr,
     rmarkdown,
     ggplot2,
     dplyr
   VignetteBuilder: knitr
   ```

4. **Set up Stan infrastructure**
   ```r
   usethis::use_directory("inst/stan")
   # Copy multidim_irt.stan to inst/stan/
   # Configure package to compile Stan models on install
   ```

5. **Initialize testing framework**
   ```r
   usethis::use_testthat()
   usethis::use_test("data-preparation")
   ```

6. **Create README.md skeleton**
   ```r
   usethis::use_readme_md()
   ```

#### Validation
- Package loads without errors: `devtools::load_all()`
- Basic structure is correct: `devtools::check()`
- Git repository is initialized and .gitignore is in place

#### Phase Completion

**Git Commit:**
```bash
git add .
git commit -m "Phase 0 Complete: Project structure initialized

- Created R package skeleton for IRTscoring
- Added .gitignore for R package development
- Configured DESCRIPTION with dependencies
- Set up testing framework with testthat
- Created directory structure for Stan models
- Initialized README
"
```

**Load Next Phase:**
```
📋 PHASE 0 COMPLETE ✓

Next Steps: Begin Phase 1 - Stan Model Integration

Phase 1 will focus on:
- Creating the Stan model file (multidim_irt.stan)
- Setting up Stan compilation infrastructure
- Creating model wrapper functions
- Testing Stan model compilation

Ready to proceed with Phase 1? [Y/n]
```

---

### Phase 1: Stan Model Integration

#### Tasks

1. **Create `inst/stan/multidim_irt.stan`**
   - Copy the finalized Stan model code
   - Ensure all variable names are clear
   - Add comprehensive comments

2. **Set up Stan compilation**
   - Configure `zzz.R` for model compilation on package load
   - Handle both `rstan` and `cmdstanr` backends
   
   ```r
   # R/zzz.R
   .onLoad <- function(libname, pkgname) {
     # Compile Stan model
     model_path <- system.file("stan", "multidim_irt.stan", 
                               package = "irtscoring")
     # Store compiled model in package environment
   }
   ```

3. **Create Stan model wrapper function**
   ```r
   # R/stan_estimation.R
   get_stan_model <- function() {
     # Return compiled Stan model
     # Handle both rstan and cmdstanr
   }
   ```

#### Testing During Development
```r
# tests/testthat/test-stan-model.R
test_that("Stan model compiles successfully", {
  expect_no_error(get_stan_model())
})

test_that("Stan model has correct parameters", {
  model <- get_stan_model()
  # Check parameter names match expectations
})
```

#### Phase Completion

**Git Commit:**
```bash
git add .
git commit -m "Phase 1 Complete: Stan model integration

- Added multidim_irt.stan model file
- Configured Stan compilation in zzz.R
- Created get_stan_model() wrapper function
- Added tests for Stan model compilation
- Verified model compiles successfully
"
```

**Load Next Phase:**
```
📋 PHASE 1 COMPLETE ✓

Next Steps: Begin Phase 2 - Data Preparation Functions

Phase 2 will focus on:
- Creating prepare_irt_data() function
- Adding data validation helpers
- Converting infinite thresholds
- Formatting data for Stan
- Comprehensive input validation tests

Ready to proceed with Phase 2? [Y/n]
```

---

### Phase 2: Data Preparation Functions

#### Tasks

1. **Create `prepare_data.R` with validation functions**

```r
#' Prepare IRT data for estimation
#' 
#' @param response_data Data frame with columns: pid, iid, response
#' @param item_loadings Matrix L x K of factor loadings
#' @param threshold_left Vector of left thresholds (use Inf for -infinity)
#' @param threshold_right Vector of right thresholds (use Inf for +infinity)
#' @param person_covariates Matrix N x J of person-level covariates (optional)
#' @param weights Vector of person weights (optional)
#' @param K Number of dimensions (auto-detected from item_loadings if NULL)
#' @return List formatted for Stan
prepare_irt_data <- function(response_data, 
                             item_loadings,
                             threshold_left,
                             threshold_right,
                             person_covariates = NULL,
                             weights = NULL,
                             K = NULL) {
  
  # Validate inputs
  validate_response_data(response_data)
  validate_item_loadings(item_loadings, K)
  validate_thresholds(threshold_left, threshold_right, 
                      nrow(response_data))
  
  # Auto-detect K if not provided
  if(is.null(K)) K <- ncol(item_loadings)
  
  # Handle infinite thresholds
  inf_code <- -999
  threshold_left[is.infinite(threshold_left)] <- inf_code
  threshold_right[is.infinite(threshold_right)] <- inf_code
  
  # Create person covariate matrix
  if(is.null(person_covariates)) {
    person_covariates <- matrix(0, nrow = length(unique(response_data$pid)), 
                                ncol = 1)
  }
  
  # Create weights
  if(is.null(weights)) {
    weights <- rep(1, length(unique(response_data$pid)))
  }
  
  # Format for Stan
  list(
    inf = inf_code,
    L = nrow(response_data),
    N = length(unique(response_data$pid)),
    J = ncol(person_covariates),
    K = K,
    pid = response_data$pid,
    dL = threshold_left,
    dR = threshold_right,
    wgt = weights,
    X = person_covariates,
    A = item_loadings
  )
}

# Validation helper functions
validate_response_data <- function(data) {
  checkmate::assert_data_frame(data)
  checkmate::assert_names(colnames(data), must.include = c("pid", "iid"))
  # Additional checks...
}

validate_item_loadings <- function(loadings, K) {
  checkmate::assert_matrix(loadings)
  if(!is.null(K)) {
    checkmate::assert_true(ncol(loadings) == K)
  }
}

validate_thresholds <- function(left, right, L) {
  checkmate::assert_numeric(left, len = L)
  checkmate::assert_numeric(right, len = L)
  # Check that left <= right (except for infinities)
}
```

#### Testing During Development
```r
# tests/testthat/test-data-preparation.R
test_that("Data preparation handles simple case", {
  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )
  loadings <- matrix(c(1, 1), ncol = 1)
  thresh_L <- c(-Inf, -1)
  thresh_R <- c(-1, Inf)
  
  result <- prepare_irt_data(response_data, loadings, 
                             thresh_L, thresh_R)
  
  expect_equal(result$K, 1)
  expect_equal(result$L, 4)
  expect_equal(result$N, 2)
  expect_true(all(result$dL[is.infinite(thresh_L)] == -999))
})

test_that("Data preparation handles multidimensional case", {
  response_data <- data.frame(
    pid = c(1, 1, 2, 2),
    iid = c(1, 2, 1, 2)
  )
  loadings <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
  thresh_L <- c(-Inf, -1, -Inf, -1)
  thresh_R <- c(-1, Inf, -1, Inf)
  
  result <- prepare_irt_data(response_data, loadings, 
                             thresh_L, thresh_R)
  
  expect_equal(result$K, 2)
  expect_equal(ncol(result$A), 2)
})

test_that("Data preparation validates inputs", {
  bad_data <- data.frame(person = 1:10)  # Missing required columns
  loadings <- matrix(1, ncol = 1)
  
  expect_error(
    prepare_irt_data(bad_data, loadings, 1, 2),
    "must.include"
  )
})
```

#### Phase Completion

**Git Commit:**
```bash
git add .
git commit -m "Phase 2 Complete: Data preparation functions

- Implemented prepare_irt_data() with full validation
- Added validate_response_data() helper
- Added validate_item_loadings() helper
- Added validate_thresholds() helper
- Implemented infinite threshold handling
- Created comprehensive validation tests
- All data preparation tests passing
"
```

**Load Next Phase:**
```
📋 PHASE 2 COMPLETE ✓

Next Steps: Begin Phase 3 - Stan-based Estimation Functions

Phase 3 will focus on:
- Implementing fit_irt_map() for MAP estimation
- Implementing fit_irt_vb() for Variational Bayes
- Implementing fit_irt_mcmc() for full Bayesian
- Supporting Pathfinder algorithm
- Creating irt_fit S3 class
- Testing all Stan estimation methods

Ready to proceed with Phase 3? [Y/n]
```

---

### Phase 3: Stan-based Estimation Functions

#### Tasks

1. **Implement MAP estimation**

```r
#' Estimate IRT model using MAP (Maximum A Posteriori)
#' 
#' @param stan_data List prepared by prepare_irt_data()
#' @param algorithm Optimization algorithm ("LBFGS" or "Newton")
#' @param ... Additional arguments passed to rstan::optimizing()
#' @return irt_fit object
fit_irt_map <- function(stan_data, algorithm = "LBFGS", ...) {
  model <- get_stan_model()
  
  fit <- rstan::optimizing(
    model,
    data = stan_data,
    algorithm = algorithm,
    ...
  )
  
  # Extract and format results
  structure(
    list(
      method = "MAP",
      estimates = fit$par,
      value = fit$value,
      return_code = fit$return_code,
      stan_data = stan_data,
      raw_fit = fit
    ),
    class = "irt_fit"
  )
}
```

2. **Implement Variational Bayes estimation**

```r
#' Estimate IRT model using Variational Bayes
#' 
#' @param stan_data List prepared by prepare_irt_data()
#' @param algorithm "meanfield" or "fullrank"
#' @param use_pathfinder Logical, use Pathfinder algorithm if available
#' @param ... Additional arguments
#' @return irt_fit object
fit_irt_vb <- function(stan_data, 
                       algorithm = "meanfield",
                       use_pathfinder = TRUE,
                       ...) {
  model <- get_stan_model()
  
  if(use_pathfinder && has_pathfinder()) {
    fit <- run_pathfinder(model, stan_data, ...)
  } else {
    fit <- rstan::vb(
      model,
      data = stan_data,
      algorithm = algorithm,
      ...
    )
  }
  
  # Extract and format results
  structure(
    list(
      method = ifelse(use_pathfinder, "Pathfinder", "VB"),
      samples = extract_vb_samples(fit),
      stan_data = stan_data,
      raw_fit = fit
    ),
    class = "irt_fit"
  )
}
```

3. **Implement Full Bayesian (MCMC) estimation**

```r
#' Estimate IRT model using Full Bayesian MCMC
#' 
#' @param stan_data List prepared by prepare_irt_data()
#' @param chains Number of MCMC chains
#' @param iter Number of iterations per chain
#' @param warmup Number of warmup iterations
#' @param ... Additional arguments passed to rstan::sampling()
#' @return irt_fit object
fit_irt_mcmc <- function(stan_data,
                         chains = 4,
                         iter = 2000,
                         warmup = 1000,
                         ...) {
  model <- get_stan_model()
  
  fit <- rstan::sampling(
    model,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    ...
  )
  
  # Extract and format results
  structure(
    list(
      method = "MCMC",
      samples = rstan::extract(fit),
      summary = rstan::summary(fit)$summary,
      stan_data = stan_data,
      raw_fit = fit
    ),
    class = "irt_fit"
  )
}
```

4. **Implement Laplace Approximation**

```r
#' Estimate IRT model using Laplace Approximation
#' 
#' Uses MAP estimate and computes Hessian-based approximation to posterior
#' 
#' @param stan_data List prepared by prepare_irt_data()
#' @param algorithm Optimization algorithm for MAP step
#' @param draws Number of draws from approximated posterior (default = 1000)
#' @param ... Additional arguments passed to rstan::optimizing()
#' @return irt_fit object
fit_irt_laplace <- function(stan_data,
                            algorithm = "LBFGS",
                            draws = 1000,
                            ...) {
  model <- get_stan_model()
  
  # First, get MAP estimate with Hessian
  fit_map <- rstan::optimizing(
    model,
    data = stan_data,
    algorithm = algorithm,
    hessian = TRUE,
    ...
  )
  
  if(fit_map$return_code != 0) {
    warning("MAP optimization did not converge (return code: ", 
            fit_map$return_code, ")")
  }
  
  # Extract mode and Hessian
  mode <- fit_map$par
  hessian <- fit_map$hessian
  
  # Compute covariance matrix (inverse negative Hessian)
  # Sigma = -H^(-1)
  tryCatch({
    cov_matrix <- solve(-hessian)
    
    # Check if covariance matrix is positive definite
    eigenvalues <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
    if(any(eigenvalues <= 0)) {
      warning("Covariance matrix from Laplace approximation is not positive definite. ",
              "Consider using a different estimation method.")
    }
    
    # Draw samples from multivariate normal approximation
    samples <- MASS::mvrnorm(n = draws, mu = mode, Sigma = cov_matrix)
    
    # Parse samples into parameter arrays
    parsed_samples <- parse_parameter_samples(samples, stan_data)
    
    structure(
      list(
        method = "Laplace",
        mode = mode,
        covariance = cov_matrix,
        samples = parsed_samples,
        stan_data = stan_data,
        raw_fit = fit_map,
        hessian = hessian
      ),
      class = "irt_fit"
    )
    
  }, error = function(e) {
    warning("Failed to compute Laplace approximation: ", e$message,
            "\nReturning MAP estimate only.")
    
    structure(
      list(
        method = "Laplace-Failed",
        estimates = mode,
        value = fit_map$value,
        return_code = fit_map$return_code,
        stan_data = stan_data,
        raw_fit = fit_map,
        error = e$message
      ),
      class = "irt_fit"
    )
  })
}

#' Parse parameter vector into named arrays
#' @keywords internal
parse_parameter_samples <- function(samples, stan_data) {
  # Get parameter dimensions
  N <- stan_data$N
  K <- stan_data$K
  J <- stan_data$J
  
  # Initialize lists for each parameter
  n_samples <- nrow(samples)
  
  # Assuming parameter order matches Stan model declaration:
  # beta0[K], b[K,J], zeta[N,K], L_Omega (lower triangle), tau[K]
  
  idx <- 1
  
  # beta0
  beta0 <- samples[, idx:(idx + K - 1), drop = FALSE]
  idx <- idx + K
  
  # b (K x J)
  b <- array(0, dim = c(n_samples, K, J))
  for(k in 1:K) {
    b[, k, ] <- samples[, idx:(idx + J - 1), drop = FALSE]
    idx <- idx + J
  }
  
  # zeta (N x K) 
  zeta <- array(0, dim = c(n_samples, N, K))
  for(n in 1:N) {
    zeta[, n, ] <- samples[, idx:(idx + K - 1), drop = FALSE]
    idx <- idx + K
  }
  
  # L_Omega (Cholesky factor, lower triangular)
  n_omega_params <- K * (K - 1) / 2
  if(K > 1) {
    L_Omega_vec <- samples[, idx:(idx + n_omega_params - 1), drop = FALSE]
    idx <- idx + n_omega_params
    # Would need to reconstruct full Cholesky factor from vector
  }
  
  # tau
  tau <- samples[, idx:(idx + K - 1), drop = FALSE]
  
  list(
    beta0 = beta0,
    b = b,
    zeta = zeta,
    tau = tau
  )
}
```

#### Testing During Development
```r
# tests/testthat/test-stan-estimation.R
test_that("MAP estimation runs on simple data", {
  # Create minimal test data
  stan_data <- create_simple_test_data(K = 1)
  
  fit <- fit_irt_map(stan_data, iter = 100)
  
  expect_s3_class(fit, "irt_fit")
  expect_equal(fit$method, "MAP")
  expect_true(!is.null(fit$estimates))
  expect_equal(fit$return_code, 0)
})

test_that("VB estimation runs on simple data", {
  stan_data <- create_simple_test_data(K = 1)
  
  fit <- fit_irt_vb(stan_data, iter = 100, use_pathfinder = FALSE)
  
  expect_s3_class(fit, "irt_fit")
  expect_true(!is.null(fit$samples))
})

test_that("MCMC estimation runs on simple data", {
  skip_on_cran()  # Too slow for CRAN
  
  stan_data <- create_simple_test_data(K = 1)
  
  fit <- fit_irt_mcmc(stan_data, chains = 1, iter = 200, warmup = 100)
  
  expect_s3_class(fit, "irt_fit")
  expect_equal(fit$method, "MCMC")
})

test_that("Laplace approximation runs on simple data", {
  stan_data <- create_simple_test_data(K = 1)
  
  fit <- fit_irt_laplace(stan_data, draws = 100)
  
  expect_s3_class(fit, "irt_fit")
  expect_true(fit$method %in% c("Laplace", "Laplace-Failed"))
  
  if(fit$method == "Laplace") {
    expect_true(!is.null(fit$samples))
    expect_true(!is.null(fit$covariance))
  }
})
```

#### Phase Completion

**Git Commit:**
```bash
git add .
git commit -m "Phase 3 Complete: Stan-based estimation methods

- Implemented fit_irt_map() for MAP estimation
- Implemented fit_irt_vb() with Pathfinder support
- Implemented fit_irt_mcmc() for full Bayesian
- Created irt_fit S3 class structure
- Added helper functions for VB and Pathfinder
- All Stan estimation tests passing
- Verified MAP, VB, and MCMC work correctly
"
```

**Load Next Phase:**
```
📋 PHASE 3 COMPLETE ✓

Next Steps: Begin Phase 4 - BridgeStan Integration

Phase 4 will focus on:
- Creating BridgeStan interface functions
- Implementing has_bridgestan() checker
- Creating objective_with_gradient() function
- Implementing fit_irt_optim() wrapper
- Implementing fit_irt_nloptr() wrapper
- Supporting multiple optimization algorithms
- Testing with and without BridgeStan installed

Ready to proceed with Phase 4? [Y/n]
```

---

### Phase 4: BridgeStan Integration

#### Tasks

1. **Create BridgeStan interface**

```r
#' Check if BridgeStan is available
#' @return Logical
has_bridgestan <- function() {
  requireNamespace("bridgestan", quietly = TRUE)
}

#' Create BridgeStan model object
#' @param stan_data Prepared data list
#' @return BridgeStan model object
create_bridgestan_model <- function(stan_data) {
  if(!has_bridgestan()) {
    stop("bridgestan package is required for this functionality")
  }
  
  model_path <- system.file("stan", "multidim_irt.stan", 
                           package = "irtscoring")
  
  # Create data JSON
  data_json <- create_stan_json(stan_data)
  
  # Initialize BridgeStan model
  bridgestan::StanModel$new(
    model = model_path,
    data = data_json
  )
}

#' Objective function for optimization
#' Returns negative log posterior and gradient
objective_with_gradient <- function(params, bs_model) {
  # Evaluate log posterior
  lp <- bs_model$log_density(params)
  
  # Evaluate gradient
  grad <- bs_model$log_density_gradient(params)
  
  # Return negative (for minimization) with gradient attribute
  nlp <- -lp
  attr(nlp, "gradient") <- -grad
  
  return(nlp)
}
```

2. **Implement optimization wrappers**

```r
#' Estimate IRT using optim() via BridgeStan
#' 
#' @param stan_data Prepared data list
#' @param method Optimization method for optim()
#' @param init Initial parameter values (optional)
#' @param ... Additional arguments passed to optim()
#' @return irt_fit object
fit_irt_optim <- function(stan_data, 
                          method = "L-BFGS-B",
                          init = NULL,
                          ...) {
  if(!has_bridgestan()) {
    stop("bridgestan package required for this method")
  }
  
  bs_model <- create_bridgestan_model(stan_data)
  
  # Get initial values
  if(is.null(init)) {
    init <- get_initial_values(stan_data)
  }
  
  # Run optimization
  fit <- optim(
    par = init,
    fn = objective_with_gradient,
    gr = function(x) attr(objective_with_gradient(x, bs_model), "gradient"),
    method = method,
    bs_model = bs_model,
    ...
  )
  
  structure(
    list(
      method = paste0("optim-", method),
      estimates = fit$par,
      value = -fit$value,  # Convert back to log posterior
      convergence = fit$convergence,
      stan_data = stan_data,
      raw_fit = fit
    ),
    class = "irt_fit"
  )
}

#' Estimate IRT using nloptr via BridgeStan
fit_irt_nloptr <- function(stan_data, 
                           algorithm = "NLOPT_LD_LBFGS",
                           init = NULL,
                           ...) {
  if(!has_bridgestan()) {
    stop("bridgestan package required")
  }
  if(!requireNamespace("nloptr", quietly = TRUE)) {
    stop("nloptr package required")
  }
  
  bs_model <- create_bridgestan_model(stan_data)
  
  if(is.null(init)) {
    init <- get_initial_values(stan_data)
  }
  
  fit <- nloptr::nloptr(
    x0 = init,
    eval_f = objective_with_gradient,
    eval_grad_f = function(x) {
      attr(objective_with_gradient(x, bs_model), "gradient")
    },
    opts = list(algorithm = algorithm, ...),
    bs_model = bs_model
  )
  
  structure(
    list(
      method = paste0("nloptr-", algorithm),
      estimates = fit$solution,
      value = -fit$objective,
      status = fit$status,
      stan_data = stan_data,
      raw_fit = fit
    ),
    class = "irt_fit"
  )
}
```

#### Testing During Development
```r
# tests/testthat/test-bridgestan.R
test_that("BridgeStan model creation works", {
  skip_if_not_installed("bridgestan")
  
  stan_data <- create_simple_test_data(K = 1)
  
  expect_no_error(create_bridgestan_model(stan_data))
})

test_that("optim via BridgeStan works", {
  skip_if_not_installed("bridgestan")
  
  stan_data <- create_simple_test_data(K = 1)
  
  fit <- fit_irt_optim(stan_data, control = list(maxit = 50))
  
  expect_s3_class(fit, "irt_fit")
  expect_true(!is.null(fit$estimates))
  expect_equal(fit$convergence, 0)
})

test_that("nloptr via BridgeStan works", {
  skip_if_not_installed("bridgestan")
  skip_if_not_installed("nloptr")
  
  stan_data <- create_simple_test_data(K = 1)
  
  fit <- fit_irt_nloptr(stan_data, maxeval = 50)
  
  expect_s3_class(fit, "irt_fit")
  expect_true(!is.null(fit$estimates))
})
```

#### Phase Completion

**Git Commit:**
```bash
git add .
git commit -m "Phase 4 Complete: BridgeStan integration

- Implemented has_bridgestan() checker
- Created create_bridgestan_model() interface
- Implemented objective_with_gradient() function
- Added fit_irt_optim() for optim() integration
- Added fit_irt_nloptr() for nloptr integration
- Created get_initial_values() helper
- All BridgeStan tests passing
- Verified gradient-based optimization works
"
```

**Load Next Phase:**
```
📋 PHASE 4 COMPLETE ✓

Next Steps: Begin Phase 5 - Main User Interface

Phase 5 will focus on:
- Creating unified fit_irt() function
- Implementing method routing logic
- Adding fit_unidim_irt() convenience wrapper
- Comprehensive argument handling
- User-friendly error messages
- Testing main interface with all methods

Ready to proceed with Phase 5? [Y/n]
```

---

### Phase 5: Main User Interface

#### Tasks

1. **Create unified `fit_ability()` function**

```r
#' Fit IRT Model for Ability Scoring
#' 
#' Main interface for fitting unidimensional or multidimensional IRT models
#' 
#' @param response_data Data frame with person and item IDs
#' @param item_loadings Matrix of factor loadings (L x K)
#' @param threshold_left Vector of left thresholds
#' @param threshold_right Vector of right thresholds
#' @param person_covariates Optional matrix of person covariates
#' @param weights Optional person weights
#' @param method Estimation method: "map", "vb", "mcmc", "optim", "nloptr", etc.
#' @param ... Additional arguments passed to specific estimation function
#' 
#' @return Object of class 'irt_fit'
#' 
#' @examples
#' # Unidimensional IRT with MAP estimation
#' fit <- fit_ability(
#'   response_data = responses,
#'   item_loadings = matrix(c(1, 1, 1), ncol = 1),
#'   threshold_left = c(-Inf, -1, 0),
#'   threshold_right = c(-1, 0, Inf),
#'   method = "map"
#' )
#' 
#' # With Laplace approximation (MAP + uncertainty)
#' fit <- fit_ability(
#'   response_data = responses,
#'   item_loadings = matrix(c(1, 1, 1), ncol = 1),
#'   threshold_left = c(-Inf, -1, 0),
#'   threshold_right = c(-1, 0, Inf),
#'   method = "laplace",
#'   draws = 1000
#' )
#' 
#' # Multidimensional IRT with MCMC
#' fit <- fit_ability(
#'   response_data = responses,
#'   item_loadings = loading_matrix,  # L x 2 matrix
#'   threshold_left = thresh_l,
#'   threshold_right = thresh_r,
#'   method = "mcmc",
#'   chains = 4,
#'   iter = 2000
#' )
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
  
  method <- match.arg(method)
  
  # Prepare data
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
    stop("Unknown method: ", method)
  )
  
  return(fit)
}
```

2. **Add wide format wrapper function `fscores()`**

```r
#' Calculate Factor Scores Using IRT Model
#' 
#' User-friendly function that accepts data in wide format (one row per person)
#' and returns IRT ability scores. Similar interface to mirt::fscores().
#' 
#' @param data Data frame in wide format with one row per person
#' @param item_params Named list of item parameters. Each element should contain:
#'   \itemize{
#'     \item \code{loadings}: numeric vector of length K (number of dimensions)
#'     \item \code{thresholds}: numeric vector of thresholds (length = n_categories - 1)
#'   }
#' @param person_id_col Name of the column containing person IDs
#' @param covariate_cols Optional vector of column names to use as covariates
#' @param method Estimation method: "map", "laplace", "vb", "mcmc", "optim", "nloptr"
#' @param ... Additional arguments passed to fit_ability()
#' 
#' @return Object of class 'irt_fit'
#' 
#' @details 
#' **IMPORTANT**: Item responses must be 0-indexed (i.e., categories start at 0, not 1).
#' For an item with 3 categories, valid responses are 0, 1, 2.
#' 
#' The thresholds define boundaries between categories:
#' - Category 0: (-Inf, threshold[1]]
#' - Category 1: (threshold[1], threshold[2]]
#' - Category 2: (threshold[2], Inf)
#' 
#' **Threshold Validation**: The function validates that sufficient thresholds are provided
#' for all observed responses. For a response value of k (0-indexed), at least k thresholds
#' must be specified. For example:
#' - Response = 0 requires 0 thresholds minimum
#' - Response = 1 requires 1 threshold minimum  
#' - Response = 2 requires 2 thresholds minimum
#' 
#' If an item has a maximum observed response of 2, you must provide at least
#' \code{thresholds = c(t1, t2)} (2 thresholds) to define the 3 categories (0, 1, 2).
#' 
#' @examples
#' # Create wide format data (responses are 0-indexed!)
#' wide_data <- data.frame(
#'   person_id = 1:100,
#'   age = rnorm(100, 30, 10),
#'   item1 = sample(0:2, 100, replace = TRUE),
#'   item2 = sample(0:2, 100, replace = TRUE)
#' )
#' 
#' # Define item parameters
#' item_params <- list(
#'   item1 = list(loadings = c(1.2, 0.1), thresholds = c(-1, 0.5)),
#'   item2 = list(loadings = c(0.8, 1.5), thresholds = c(-0.5, 1.0))
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
#' # Get ability scores
#' abilities <- predict(fit, type = "abilities")
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
  if(!person_id_col %in% colnames(data)) {
    stop("person_id_col '", person_id_col, "' not found in data")
  }
  
  # 1. Identify item columns
  all_cols <- colnames(data)
  item_cols <- setdiff(all_cols, c(person_id_col, covariate_cols))
  
  # Check that all item_cols have parameters
  missing_params <- setdiff(item_cols, names(item_params))
  if(length(missing_params) > 0) {
    stop("Missing item parameters for: ", paste(missing_params, collapse = ", "))
  }
  
  # Get number of dimensions from first item
  K <- length(item_params[[item_cols[1]]]$loadings)
  
  # 2. Convert wide to long format
  response_list <- list()
  loading_list <- list()
  thresh_left_list <- list()
  thresh_right_list <- list()
  
  for(i in 1:nrow(data)) {
    person_id <- data[[person_id_col]][i]
    
    for(item_name in item_cols) {
      response_value <- data[[item_name]][i]
      
      if(!is.na(response_value)) {
        # Validate response is 0-indexed
        if(response_value < 0) {
          stop("Response values must be >= 0 (0-indexed). Found ", 
               response_value, " for person ", person_id, ", item ", item_name)
        }
        
        # Get item parameters
        item_info <- item_params[[item_name]]
        
        # Validate sufficient thresholds for observed response
        # For response value k (0-indexed), we need at least k thresholds
        # Categories: 0, 1, ..., k means (k+1) categories, needing k thresholds
        n_thresholds_needed <- response_value
        
        if(length(item_info$thresholds) < n_thresholds_needed) {
          stop("Item '", item_name, "' has response value ", response_value,
               " (category ", response_value, "), which requires at least ", 
               n_thresholds_needed, " threshold(s), but only ", 
               length(item_info$thresholds), " threshold(s) were provided. ",
               "Response categories 0 through ", response_value, 
               " require ", response_value, " threshold(s).")
        }
        
        # Create threshold boundaries: (-Inf, t1, t2, ..., tn, Inf)
        thresholds <- c(-Inf, item_info$thresholds, Inf)
        
        # For 0-indexed response k, category spans (thresholds[k+1], thresholds[k+2]]
        # (adding 1 because R is 1-indexed)
        
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
  
  # 3. Combine into matrices
  response_data <- do.call(rbind, response_list)
  item_loadings <- do.call(rbind, loading_list)
  threshold_left <- unlist(thresh_left_list)
  threshold_right <- unlist(thresh_right_list)
  
  # 4. Prepare covariates
  if(!is.null(covariate_cols)) {
    # Create design matrix with intercept
    covariate_data <- data[, c(person_id_col, covariate_cols), drop = FALSE]
    
    # Create model matrix (handles factors automatically)
    formula_str <- paste("~", paste(covariate_cols, collapse = " + "))
    person_covariates <- model.matrix(as.formula(formula_str), 
                                      data = covariate_data)
  } else {
    # Just intercept
    N <- length(unique(response_data$pid))
    person_covariates <- matrix(1, nrow = N, ncol = 1)
  }
  
  # 5. Prepare weights
  if(is.null(weights)) {
    weights <- rep(1, length(unique(response_data$pid)))
  }
  
  # 6. Call the main function
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
```

3. **Add helper function for unidimensional case**

```r
#' Fit Unidimensional IRT Model
#' 
#' Convenience wrapper for unidimensional IRT
#' 
#' @param response_data Data frame with person and item IDs
#' @param item_discrimination Vector of discrimination parameters (length L)
#' @param threshold_left Vector of left thresholds
#' @param threshold_right Vector of right thresholds
#' @param ... Other arguments passed to fit_ability()
#' 
#' @export
fit_unidim_ability <- function(response_data,
                           item_discrimination,
                           threshold_left,
                           threshold_right,
                           ...) {
  
  # Convert to L x 1 loading matrix
  item_loadings <- matrix(item_discrimination, ncol = 1)
  
  fit_ability(
    response_data = response_data,
    item_loadings = item_loadings,
    threshold_left = threshold_left,
    threshold_right = threshold_right,
    ...
  )
}
```

#### Testing During Development
```r
# tests/testthat/test-main-interface.R
test_that("fit_ability works with MAP method", {
  data <- create_test_response_data()
  loadings <- matrix(rep(1, nrow(data)), ncol = 1)
  
  fit <- fit_ability(
    response_data = data,
    item_loadings = loadings,
    threshold_left = rep(-Inf, nrow(data)),
    threshold_right = rep(Inf, nrow(data)),
    method = "map"
  )
  
  expect_s3_class(fit, "irt_fit")
})

test_that("fit_unidim_ability works", {
  data <- create_test_response_data()
  
  fit <- fit_unidim_ability(
    response_data = data,
    item_discrimination = rep(1, nrow(data)),
    threshold_left = rep(-Inf, nrow(data)),
    threshold_right = rep(Inf, nrow(data)),
    method = "map"
  )
  
  expect_s3_class(fit, "irt_fit")
})

test_that("fscores works with 0-indexed responses", {
  # Create wide format data with 0-indexed responses
  wide_data <- data.frame(
    person_id = 1:10,
    age = rnorm(10, 30, 10),
    item1 = sample(0:2, 10, replace = TRUE),  # 0-indexed!
    item2 = sample(0:1, 10, replace = TRUE)   # 0-indexed!
  )
  
  item_params <- list(
    item1 = list(loadings = c(1.2), thresholds = c(-1, 0.5)),
    item2 = list(loadings = c(1.5), thresholds = c(0.0))
  )
  
  fit <- fscores(
    data = wide_data,
    item_params = item_params,
    person_id_col = "person_id",
    covariate_cols = "age",
    method = "map"
  )
  
  expect_s3_class(fit, "irt_fit")
})

test_that("fscores rejects 1-indexed responses", {
  wide_data <- data.frame(
    person_id = 1:10,
    item1 = sample(1:3, 10, replace = TRUE)  # Wrong! Should be 0:2
  )
  
  item_params <- list(
    item1 = list(loadings = c(1.2), thresholds = c(-1, 0.5))
  )
  
  expect_error(
    fscores(
      data = wide_data,
      item_params = item_params,
      person_id_col = "person_id",
      method = "map"
    ),
    "must be >= 0"
  )
})

test_that("fscores validates sufficient thresholds", {
  # Response value of 2 but only 1 threshold provided
  wide_data <- data.frame(
    person_id = 1:5,
    item1 = c(0, 1, 2, 1, 0)  # Max response is 2
  )
  
  item_params <- list(
    item1 = list(
      loadings = c(1.2), 
      thresholds = c(-1.0)  # Only 1 threshold, but need 2 for response value 2
    )
  )
  
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

test_that("fscores accepts correct number of thresholds", {
  # Response value of 2 with 2 thresholds provided
  wide_data <- data.frame(
    person_id = 1:5,
    item1 = c(0, 1, 2, 1, 0)  # Max response is 2
  )
  
  item_params <- list(
    item1 = list(
      loadings = c(1.2), 
      thresholds = c(-1.0, 0.5)  # 2 thresholds for 3 categories (0, 1, 2)
    )
  )
  
  expect_no_error(
    fscores(
      data = wide_data,
      item_params = item_params,
      person_id_col = "person_id",
      method = "map"
    )
  )
})
```

#### Phase Completion

**Git Commit:**
```bash
git add .
git commit -m "Phase 5 Complete: Main user interface

- Implemented unified fit_irt() function
- Added method routing to all estimation functions
- Created fit_unidim_irt() convenience wrapper
- Added comprehensive argument validation
- Implemented clear error messages
- All main interface tests passing
- User-facing API complete and tested
"
```

**Load Next Phase:**
```
📋 PHASE 5 COMPLETE ✓

Next Steps: Begin Phase 6 - S3 Methods and Utilities

Phase 6 will focus on:
- Implementing print.irt_fit() method
- Implementing summary.irt_fit() method
- Implementing coef.irt_fit() method
- Implementing predict.irt_fit() method
- Creating utility functions for extraction
- Testing all S3 methods
- Ensuring consistent output formatting

Ready to proceed with Phase 6? [Y/n]
```

---

### Phase 6: S3 Methods and Utilities

#### Tasks

1. **Implement `print.irt_fit()`**

```r
#' @export
print.irt_fit <- function(x, ...) {
  cat("IRT Model Fit\n")
  cat("=============\n\n")
  cat("Method:", x$method, "\n")
  cat("Dimensions:", x$stan_data$K, "\n")
  cat("Number of persons:", x$stan_data$N, "\n")
  cat("Number of responses:", x$stan_data$L, "\n")
  
  if(x$method == "MAP") {
    cat("Optimization converged:", x$return_code == 0, "\n")
    cat("Log posterior:", round(x$value, 2), "\n")
  } else if(x$method == "Laplace") {
    cat("Mode found successfully\n")
    cat("Log posterior at mode:", round(x$raw_fit$value, 2), "\n")
    cat("Posterior samples drawn:", nrow(x$samples$beta0), "\n")
  } else if(x$method == "Laplace-Failed") {
    cat("Warning: Laplace approximation failed\n")
    cat("Returning MAP estimate only\n")
    cat("Error:", x$error, "\n")
  } else if(x$method == "MCMC") {
    cat("Chains:", attr(x$raw_fit, "stan_args")[[1]]$chain_id, "\n")
    cat("Iterations:", attr(x$raw_fit, "stan_args")[[1]]$iter, "\n")
  }
  
  invisible(x)
}
```

2. **Implement `summary.irt_fit()`**

```r
#' @export
summary.irt_fit <- function(object, ...) {
  cat("IRT Model Summary\n")
  cat("=================\n\n")
  
  # Print basic info
  print(object)
  
  cat("\nPerson Ability Estimates:\n")
  abilities <- get_ability_estimates(object)
  print(summary(abilities))
  
  # Show score decomposition if covariates present
  if(object$stan_data$J > 1 || any(object$stan_data$X != 0)) {
    cat("\nScore Decomposition:\n")
    cat("Total Ability = Covariate Effects + Individual Effects\n")
    
    individual_effects <- get_individual_effects(object)
    covariate_effects <- get_covariate_effects(object)
    
    if(object$method %in% c("MAP", "optim", "nloptr")) {
      cat("\nCovariate Effects Summary:\n")
      print(summary(covariate_effects))
      cat("\nIndividual Effects Summary:\n")
      print(summary(individual_effects))
    }
  }
  
  if(object$stan_data$K > 1) {
    cat("\nDimension Correlations:\n")
    print(get_dimension_correlations(object))
  }
  
  invisible(object)
}
```

3. **Implement `coef.irt_fit()`**

```r
#' @export
coef.irt_fit <- function(object, 
                         pars = c("beta0", "beta", "tau", "Omega"),
                         ...) {
  
  if(object$method %in% c("MAP", "optim", "nloptr")) {
    # Point estimates
    extract_point_estimates(object, pars)
  } else if(object$method == "Laplace-Failed") {
    # Only MAP available
    extract_point_estimates(object, pars)
  } else {
    # Posterior summaries (VB, MCMC, Laplace)
    extract_posterior_summaries(object, pars)
  }
}
```

4. **Implement `predict.irt_fit()`** (get person abilities)

```r
#' @export
predict.irt_fit <- function(object, 
                           type = c("abilities", "individual_effects", 
                                    "covariate_effects", "response_probs"),
                           summary = TRUE,
                           ...) {
  
  type <- match.arg(type)
  
  if(type == "abilities") {
    # Return full theta: beta0 + X*b + eta
    # This is the complete ability score including both covariates and individual effects
    abilities <- get_ability_estimates(object)
    
    if(summary && object$method %in% c("vb", "mcmc", "laplace")) {
      # Return posterior summaries
      summarize_abilities(abilities)
    } else {
      abilities
    }
    
  } else if(type == "individual_effects") {
    # Return just eta (person-specific random effects/deviations)
    individual_effects <- get_individual_effects(object)
    
    if(summary && object$method %in% c("vb", "mcmc", "laplace")) {
      summarize_abilities(individual_effects)
    } else {
      individual_effects
    }
    
  } else if(type == "covariate_effects") {
    # Return just beta0 + X*b (contribution from covariates)
    get_covariate_effects(object)
    
  } else if(type == "response_probs") {
    # Calculate expected response probabilities
    calculate_response_probabilities(object)
  }
}
```

5. **Create utility functions**

```r
#' Extract ability estimates (full theta = beta0 + X*b + eta)
#' @keywords internal
get_ability_estimates <- function(fit) {
  if(fit$method %in% c("MAP", "optim", "nloptr", "Laplace-Failed")) {
    # For point estimates, need to reconstruct theta from components
    
    # Extract beta0
    beta0_names <- grep("^beta0\\[", names(fit$estimates), value = TRUE)
    beta0 <- fit$estimates[beta0_names]
    
    # Extract b (covariate coefficients)
    b_list <- list()
    for(k in 1:fit$stan_data$K) {
      b_names <- grep(paste0("^b\\[", k, ","), names(fit$estimates), value = TRUE)
      b_list[[k]] <- fit$estimates[b_names]
    }
    
    # Extract eta (individual effects)
    eta_names <- grep("^eta\\[", names(fit$estimates), value = TRUE)
    eta_values <- fit$estimates[eta_names]
    eta <- matrix(eta_values, ncol = fit$stan_data$K)
    
    # Reconstruct theta = beta0 + X*b + eta
    theta <- matrix(0, nrow = fit$stan_data$N, ncol = fit$stan_data$K)
    for(k in 1:fit$stan_data$K) {
      theta[, k] <- beta0[k] + fit$stan_data$X %*% b_list[[k]] + eta[, k]
    }
    
    return(theta)
    
  } else {
    # For MCMC/VB/Laplace, extract theta directly if available, otherwise reconstruct
    if("theta" %in% names(fit$samples)) {
      return(fit$samples$theta)
    } else {
      # Reconstruct from components
      beta0 <- fit$samples$beta0
      b <- fit$samples$b
      eta <- fit$samples$eta
      
      # theta will be S x N x K (samples x persons x dimensions)
      S <- dim(eta)[1]
      N <- fit$stan_data$N
      K <- fit$stan_data$K
      
      theta <- array(0, dim = c(S, N, K))
      for(k in 1:K) {
        for(s in 1:S) {
          theta[s, , k] <- beta0[s, k] + fit$stan_data$X %*% b[s, k, ] + eta[s, , k]
        }
      }
      
      return(theta)
    }
  }
}

#' Extract individual effects only (eta)
#' @keywords internal
get_individual_effects <- function(fit) {
  if(fit$method %in% c("MAP", "optim", "nloptr", "Laplace-Failed")) {
    # Extract eta parameters
    eta_names <- grep("^eta\\[", names(fit$estimates), value = TRUE)
    eta_values <- fit$estimates[eta_names]
    
    # Reshape to N x K matrix
    matrix(eta_values, ncol = fit$stan_data$K)
    
  } else {
    # Extract from posterior samples
    fit$samples$eta
  }
}

#' Extract covariate effects (beta0 + X*b)
#' @keywords internal
get_covariate_effects <- function(fit) {
  if(fit$method %in% c("MAP", "optim", "nloptr")) {
    # Extract beta0
    beta0_names <- grep("^beta0\\[", names(fit$estimates), value = TRUE)
    beta0 <- fit$estimates[beta0_names]
    
    # Extract b (covariate coefficients)
    b_list <- list()
    for(k in 1:fit$stan_data$K) {
      b_names <- grep(paste0("^b\\[", k, ","), names(fit$estimates), value = TRUE)
      b_list[[k]] <- fit$estimates[b_names]
    }
    
    # Calculate beta0 + X*b for each dimension
    covariate_effects <- matrix(0, nrow = fit$stan_data$N, ncol = fit$stan_data$K)
    for(k in 1:fit$stan_data$K) {
      covariate_effects[, k] <- beta0[k] + fit$stan_data$X %*% b_list[[k]]
    }
    
    return(covariate_effects)
    
  } else {
    # For MCMC/VB, compute for each sample
    beta0 <- fit$samples$beta0
    b <- fit$samples$b
    
    S <- dim(beta0)[1]
    N <- fit$stan_data$N
    K <- fit$stan_data$K
    
    covariate_effects <- array(0, dim = c(S, N, K))
    for(k in 1:K) {
      for(s in 1:S) {
        covariate_effects[s, , k] <- beta0[s, k] + fit$stan_data$X %*% b[s, k, ]
      }
    }
    
    return(covariate_effects)
  }
}

#' Extract dimension correlations
#' @keywords internal
get_dimension_correlations <- function(fit) {
  if(fit$stan_data$K == 1) {
    return(matrix(1, 1, 1))
  }
  
  if(fit$method %in% c("MAP", "optim", "nloptr")) {
    # Extract Omega from point estimate
    omega_names <- grep("^Omega\\[", names(fit$estimates), value = TRUE)
    omega_values <- fit$estimates[omega_names]
    matrix(omega_values, nrow = fit$stan_data$K)
    
  } else {
    # Return posterior mean
    apply(fit$samples$Omega, c(2, 3), mean)
  }
}
```

#### Testing During Development
```r
# tests/testthat/test-methods.R
test_that("print method works", {
  fit <- create_simple_fit_object()
  expect_output(print(fit), "IRT Model Fit")
})

test_that("summary method works", {
  fit <- create_simple_fit_object()
  expect_output(summary(fit), "Person Ability")
})

test_that("coef method extracts parameters", {
  fit <- create_simple_fit_object()
  coeffs <- coef(fit)
  expect_type(coeffs, "list")
})

test_that("predict returns abilities", {
  fit <- create_simple_fit_object()
  abilities <- predict(fit, type = "abilities")
  expect_true(is.matrix(abilities) || is.array(abilities))
})

test_that("predict returns individual effects", {
  fit <- create_simple_fit_object()
  ind_effects <- predict(fit, type = "individual_effects")
  expect_true(is.matrix(ind_effects) || is.array(ind_effects))
})

test_that("predict returns covariate effects", {
  fit <- create_simple_fit_object()
  cov_effects <- predict(fit, type = "covariate_effects")
  expect_true(is.matrix(cov_effects) || is.array(cov_effects))
})

test_that("abilities equal sum of components", {
  fit <- create_simple_fit_object()
  
  abilities <- predict(fit, type = "abilities")
  ind_effects <- predict(fit, type = "individual_effects")
  cov_effects <- predict(fit, type = "covariate_effects")
  
  # For point estimates: theta = covariate_effects + individual_effects
  expect_equal(abilities, cov_effects + ind_effects, tolerance = 1e-6)
})
```

#### Phase Completion

**Git Commit:**
```bash
git add .
git commit -m "Phase 6 Complete: S3 methods and utilities

- Implemented print.irt_fit() for clean output
- Implemented summary.irt_fit() with detailed info
- Implemented coef.irt_fit() for parameter extraction
- Implemented predict.irt_fit() for ability scores
- Created get_ability_estimates() utility
- Created get_dimension_correlations() utility
- Added response probability calculations
- All S3 method tests passing
"
```

**Load Next Phase:**
```
📋 PHASE 6 COMPLETE ✓

Next Steps: Begin Phase 7 - Documentation and Vignettes

Phase 7 will focus on:
- Writing roxygen2 documentation for all functions
- Creating comprehensive README.md
- Writing quickstart vignette
- Writing unidimensional IRT vignette
- Writing multidimensional IRT vignette
- Writing estimation methods comparison vignette
- Generating .Rd files with devtools::document()

Ready to proceed with Phase 7? [Y/n]
```

---

### Phase 7: Documentation and Vignettes

#### Tasks

1. **Write comprehensive function documentation**
   - Add roxygen2 comments to all exported functions
   - Include parameter descriptions, return values, examples
   - Run `devtools::document()` to generate .Rd files

2. **Create README.md with quick examples**

```markdown
# IRTscoring: Flexible IRT Scoring in R

## Installation

```r
# Install from GitHub
devtools::install_github("yourusername/IRTscoring")
```

## Quick Start

### Unidimensional IRT

```r
library(IRTscoring)

# Fit model with MAP estimation (responses are 0-indexed!)
fit <- fit_unidim_ability(
  response_data = my_responses,
  item_discrimination = rep(1, nrow(my_responses)),
  threshold_left = c(-Inf, -1, 0),
  threshold_right = c(-1, 0, Inf),
  method = "map"
)

# Get ability estimates
abilities <- predict(fit)
```

### Multidimensional IRT

```r
# Define 2D loading structure
loadings <- matrix(c(
  1, 0,
  1, 0,
  0, 1,
  0, 1
), ncol = 2, byrow = TRUE)

fit <- fit_ability(
  response_data = my_responses,
  item_loadings = loadings,
  threshold_left = thresholds_left,
  threshold_right = thresholds_right,
  method = "mcmc",
  chains = 4
)
```

### Wide Format (User-Friendly)

```r
# Data in wide format (one row per person)
# IMPORTANT: Responses must be 0-indexed (0, 1, 2, ...)
wide_data <- data.frame(
  person_id = 1:100,
  age = rnorm(100, 30, 10),
  item1 = sample(0:2, 100, replace = TRUE),
  item2 = sample(0:2, 100, replace = TRUE)
)

# Define item parameters
item_params <- list(
  item1 = list(loadings = c(1.2, 0.1), thresholds = c(-1, 0.5)),
  item2 = list(loadings = c(0.8, 1.5), thresholds = c(-0.5, 1.0))
)

# Fit model
fit <- fscores(
  data = wide_data,
  item_params = item_params,
  person_id_col = "person_id",
  covariate_cols = "age",
  method = "map"
)

# Get ability scores
abilities <- predict(fit, type = "abilities")
```
```

3. **Create vignettes**

```r
usethis::use_vignette("quickstart")
usethis::use_vignette("unidimensional-irt")
usethis::use_vignette("multidimensional-irt")
usethis::use_vignette("estimation-methods")
```

#### Example Vignette Structure

```markdown
---
title: "Multidimensional IRT with irtscoring"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Multidimensional IRT with irtscoring}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction

This vignette demonstrates fitting multidimensional IRT models...

## Data Setup

## Model Specification

## Estimation

## Interpreting Results

## Visualization
```

#### Phase Completion

**Git Commit:**
```bash
git add .
git commit -m "Phase 7 Complete: Documentation and vignettes

- Added roxygen2 documentation to all exported functions
- Created comprehensive README.md with examples
- Wrote quickstart vignette
- Wrote unidimensional IRT vignette
- Wrote multidimensional IRT vignette
- Wrote estimation methods vignette
- Generated all .Rd manual files
- Documentation complete and ready for users
"
```

**Load Next Phase:**
```
📋 PHASE 7 COMPLETE ✓

Next Steps: Final Testing and Validation

Final phase will focus on:
- Running comprehensive unit test suite
- Executing integration tests
- Running performance benchmarks
- Conducting simulation study vs. mirt package
- Running R CMD check
- Generating code coverage report
- Validating against known results
- Creating final validation report

Ready to proceed with Final Testing? [Y/n]
```

---

## Final Testing and Validation

### Comprehensive Test Suite

#### Unit Tests (Already created during phases)
- Data preparation functions
- Stan estimation methods
- BridgeStan integration
- S3 methods
- Utility functions

#### Integration Tests

```r
# tests/testthat/test-integration.R

test_that("Full workflow: unidimensional MAP", {
  # Create synthetic data
  set.seed(123)
  data <- simulate_irt_data(N = 100, I = 10, K = 1)
  
  # Fit model
  fit <- fit_unidim_irt(
    response_data = data$responses,
    item_discrimination = data$true_loadings,
    threshold_left = data$thresh_left,
    threshold_right = data$thresh_right,
    method = "map"
  )
  
  # Check recovery
  abilities_true <- data$true_abilities
  abilities_est <- predict(fit)
  
  correlation <- cor(abilities_true, abilities_est)
  expect_gt(correlation, 0.8)
})

test_that("Full workflow: multidimensional MCMC", {
  skip_on_cran()
  
  set.seed(123)
  data <- simulate_irt_data(N = 50, I = 20, K = 2)
  
  fit <- fit_irt(
    response_data = data$responses,
    item_loadings = data$true_loadings,
    threshold_left = data$thresh_left,
    threshold_right = data$thresh_right,
    method = "mcmc",
    chains = 2,
    iter = 500,
    warmup = 250
  )
  
  expect_s3_class(fit, "irt_fit")
  
  # Check correlation matrix recovery
  omega_true <- data$true_correlation
  omega_est <- get_dimension_correlations(fit)
  
  expect_equal(dim(omega_est), c(2, 2))
})
```

#### Performance Tests

```r
# tests/testthat/test-performance.R

test_that("Large dataset performance is acceptable", {
  skip_on_cran()
  
  # Create large dataset
  data <- simulate_irt_data(N = 1000, I = 50, K = 2)
  
  # Time the estimation
  time_taken <- system.time({
    fit <- fit_irt(
      response_data = data$responses,
      item_loadings = data$true_loadings,
      threshold_left = data$thresh_left,
      threshold_right = data$thresh_right,
      method = "map"
    )
  })
  
  # Should complete in reasonable time (adjust threshold as needed)
  expect_lt(time_taken["elapsed"], 60)  # 60 seconds
})
```

#### Comparison Tests

```r
# tests/testthat/test-comparison.R

test_that("Different estimation methods give similar results", {
  skip_on_cran()
  
  set.seed(123)
  data <- simulate_irt_data(N = 100, I = 10, K = 1)
  
  # MAP
  fit_map <- fit_irt(
    response_data = data$responses,
    item_loadings = matrix(data$true_loadings, ncol = 1),
    threshold_left = data$thresh_left,
    threshold_right = data$thresh_right,
    method = "map"
  )
  
  # VB
  fit_vb <- fit_irt(
    response_data = data$responses,
    item_loadings = matrix(data$true_loadings, ncol = 1),
    threshold_left = data$thresh_left,
    threshold_right = data$thresh_right,
    method = "vb",
    use_pathfinder = FALSE,
    iter = 10000
  )
  
  abilities_map <- predict(fit_map)
  abilities_vb <- predict(fit_vb, summary = TRUE)[, "mean"]
  
  # Should be highly correlated
  expect_gt(cor(abilities_map, abilities_vb), 0.95)
})
```

### Package Validation

#### R CMD check
```r
# Run comprehensive checks
devtools::check()

# Should pass with:
# - 0 errors
# - 0 warnings
# - 0 notes (or only acceptable notes)
```

#### Code Coverage
```r
# Install covr
install.packages("covr")

# Generate coverage report
covr::package_coverage()

# Aim for >80% code coverage
```

#### Reverse Dependency Checks
```r
# Once package is on CRAN, check reverse dependencies
devtools::revdep_check()
```

### Validation Against Known Results

```r
# tests/testthat/test-validation.R

test_that("Results match published examples", {
  # Load published dataset (if available)
  data("example_data", package = "IRTscoring")
  
  # Fit using known parameters
  fit <- fit_irt(
    response_data = example_data$responses,
    item_loadings = example_data$loadings,
    threshold_left = example_data$thresh_left,
    threshold_right = example_data$thresh_right,
    method = "map"
  )
  
  # Compare to published results
  published_abilities <- example_data$published_abilities
  estimated_abilities <- predict(fit)
  
  # Allow for small numerical differences
  expect_equal(estimated_abilities, published_abilities, tolerance = 0.01)
})
```

### Simulation Study: Comparison with mirt Package

A critical validation is comparing our MAP estimates with the established `mirt` package's `fscores()` function.

#### Create Simulation Data Using mirt

```r
# tests/testthat/test-mirt-comparison.R

library(mirt)
library(IRTscoring)

#' Generate simulated IRT data using mirt
#' @param N Number of persons
#' @param I Number of items
#' @param K Number of dimensions
#' @param categories Number of response categories per item
#' @return List with response data and true parameters
simulate_mirt_data <- function(N = 500, I = 20, K = 1, categories = 3) {
  
  # Define item parameters for graded response model
  # a = discrimination, d = difficulty/threshold parameters
  set.seed(123)
  
  if(K == 1) {
    # Unidimensional: a parameters (discrimination)
    a_params <- matrix(rlnorm(I, meanlog = 0, sdlog = 0.3), ncol = 1)
    
    # d parameters (thresholds) - need (categories - 1) per item
    d_params <- t(sapply(1:I, function(i) {
      sort(rnorm(categories - 1, mean = 0, sd = 1))
    }))
    
    # Combine into parameter matrix for mirt
    item_params <- cbind(a_params, d_params)
    
    # Specify model
    model <- mirt::mirt.model('F1 = 1-' %s% I)
    
  } else if(K == 2) {
    # Multidimensional: half items load on F1, half on F2
    a_params <- matrix(0, nrow = I, ncol = 2)
    a_params[1:(I/2), 1] <- rlnorm(I/2, meanlog = 0, sdlog = 0.3)
    a_params[(I/2 + 1):I, 2] <- rlnorm(I/2, meanlog = 0, sdlog = 0.3)
    
    # d parameters
    d_params <- t(sapply(1:I, function(i) {
      sort(rnorm(categories - 1, mean = 0, sd = 1))
    }))
    
    item_params <- cbind(a_params, d_params)
    
    # Specify model
    model <- mirt::mirt.model('
      F1 = 1-' %s% (I/2) %s% '
      F2 = ' %s% (I/2 + 1) %s% '-' %s% I %s% '
      COV = F1*F2
    ')
    
  } else {
    stop("K > 2 not implemented in this helper")
  }
  
  # Generate ability scores
  if(K == 1) {
    theta <- matrix(rnorm(N), ncol = 1)
  } else {
    # Correlated abilities
    Sigma <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
    theta <- MASS::mvrnorm(N, mu = rep(0, K), Sigma = Sigma)
  }
  
  # Simulate responses
  responses <- mirt::simdata(
    a = a_params,
    d = d_params,
    N = N,
    itemtype = 'graded',
    Theta = theta
  )
  
  list(
    responses = responses,
    theta_true = theta,
    a_params = a_params,
    d_params = d_params,
    item_params = item_params,
    K = K,
    categories = categories
  )
}

#' Convert mirt response data to irtscoring format
#' @param mirt_responses Matrix from mirt simulation (1-indexed: categories 1, 2, 3, ...)
#' @param a_params Item discrimination parameters
#' @param d_params Threshold parameters
#' @return List with response_data, loadings, thresholds
#' @note This function converts from mirt's 1-indexed responses to IRTscoring's 0-indexed responses
convert_mirt_to_irtscoring <- function(mirt_responses, a_params, d_params) {
  
  N <- nrow(mirt_responses)
  I <- ncol(mirt_responses)
  K <- ncol(a_params)
  
  # Create long-format response data
  response_list <- list()
  loading_list <- list()
  thresh_left_list <- list()
  thresh_right_list <- list()
  
  for(person in 1:N) {
    for(item in 1:I) {
      response_cat <- mirt_responses[person, item]
      
      # IMPORTANT: mirt uses 1-indexed responses (1, 2, 3, ...)
      # IRTscoring uses 0-indexed responses (0, 1, 2, ...)
      # Convert: mirt category k becomes IRTscoring category k-1
      response_cat_0indexed <- response_cat - 1
      
      # For each response category, create a row
      # Category coding: if response = k (0-indexed), then we observed category k
      # We need to set up thresholds for the probability of that category
      
      # In graded response model:
      # P(Y >= k) = invlogit(a * theta - d[k])
      # P(Y = k) = P(Y >= k) - P(Y >= k+1)
      
      # For our model: P = invlogit(theta + dL) - invlogit(theta + dR)
      # So we need to convert mirt's d parameters to our threshold format
      
      # mirt: P(Y >= k) = invlogit(a * theta - d[k])
      # ours: P = invlogit(a'*theta + threshold)
      # So threshold = -d in mirt parameterization (when a is in the loading)
      
      response_list <- c(response_list, list(
        data.frame(
          pid = person,
          iid = item,
          response = response_cat_0indexed  # 0-indexed!
        )
      ))
      
      # Loadings (L x K matrix where L is number of response observations)
      loading_list <- c(loading_list, list(a_params[item, , drop = FALSE]))
      
      # Thresholds for this response category
      n_cats <- length(d_params[item, ]) + 1
      
      # Use 1-indexed response_cat for threshold lookup (mirt convention)
      if(response_cat == 1) {
        # Lowest category: (-Inf, d1]
        thresh_left_list <- c(thresh_left_list, -Inf)
        thresh_right_list <- c(thresh_right_list, -d_params[item, 1])
      } else if(response_cat == n_cats) {
        # Highest category: (d[K-1], Inf)
        thresh_left_list <- c(thresh_left_list, -d_params[item, n_cats - 1])
        thresh_right_list <- c(thresh_right_list, Inf)
      } else {
        # Middle categories: (d[k-1], d[k]]
        thresh_left_list <- c(thresh_left_list, -d_params[item, response_cat - 1])
        thresh_right_list <- c(thresh_right_list, -d_params[item, response_cat])
      }
    }
  }
  
  response_data <- do.call(rbind, response_list)
  loadings <- do.call(rbind, loading_list)
  thresh_left <- unlist(thresh_left_list)
  thresh_right <- unlist(thresh_right_list)
  
  list(
    response_data = response_data,
    item_loadings = loadings,
    threshold_left = thresh_left,
    threshold_right = thresh_right
  )
}

test_that("Unidimensional MAP estimates match mirt::fscores", {
  skip_if_not_installed("mirt")
  
  # Generate data using mirt
  sim_data <- simulate_mirt_data(N = 200, I = 10, K = 1, categories = 3)
  
  # Fit with mirt
  mirt_fit <- mirt::mirt(
    data = sim_data$responses,
    model = 1,  # unidimensional
    itemtype = 'graded',
    verbose = FALSE
  )
  
  # Get MAP scores from mirt
  mirt_scores <- mirt::fscores(mirt_fit, method = "MAP", full.scores = FALSE)
  
  # Convert to irtscoring format
  irtscoring_data <- convert_mirt_to_irtscoring(
    sim_data$responses,
    sim_data$a_params,
    sim_data$d_params
  )
  
  # Fit with irtscoring
  irt_fit <- fit_ability(
    response_data = irtscoring_data$response_data,
    item_loadings = irtscoring_data$item_loadings,
    threshold_left = irtscoring_data$threshold_left,
    threshold_right = irtscoring_data$threshold_right,
    method = "map"
  )
  
  # Get scores from irtscoring
  irt_scores <- predict(irt_fit, type = "abilities")
  
  # Compare
  correlation <- cor(mirt_scores[, 1], irt_scores[, 1])
  
  # Scores should be highly correlated (> 0.95)
  expect_gt(correlation, 0.95)
  
  # RMSE should be small
  rmse <- sqrt(mean((mirt_scores[, 1] - irt_scores[, 1])^2))
  expect_lt(rmse, 0.2)
  
  # Print diagnostic info
  cat("\nUnidimensional Comparison:\n")
  cat("Correlation:", round(correlation, 4), "\n")
  cat("RMSE:", round(rmse, 4), "\n")
})

test_that("Multidimensional MAP estimates match mirt::fscores", {
  skip_if_not_installed("mirt")
  skip_on_cran()  # Takes longer
  
  # Generate 2D data using mirt
  sim_data <- simulate_mirt_data(N = 300, I = 20, K = 2, categories = 3)
  
  # Fit with mirt
  mirt_model <- mirt::mirt.model('
    F1 = 1-10
    F2 = 11-20
    COV = F1*F2
  ')
  
  mirt_fit <- mirt::mirt(
    data = sim_data$responses,
    model = mirt_model,
    itemtype = 'graded',
    verbose = FALSE
  )
  
  # Get MAP scores from mirt
  mirt_scores <- mirt::fscores(mirt_fit, method = "MAP", full.scores = FALSE)
  
  # Convert to irtscoring format
  irtscoring_data <- convert_mirt_to_irtscoring(
    sim_data$responses,
    sim_data$a_params,
    sim_data$d_params
  )
  
  # Fit with irtscoring
  irt_fit <- fit_ability(
    response_data = irtscoring_data$response_data,
    item_loadings = irtscoring_data$item_loadings,
    threshold_left = irtscoring_data$threshold_left,
    threshold_right = irtscoring_data$threshold_right,
    method = "map"
  )
  
  # Get scores from irtscoring
  irt_scores <- predict(irt_fit, type = "abilities")
  
  # Compare dimension 1
  cor_dim1 <- cor(mirt_scores[, 1], irt_scores[, 1])
  rmse_dim1 <- sqrt(mean((mirt_scores[, 1] - irt_scores[, 1])^2))
  
  # Compare dimension 2
  cor_dim2 <- cor(mirt_scores[, 2], irt_scores[, 2])
  rmse_dim2 <- sqrt(mean((mirt_scores[, 2] - irt_scores[, 2])^2))
  
  # Both dimensions should match well
  expect_gt(cor_dim1, 0.95)
  expect_gt(cor_dim2, 0.95)
  expect_lt(rmse_dim1, 0.2)
  expect_lt(rmse_dim2, 0.2)
  
  # Print diagnostic info
  cat("\nMultidimensional Comparison:\n")
  cat("Dimension 1 - Correlation:", round(cor_dim1, 4), "RMSE:", round(rmse_dim1, 4), "\n")
  cat("Dimension 2 - Correlation:", round(cor_dim2, 4), "RMSE:", round(rmse_dim2, 4), "\n")
  
  # Compare correlation between dimensions
  mirt_cor <- cor(mirt_scores[, 1], mirt_scores[, 2])
  irt_cor <- cor(irt_scores[, 1], irt_scores[, 2])
  
  cat("Between-dimension correlation - mirt:", round(mirt_cor, 3), 
      "irtscoring:", round(irt_cor, 3), "\n")
})

test_that("Parameter recovery matches across methods", {
  skip_if_not_installed("mirt")
  skip_on_cran()
  
  # Generate data with known parameters
  sim_data <- simulate_mirt_data(N = 500, I = 15, K = 1, categories = 4)
  
  # Fit with mirt
  mirt_fit <- mirt::mirt(
    data = sim_data$responses,
    model = 1,
    itemtype = 'graded',
    verbose = FALSE
  )
  
  # Extract item parameters from mirt
  mirt_coefs <- coef(mirt_fit, simplify = TRUE)
  mirt_items <- mirt_coefs$items
  
  # Convert to irtscoring format and fit
  irtscoring_data <- convert_mirt_to_irtscoring(
    sim_data$responses,
    sim_data$a_params,
    sim_data$d_params
  )
  
  irt_fit <- fit_ability(
    response_data = irtscoring_data$response_data,
    item_loadings = irtscoring_data$item_loadings,
    threshold_left = irtscoring_data$threshold_left,
    threshold_right = irtscoring_data$threshold_right,
    method = "map"
  )
  
  # Compare to true parameters used in simulation
  true_a <- sim_data$a_params[, 1]
  true_d <- sim_data$d_params
  
  # mirt estimates
  mirt_a <- mirt_items[, 1]
  mirt_d <- mirt_items[, 2:ncol(mirt_items)]
  
  # Check recovery
  cor_a_mirt <- cor(true_a, mirt_a)
  cor_d_mirt <- cor(as.vector(true_d), as.vector(mirt_d))
  
  expect_gt(cor_a_mirt, 0.90)
  expect_gt(cor_d_mirt, 0.90)
  
  cat("\nParameter Recovery:\n")
  cat("Discrimination (a) correlation with true:", round(cor_a_mirt, 4), "\n")
  cat("Threshold (d) correlation with true:", round(cor_d_mirt, 4), "\n")
})
```

#### Comprehensive Simulation Study (for vignette/paper)

```r
# vignettes/simulation-study.Rmd or scripts/run_simulation_study.R

#' Run comprehensive comparison simulation
#' 
#' Tests multiple conditions:
#' - Sample sizes: N = 100, 200, 500, 1000
#' - Number of items: I = 10, 20, 30
#' - Dimensions: K = 1, 2, 3
#' - Categories: 2, 3, 4, 5
#' 
#' @param n_replications Number of replications per condition
#' @return Data frame with results
run_simulation_study <- function(n_replications = 100) {
  
  conditions <- expand.grid(
    N = c(100, 200, 500, 1000),
    I = c(10, 20, 30),
    K = c(1, 2, 3),
    categories = c(2, 3, 4, 5),
    replication = 1:n_replications
  )
  
  results <- data.frame()
  
  for(i in 1:nrow(conditions)) {
    cond <- conditions[i, ]
    
    tryCatch({
      # Generate data
      if(cond$K <= 2) {
        sim_data <- simulate_mirt_data(
          N = cond$N,
          I = cond$I,
          K = cond$K,
          categories = cond$categories
        )
        
        # Fit with mirt
        mirt_start <- Sys.time()
        mirt_fit <- mirt::mirt(
          data = sim_data$responses,
          model = cond$K,
          itemtype = 'graded',
          verbose = FALSE
        )
        mirt_time <- as.numeric(Sys.time() - mirt_start)
        mirt_scores <- mirt::fscores(mirt_fit, method = "MAP")
        
        # Fit with irtscoring
        irtscoring_data <- convert_mirt_to_irtscoring(
          sim_data$responses,
          sim_data$a_params,
          sim_data$d_params
        )
        
        irt_start <- Sys.time()
        irt_fit <- fit_ability(
          response_data = irtscoring_data$response_data,
          item_loadings = irtscoring_data$item_loadings,
          threshold_left = irtscoring_data$threshold_left,
          threshold_right = irtscoring_data$threshold_right,
          method = "map"
        )
        irt_time <- as.numeric(Sys.time() - irt_start)
        irt_scores <- predict(irt_fit, type = "abilities")
        
        # Calculate agreement metrics
        correlations <- numeric(cond$K)
        rmses <- numeric(cond$K)
        
        for(k in 1:cond$K) {
          correlations[k] <- cor(mirt_scores[, k], irt_scores[, k])
          rmses[k] <- sqrt(mean((mirt_scores[, k] - irt_scores[, k])^2))
        }
        
        # Store results
        result_row <- data.frame(
          N = cond$N,
          I = cond$I,
          K = cond$K,
          categories = cond$categories,
          replication = cond$replication,
          mean_correlation = mean(correlations),
          min_correlation = min(correlations),
          mean_rmse = mean(rmses),
          max_rmse = max(rmses),
          mirt_time = mirt_time,
          irtscoring_time = irt_time,
          convergence = TRUE
        )
        
        results <- rbind(results, result_row)
      }
      
    }, error = function(e) {
      # Log failures
      result_row <- data.frame(
        N = cond$N,
        I = cond$I,
        K = cond$K,
        categories = cond$categories,
        replication = cond$replication,
        mean_correlation = NA,
        min_correlation = NA,
        mean_rmse = NA,
        max_rmse = NA,
        mirt_time = NA,
        irtscoring_time = NA,
        convergence = FALSE
      )
      results <- rbind(results, result_row)
    })
    
    # Progress
    if(i %% 10 == 0) {
      cat("Completed", i, "of", nrow(conditions), "conditions\n")
    }
  }
  
  return(results)
}

# Run the study
simulation_results <- run_simulation_study(n_replications = 100)

# Save results
saveRDS(simulation_results, "simulation_results.rds")

# Summarize results
library(dplyr)

summary_table <- simulation_results %>%
  group_by(N, I, K, categories) %>%
  summarise(
    mean_corr = mean(mean_correlation, na.rm = TRUE),
    sd_corr = sd(mean_correlation, na.rm = TRUE),
    mean_rmse = mean(mean_rmse, na.rm = TRUE),
    convergence_rate = mean(convergence),
    relative_time = mean(irtscoring_time / mirt_time, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_table)

# Create plots
library(ggplot2)

ggplot(simulation_results, aes(x = factor(N), y = mean_correlation)) +
  geom_boxplot() +
  facet_grid(K ~ categories, labeller = label_both) +
  labs(title = "Agreement between irtscoring and mirt",
       x = "Sample Size", y = "Correlation") +
  theme_bw()
```

---

## Example Usage Code

### Understanding Score Components

The package provides three types of scores:

1. **Full Ability Scores** (`type = "abilities"`): θ = β₀ + X·b + η
   - This is the **default** and most commonly used
   - Includes both covariate effects and individual-specific effects
   - This is what you typically want for scoring individuals

2. **Individual Effects** (`type = "individual_effects"`): η only
   - The person-specific random effects/deviations
   - Captures unique ability after controlling for covariates
   - Useful for understanding individual variation beyond covariates

3. **Covariate Effects** (`type = "covariate_effects"`): β₀ + X·b only
   - The predicted ability based solely on observed characteristics
   - Useful for understanding how covariates influence ability

### Basic Unidimensional Example

```r
library(IRTscoring)

# Create sample data with 0-indexed responses
set.seed(42)
N <- 200  # people
I <- 10   # items

# Person IDs and item IDs
responses <- data.frame(
  pid = rep(1:N, each = I),
  iid = rep(1:I, times = N)
)

# Item discriminations (all equal for simplicity)
discrimination <- rep(1.5, nrow(responses))

# Thresholds for 3 categories (0, 1, 2)
thresh_left <- rep(c(-Inf, -0.5), times = I * N / 2)
thresh_right <- rep(c(-0.5, Inf), times = I * N / 2)

# Fit model with MAP
fit <- fit_unidim_ability(
  response_data = responses,
  item_discrimination = discrimination,
  threshold_left = thresh_left,
  threshold_right = thresh_right,
  method = "map"
)

# View results
print(fit)
summary(fit)

# Extract full abilities (default - includes both covariates and individual effects)
abilities <- predict(fit, type = "abilities")
hist(abilities, main = "Estimated Ability Distribution")

# Extract just individual effects (person-specific deviations)
individual_effects <- predict(fit, type = "individual_effects")

# Extract just covariate effects (predicted ability from covariates)
covariate_effects <- predict(fit, type = "covariate_effects")

# Verify: abilities = covariate_effects + individual_effects
all.equal(abilities, covariate_effects + individual_effects)

# Extract parameters
params <- coef(fit)
print(params$beta0)
```

### Multidimensional Example with Covariates

```r
library(IRTscoring)

# Create 2-dimensional data with 0-indexed responses
set.seed(42)
N <- 200
I <- 20
L <- N * I

responses <- data.frame(
  pid = rep(1:N, each = I),
  iid = rep(1:I, times = N)
)

# Define loading structure (10 items per dimension)
loadings <- matrix(0, nrow = L, ncol = 2)
for(i in 1:L) {
  item <- responses$iid[i]
  if(item <= 10) {
    loadings[i, 1] <- 1.2
  } else {
    loadings[i, 2] <- 1.2
  }
}

# Create thresholds
thresh_left <- runif(L, -2, 0)
thresh_right <- runif(L, 0, 2)

# Add person covariates
person_covariates <- cbind(
  intercept = 1,
  age = rnorm(N),
  gender = rbinom(N, 1, 0.5)
)

# Fit with MCMC
fit <- fit_ability(
  response_data = responses,
  item_loadings = loadings,
  threshold_left = thresh_left,
  threshold_right = thresh_right,
  person_covariates = person_covariates,
  method = "mcmc",
  chains = 4,
  iter = 2000,
  warmup = 1000
)

# View results
print(fit)
summary(fit)

# Get dimension correlations
omega <- get_dimension_correlations(fit)
print(omega)

# Extract full abilities (includes both covariate and individual effects)
abilities <- predict(fit, type = "abilities", summary = TRUE)

# Extract individual effects only (person-specific deviations)
individual_effects <- predict(fit, type = "individual_effects", summary = TRUE)

# Extract covariate effects only (predicted from age and gender)
covariate_effects <- predict(fit, type = "covariate_effects")

# Plot full abilities in 2D space
plot(abilities[, 1], abilities[, 2], 
     xlab = "Dimension 1", ylab = "Dimension 2",
     main = "Person Abilities in 2D Space (Full Scores)")

# Plot individual effects in 2D space (after controlling for covariates)
plot(individual_effects[, 1], individual_effects[, 2],
     xlab = "Dimension 1", ylab = "Dimension 2",
     main = "Individual Effects (Residual Abilities)")

# Compare how much covariates explain
cor(abilities[, 1], covariate_effects[, 1])  # Correlation for dimension 1
```

### Using Alternative Optimizers

```r
library(IRTscoring)

# Same data setup as before...

# Fit using optim() via BridgeStan
fit_optim <- fit_ability(
  response_data = responses,
  item_loadings = loadings,
  threshold_left = thresh_left,
  threshold_right = thresh_right,
  method = "optim",
  control = list(maxit = 500)
)

# Fit using nloptr via BridgeStan
fit_nloptr <- fit_ability(
  response_data = responses,
  item_loadings = loadings,
  threshold_left = thresh_left,
  threshold_right = thresh_right,
  method = "nloptr",
  algorithm = "NLOPT_LD_LBFGS",
  maxeval = 1000
)

# Compare results
cor(predict(fit_optim), predict(fit_nloptr))
```

---

## Success Criteria

### Phase Completion Checklist

- [ ] Phase 0: Package structure initialized, dependencies configured
- [ ] Phase 1: Stan model compiles and is accessible
- [ ] Phase 2: Data preparation functions work with validation
- [ ] Phase 3: All Stan estimation methods (MAP, VB, MCMC) functional
- [ ] Phase 4: BridgeStan integration works with multiple optimizers
- [ ] Phase 5: Main user interface is intuitive and documented
- [ ] Phase 6: S3 methods provide useful output
- [ ] Phase 7: Documentation and vignettes are complete

### Quality Metrics

- [ ] All unit tests pass
- [ ] Integration tests pass
- [ ] Code coverage > 80%
- [ ] `R CMD check` passes with 0 errors, 0 warnings
- [ ] Documentation is comprehensive
- [ ] Vignettes provide clear usage examples
- [ ] Performance is acceptable for realistic datasets

### Final Deliverables

1. **Working R package** that can be installed and loaded
2. **Comprehensive test suite** with >80% coverage
3. **Documentation** for all exported functions
4. **Vignettes** demonstrating key use cases
5. **README** with installation and quick start
6. **NEWS.md** documenting version changes

## Notes for Implementation

### Dependencies to Monitor
- **Stan**: Both `rstan` and `cmdstanr` should be supported
- **BridgeStan**: Optional dependency, check with `has_bridgestan()`
- **Optimization packages**: All optional, graceful handling if missing

### Common Pitfalls
1. **Stan compilation**: Ensure models compile on package install
2. **Infinite threshold handling**: Test edge cases thoroughly
3. **Large datasets**: Monitor memory usage and performance
4. **Convergence**: Provide diagnostic tools for MCMC
5. **Parameter constraints**: Ensure tau > 0, proper correlation matrix

### Extension Points
- Additional estimation methods beyond Laplace approximation
- Diagnostic plots and tools
- Model comparison functions (DIC, WAIC, LOO-CV)
- Simulation functions for power analysis
- Support for polytomous response models beyond graded response

---

## Timeline Estimate

- **Phase 0**: 1-2 days
- **Phase 1**: 2-3 days
- **Phase 2**: 3-4 days
- **Phase 3**: 4-5 days
- **Phase 4**: 3-4 days
- **Phase 5**: 2-3 days
- **Phase 6**: 3-4 days
- **Phase 7**: 3-5 days

**Total**: 3-4 weeks for full implementation and testing