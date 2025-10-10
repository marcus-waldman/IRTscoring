# IRTscoring: Flexible IRT Scoring in R

<!-- badges: start -->
<!-- badges: end -->

**IRTscoring** provides flexible estimation methods for unidimensional and multidimensional Item Response Theory (IRT) models with support for multiple estimation approaches including MAP, Laplace Approximation, Variational Bayes (including Pathfinder), and full Bayesian MCMC via Stan.

## Key Features

- **Multidimensional IRT**: Support for K = 1 (unidimensional) through any positive integer (multidimensional)
- **Multiple Estimation Methods**: MAP, Laplace Approximation, Variational Bayes, MCMC, and alternative optimizers via BridgeStan
- **Flexible Model Structure**: Response-level factor loadings, person-level covariates, and dimension correlations
- **User-Friendly Interface**: `fscores()` function accepts data in wide format (similar to mirt)
- **0-Indexed Responses**: Responses must start at 0 (not 1)
- **Automatic Validation**: Threshold validation ensures sufficient thresholds for observed responses

## Installation

```r
# Install from GitHub (once available)
# devtools::install_github("marcus-waldman/IRTscoring")
```

## Quick Start

### Unidimensional IRT with Wide Format Data

```r
library(IRTscoring)

# IMPORTANT: Item responses must be 0-indexed (start at 0, not 1)
wide_data <- data.frame(
  person_id = 1:100,
  age = rnorm(100, 30, 10),
  item1 = sample(0:2, 100, replace = TRUE),  # 0, 1, 2 (NOT 1, 2, 3)
  item2 = sample(0:2, 100, replace = TRUE)
)

# Define item parameters
item_params <- list(
  item1 = list(loadings = c(1.2), thresholds = c(-1, 0.5)),
  item2 = list(loadings = c(1.5), thresholds = c(-0.5, 1.0))
)

# Fit model with MAP estimation
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

### Multidimensional IRT

```r
# For K=2 dimensions
item_params_2d <- list(
  item_Q1 = list(
    loadings = c(1.2, 0.1),      # Loads on both dimensions
    thresholds = c(-1.0, 0.5)    # 2 thresholds for 3 categories
  ),
  item_Q2 = list(
    loadings = c(0.1, 1.3),      # Loads primarily on dimension 2
    thresholds = c(1.0)          # 1 threshold for 2 categories (binary)
  )
)

fit_2d <- fscores(
  data = wide_data,
  item_params = item_params_2d,
  person_id_col = "person_id",
  method = "map"
)

# Get dimension correlations
omega <- get_dimension_correlations(fit_2d)
```

## Development Status

**Current Status**: Phase 0 - Project Setup Complete

This package is currently under active development. The implementation follows a structured 6-phase plan:

- [x] Phase 0: Project Setup
- [ ] Phase 1: Stan Model Integration
- [ ] Phase 2: Data Preparation Functions
- [ ] Phase 3: Stan-based Estimation Functions
- [ ] Phase 4: BridgeStan Integration
- [ ] Phase 5: Main User Interface
- [ ] Phase 6: S3 Methods and Utilities

## License

MIT License - see LICENSE file for details

## Author

Marcus Waldman
