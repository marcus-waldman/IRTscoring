# IRTscoring

<!-- badges: start -->
<!-- badges: end -->

Flexible Item Response Theory (IRT) scoring with multiple estimation methods.

## Overview

**IRTscoring** provides comprehensive tools for estimating person abilities in unidimensional and multidimensional IRT models. The package supports multiple estimation backends (Stan, BridgeStan) and methods (MAP, Laplace, VB, MCMC, optimization), making it suitable for both research and operational scoring applications.

## Features

- **Multiple Estimation Methods**: MAP, Laplace Approximation, Variational Bayes (including Pathfinder), Full Bayesian MCMC, and optimization-based methods
- **Flexible Stan Integration**: Supports both `cmdstanr` (recommended) and `rstan` backends
- **Experimental BridgeStan Support**: Alternative optimization via `nloptr` and `optim` with gradient-based methods
- **Multidimensional IRT**: Native support for correlated latent dimensions
- **Response-Level Factor Loadings**: Full flexibility with L × K loading matrices (L = responses, K = dimensions)
- **Categorical Responses**: Handles ordered categorical data with threshold-based models
- **Score Decomposition**: Ability scores decomposed into covariate effects and individual deviations
- **Flexible Data Formats**: Long format (person-item-response) and wide format (person × items matrix)

## Installation

```r
# Install from GitHub
devtools::install_github("username/IRTscoring")

# For cmdstanr support (recommended):
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

# For rstan support:
install.packages("rstan")

# For BridgeStan support (experimental, optional):
# Follow instructions at https://github.com/roualdes/bridgestan
```

## ⚠️ CRITICAL: Response Categories Must Be 0-Indexed

**Response categories MUST start at 0, not 1.** This is a fundamental requirement of the package.

- ✅ **Correct**: Binary responses coded as `0, 1`
- ✅ **Correct**: 3-category responses coded as `0, 1, 2`
- ❌ **Incorrect**: Binary responses coded as `1, 2`
- ❌ **Incorrect**: 3-category responses coded as `1, 2, 3`

The package will detect 1-indexed responses and throw an error. Always check your data coding before estimation.

## Quick Start

### Unidimensional IRT Scoring

```r
library(IRTscoring)

# Create response data (long format)
# Categories MUST be 0-indexed!
response_data <- data.frame(
  pid = c(1, 1, 1, 2, 2, 2, 3, 3, 3),  # Person IDs
  iid = c(1, 2, 3, 1, 2, 3, 1, 2, 3),  # Item IDs (response-level)
  response = c(0, 1, 0, 1, 1, 0, 0, 0, 1)  # 0-indexed!
)

# Response-level factor loadings (L × K matrix)
# Each response gets its own loading
item_loadings <- matrix(c(
  1.2,  # Response 1 loading
  1.5,  # Response 2 loading
  1.8   # Response 3 loading
), ncol = 1)

# Threshold parameters for ordered categories
# For binary: one threshold separates categories 0 and 1
threshold_left <- c(-Inf, -Inf, -Inf)
threshold_right <- c(-0.5, 0.0, 0.5)

# Fit model using MAP estimation
fit <- fit_ability(
  response_data = response_data,
  item_loadings = item_loadings,
  threshold_left = threshold_left,
  threshold_right = threshold_right,
  method = "map",
  backend = "cmdstanr"  # or "rstan"
)

# View results
print(fit)
summary(fit)

# Extract ability estimates
theta <- predict(fit)
print(theta)
```

### Multidimensional IRT Scoring

```r
# Response data with 2 dimensions
response_data <- data.frame(
  pid = c(1, 1, 1, 2, 2, 2),
  iid = c(1, 2, 3, 1, 2, 3),
  response = c(0, 1, 0, 1, 0, 1)
)

# Response-level loadings: L × K matrix (6 responses × 2 dimensions)
item_loadings <- matrix(c(
  1.0, 0.2,  # Response 1: loads on both dimensions
  1.2, 0.3,  # Response 2
  0.8, 0.1,  # Response 3
  1.0, 0.2,  # Response 4 (person 2)
  1.2, 0.3,  # Response 5
  0.8, 0.1   # Response 6
), ncol = 2, byrow = TRUE)

threshold_left <- rep(-Inf, 6)
threshold_right <- rep(-0.5, 6)

# Fit multidimensional model
fit <- fit_ability(
  response_data = response_data,
  item_loadings = item_loadings,
  threshold_left = threshold_left,
  threshold_right = threshold_right,
  method = "map"
)

# Get dimension correlations
Omega <- coef(fit, pars = "Omega")
print(Omega)  # 2×2 correlation matrix

# Extract 2D ability estimates
theta <- predict(fit)  # N × 2 matrix
print(theta)
```

### Wide Format (Factor Scores Interface)

```r
# Wide format data: persons × items
wide_data <- data.frame(
  person_id = 1:3,
  item1 = c(0, 1, 0),  # Still 0-indexed!
  item2 = c(1, 1, 0),
  item3 = c(0, 0, 1)
)

# Item parameters (one entry per item)
item_params <- list(
  item1 = list(loadings = c(1.2), thresholds = c(-0.5)),
  item2 = list(loadings = c(1.5), thresholds = c(0.0)),
  item3 = list(loadings = c(1.8), thresholds = c(0.5))
)

# Simplified interface for wide data
scores <- fscores(
  data = wide_data,
  item_params = item_params,
  method = "map"
)

print(scores)
```

## Estimation Methods Comparison

| Method | Speed | Uncertainty | Use Case |
|--------|-------|-------------|----------|
| **MAP** | ⚡⚡⚡ Fastest | Point estimate only | Operational scoring, large-scale |
| **Laplace** | ⚡⚡ Fast | Asymptotic SE | Quick uncertainty estimates |
| **VB (ADVI)** | ⚡⚡ Fast | Approximate posterior | Balance speed/uncertainty |
| **VB (Pathfinder)** | ⚡⚡⚡ Fastest Bayesian | Approximate posterior | Modern alternative to ADVI |
| **MCMC** | ⚡ Slow | Full posterior | Research, small samples |
| **optim** | ⚡⚡⚡ Very fast | Point estimate only | Experimental, requires BridgeStan |
| **nloptr** | ⚡⚡ Fast | Point estimate only | Experimental, requires BridgeStan |

```r
# Different estimation methods
fit_map <- fit_ability(..., method = "map")
fit_laplace <- fit_ability(..., method = "laplace")
fit_vb <- fit_ability(..., method = "vb")
fit_mcmc <- fit_ability(..., method = "mcmc", mcmc_chains = 4, mcmc_iter = 2000)

# Using specific Stan VB algorithms
fit_pathfinder <- fit_ability(..., method = "vb", vb_algorithm = "pathfinder")
fit_advi <- fit_ability(..., method = "vb", vb_algorithm = "meanfield")

# Experimental BridgeStan methods (requires bridgestan package)
fit_optim <- fit_ability(..., method = "optim", backend = "bridgestan")
fit_nloptr <- fit_ability(..., method = "nloptr", backend = "bridgestan")
```

## Score Decomposition

Ability scores are decomposed into fixed and random components:

**θ = β₀ + X·b + η**

Where:
- **θ** = Total ability scores (what you observe)
- **β₀** = Dimension intercepts (population means)
- **X·b** = Covariate effects (if covariates included)
- **η** = Individual-specific deviations (person effects)

```r
# Extract components
abilities <- predict(fit)                    # θ (total scores)
cov_effects <- get_covariate_effects(fit)    # β₀ + X·b
ind_effects <- get_individual_effects(fit)   # η

# Verify decomposition
all.equal(abilities, cov_effects + ind_effects)  # TRUE
```

## Stan Backends

The package supports both `cmdstanr` (recommended) and `rstan`:

```r
# Using cmdstanr (recommended - faster, more features)
fit <- fit_ability(..., backend = "cmdstanr")

# Using rstan (fallback)
fit <- fit_ability(..., backend = "rstan")

# Automatic detection
fit <- fit_ability(...)  # Uses cmdstanr if available, else rstan
```

**Recommendation**: Use `cmdstanr` for better performance and access to modern algorithms (Pathfinder).

## BridgeStan Support (Experimental)

BridgeStan provides an alternative interface for optimization-based methods:

```r
# Requires: devtools::install_github("roualdes/bridgestan")

# Check availability
has_bridgestan()  # TRUE if installed

# Use with gradient-based optimizers
fit_optim <- fit_ability(
  ...,
  method = "optim",
  backend = "bridgestan",
  optim_method = "BFGS"
)

fit_nloptr <- fit_ability(
  ...,
  method = "nloptr",
  backend = "bridgestan",
  nloptr_algorithm = "NLOPT_LD_LBFGS"
)
```

**Status**: Experimental. Stan-based methods (MAP, Laplace, VB, MCMC) are recommended for production use.

## S3 Methods

```r
# Print summary
print(fit)

# Detailed summary with all parameters
summary(fit)

# Extract specific parameters
coef(fit, pars = "beta0")       # Intercepts
coef(fit, pars = "tau")         # SDs
coef(fit, pars = "Omega")       # Dimension correlations
coef(fit, pars = "theta")       # Ability scores

# Predict abilities (alias for coef(fit, "theta"))
predict(fit)
predict(fit, type = "abilities")
```

## Development Status

This package is complete and ready for use. All phases implemented:

- [x] Phase 0: Project Setup
- [x] Phase 1: Stan Model Integration
- [x] Phase 2: Data Preparation Functions
- [x] Phase 3: Stan-based Estimation Functions
- [x] Phase 4: BridgeStan Integration
- [x] Phase 5: Main User Interface
- [x] Phase 6: S3 Methods and Utilities
- [x] Phase 7: Documentation and Vignettes

## License

MIT License - see LICENSE file for details

## Author

Marcus Waldman
