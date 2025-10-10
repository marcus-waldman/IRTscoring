# IRTscoring 0.1.0

## Major Features

* **Multiple Estimation Methods**: MAP, Laplace Approximation, Variational Bayes (ADVI and Pathfinder), Full Bayesian MCMC, and optimization-based methods (optim, nloptr)

* **Flexible Stan Integration**: Supports both `cmdstanr` (recommended) and `rstan` backends with automatic detection

* **Multidimensional IRT**: Native support for K-dimensional latent spaces with correlated dimensions (Omega correlation matrix)

* **Response-Level Factor Loadings**: Full L × K loading matrices where L = number of responses, K = number of dimensions

* **Categorical Response Models**: Threshold-based ordered categorical responses with automatic validation

* **Flexible Data Formats**:
  - Long format via `fit_ability()` (person-item-response structure)
  - Wide format via `fscores()` (person × items matrix, mirt-compatible)
  - Simplified unidimensional interface via `fit_unidim_ability()`

* **Score Decomposition**: Ability scores decomposed as θ = β₀ + X·b + η
  - Fixed effects (intercepts + covariate effects)
  - Random effects (individual deviations)
  - Utility functions: `get_covariate_effects()`, `get_individual_effects()`

* **Comprehensive S3 Methods**:
  - `print.irt_fit()`: Clean summary output
  - `summary.irt_fit()`: Detailed model information
  - `coef.irt_fit()`: Flexible parameter extraction
  - `predict.irt_fit()`: Ability score predictions

* **Experimental BridgeStan Support**: Alternative optimization interface with gradient-based methods

## Data Requirements

* **CRITICAL**: Response categories MUST be 0-indexed (0, 1, 2, ... not 1, 2, 3, ...)
* Automatic validation detects 1-indexed responses and provides clear error messages
* Threshold validation ensures sufficient thresholds for observed response categories

## Backend Support

* **cmdstanr** (recommended): Faster compilation, modern algorithms (Pathfinder), better diagnostics
* **rstan**: Fallback support for compatibility
* **bridgestan** (experimental): Gradient-based optimization methods

## Package Architecture

Complete implementation across 7 phases:

1. **Phase 0**: Project setup and structure
2. **Phase 1**: Stan model integration (multidim_irt.stan)
3. **Phase 2**: Data preparation and validation functions
4. **Phase 3**: Stan-based estimation (MAP, Laplace, VB, MCMC)
5. **Phase 4**: BridgeStan integration (optim, nloptr)
6. **Phase 5**: Main user interfaces (fit_ability, fscores, fit_unidim_ability)
7. **Phase 6**: S3 methods and utility functions
8. **Phase 7**: Documentation and vignettes

## Testing

* Comprehensive test suite with 100+ tests
* Tests for all estimation methods across backends
* Score decomposition validation
* Edge case handling (single person, multidimensional, etc.)

## Documentation

* Complete roxygen2 documentation for all exported functions
* Comprehensive README with Quick Start examples
* Examples for unidimensional, multidimensional, and wide format usage
* Estimation methods comparison table
* Score decomposition explanation
