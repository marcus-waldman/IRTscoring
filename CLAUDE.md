# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**IRTscoring** is an R package for flexible Item Response Theory (IRT) scoring that supports both unidimensional and multidimensional models with multiple estimation methods.

### Key Features
- Multidimensional IRT with full factor loadings matrix (L × K)
- Multiple estimation methods: MAP, Laplace Approximation, Variational Bayes (Pathfinder), MCMC
- Stan-based estimation with optional BridgeStan integration for alternative optimizers
- Person-level covariates and dimension correlations
- **0-indexed responses required** (responses must start at 0, not 1)
- Automatic threshold validation

## Critical Design Decisions

### Response Indexing
**IMPORTANT**: All item responses must be 0-indexed:
- Valid responses: 0, 1, 2, ... (NOT 1, 2, 3, ...)
- For 3 categories: use 0, 1, 2
- For binary items: use 0, 1

### Threshold Specification
Thresholds define category boundaries:
- For response value k (0-indexed), at least k thresholds must be provided
- Categories are defined by intervals: (-∞, t₁], (t₁, t₂], ..., (tₙ, ∞)
- The package validates sufficient thresholds for all observed responses

### Score Decomposition
The package decomposes ability scores (θ) into components for each dimension k:
- **Full abilities**: θ[n,k] = β₀[k] + Σⱼ X[n,j]·b[k,j] + η[n,k]
- **Individual effects**: η (N × K matrix of person-specific deviations)
- **Covariate effects**: β₀ + X·b (N × K matrix predicted from covariates)
- Each dimension has its own intercept β₀[k] and covariate coefficients b[k,j]

## Development Commands

### Package Setup
```r
# Load package for development
devtools::load_all()

# Build and check
devtools::check()

# Run tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-data-preparation.R")

# Generate documentation
devtools::document()

# Install locally
devtools::install()
```

### Stan Model Compilation
Stan models are located in `inst/stan/multidim_irt.stan`. The model compiles on package installation. During development:
```r
# Manually compile Stan model
model <- cmdstanr::cmdstan_model("inst/stan/multidim_irt.stan")
```

### Testing Against mirt Package
The package should match results from the `mirt` package:
```r
# Run comparison tests
testthat::test_file("tests/testthat/test-mirt-comparison.R")
```

**Note**: When comparing with mirt, convert from mirt's 1-indexed responses to IRTscoring's 0-indexed responses.

## Architecture

### Package Structure
```
IRTscoring/
├── R/
│   ├── fscores.R              # User-facing wide format function
│   ├── fit_ability.R          # Main long format function
│   ├── prepare_data.R         # Data validation/formatting
│   ├── stan_estimation.R      # MAP, VB, MCMC, Laplace methods
│   ├── bridgestan_estimation.R # Alternative optimizers
│   ├── methods.R              # S3 methods (print, summary, coef, predict)
│   ├── utilities.R            # Helper functions
│   └── zzz.R                  # Package initialization
├── inst/stan/
│   └── multidim_irt.stan      # Stan model
├── tests/testthat/
│   ├── test-data-preparation.R
│   ├── test-stan-estimation.R
│   ├── test-mirt-comparison.R
│   └── test-methods.R
└── vignettes/
```

### Function Hierarchy

**User-facing functions**:
- `fscores()`: Wide format wrapper (easiest for users)
- `fit_ability()`: Long format main function
- `fit_unidim_ability()`: Unidimensional convenience wrapper

**Internal functions**:
- `prepare_irt_data()`: Converts to Stan data format
- `fit_irt_map()`, `fit_irt_laplace()`, `fit_irt_vb()`, `fit_irt_mcmc()`: Estimation methods
- `fit_irt_optim()`, `fit_irt_nloptr()`: BridgeStan-based optimizers

**S3 methods** on `irt_fit` class:
- `print()`: Summary output
- `summary()`: Detailed results
- `coef()`: Parameter extraction
- `predict()`: Ability scores with `type` argument

### Data Flow

1. **Wide format input** → `fscores()` → converts to long format
2. **Long format** → `prepare_irt_data()` → Stan data list
3. **Stan data** → `fit_irt_*()` → estimation
4. **Estimation result** → `irt_fit` object
5. **irt_fit** → S3 methods → user results

### Stan Data Structure

The prepared data list for Stan includes:
- `inf`: Code for infinity (-999)
- `L`: Number of response observations
- `N`: Number of persons
- `J`: Number of covariates
- `K`: Number of dimensions
- `pid`: Person IDs (length L)
- `dL`, `dR`: Left/right thresholds (length L)
- `wgt`: Person weights (length N)
- `X`: Covariate matrix (N × J)
- `A`: Loading matrix (L × K)

## Implementation Status

**Project Complete**: All phases implemented and tested.

The project followed an 8-phase implementation plan:

- [x] **Phase 0**: Project setup (package structure, dependencies)
- [x] **Phase 1**: Stan model integration
- [x] **Phase 2**: Data preparation functions
- [x] **Phase 3**: Stan-based estimation (MAP, Laplace, VB, MCMC)
- [x] **Phase 4**: BridgeStan integration (optional optimizers)
- [x] **Phase 5**: Main user interface (`fscores()`, `fit_ability()`)
- [x] **Phase 6**: S3 methods and utilities
- [x] **Phase 7**: Documentation and package finalization

The package is production-ready with comprehensive documentation and testing.

## R Environment

- R.exe and Rscript.exe are located in `C:\Program Files\R\R-4.5.1\bin`
- The `%||%` operator does not exist in base R - use if statements instead
- Always ask before taking simpler approaches when implementation barriers are encountered

## Key Dependencies

**Required**:
- `rstan` (>= 2.26.0) or `cmdstanr`: Stan estimation
- `checkmate`: Input validation

**Suggested**:
- `bridgestan`: Alternative optimizers (experimental)
- `optimx`, `trustOptim`, `nloptr`: Optimization algorithms
- `mirt`: Validation comparisons
- `MASS`: Multivariate normal sampling (Laplace approximation)
- `testthat` (>= 3.0.0): Testing
- `devtools`: Development workflow

## Testing Strategy

1. **Unit tests**: Test individual functions in isolation
2. **Integration tests**: Full workflow from data to results
3. **Comparison tests**: Validate against `mirt` package
4. **Performance tests**: Ensure reasonable computation time
5. **Validation tests**: Match published examples

### Critical Test Cases
- 0-indexed response validation
- Threshold sufficiency validation
- Score decomposition correctness (abilities = covariate_effects + individual_effects)
- Agreement with mirt::fscores() for MAP estimates (correlation > 0.95, RMSE < 0.2)
