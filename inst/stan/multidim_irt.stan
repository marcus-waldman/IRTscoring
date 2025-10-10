// Multidimensional IRT Model with Person Covariates
//
// This Stan model implements a flexible multidimensional IRT model with:
// - K latent dimensions (unidimensional when K=1)
// - Response-level factor loadings (L x K matrix)
// - Person-level covariates affecting each dimension
// - Correlated dimensions via Cholesky-factored correlation matrix
// - Ordered categorical responses with threshold boundaries
// - Infinite threshold handling: -Inf→-999, +Inf→+999 (done in R)

data {
  // Note: Infinite thresholds are handled in R before passing to Stan
  // -Inf is mapped to -999, +Inf is mapped to +999
  real inf;                    // Placeholder (999), not actually used

  // Dimensions
  int<lower=1> L;              // Number of response observations
  int<lower=1> N;              // Number of persons
  int<lower=1> J;              // Number of covariates (including intercept)
  int<lower=1> K;              // Number of latent dimensions

  // Response data
  array[L] int<lower=1,upper=N> pid;  // Person ID for each response
  vector[L] dL;                // Left threshold for each response
  vector[L] dR;                // Right threshold for each response

  // Person-level data
  vector[N] wgt;               // Person weights
  matrix[N, J] X;              // Covariate matrix (N x J)

  // Item parameters (known/fixed)
  matrix[L, K] A;              // Factor loadings matrix (L x K)
}

transformed data {
  // No transformation needed - infinities already handled in R
  // -Inf mapped to -999, +Inf mapped to +999
}

parameters {
  // Dimension intercepts
  vector[K] beta0;

  // Covariate effects on each dimension
  matrix[K, J] b;

  // Individual effects (deviations from covariate predictions)
  matrix[N, K] zeta;

  // Dimension correlation structure (Cholesky factor)
  cholesky_factor_corr[K] L_Omega;

  // Dimension standard deviations
  vector<lower=0>[K] tau;
}

transformed parameters {
  // Full ability parameters: theta = beta0 + X*b + eta
  // where eta = zeta * diag(tau) * L_Omega'
  matrix[N, K] eta;
  matrix[N, K] theta;

  // Compute eta from standardized zeta
  eta = zeta * diag_pre_multiply(tau, L_Omega)';

  // Compute full theta (ability scores)
  theta = rep_matrix(beta0', N) + X * b' + eta;
}

model {
  // Priors

  // Dimension intercepts: weakly informative
  beta0 ~ normal(0, 2);

  // Covariate effects: weakly informative
  to_vector(b) ~ normal(0, 2);

  // Individual effects: standard normal (before scaling)
  to_vector(zeta) ~ std_normal();

  // Dimension correlations: LKJ prior
  L_Omega ~ lkj_corr_cholesky(1);

  // Dimension standard deviations: half-normal
  tau ~ normal(0, 1);

  // Likelihood
  for (l in 1:L) {
    int person = pid[l];
    real eta_l;
    real pr;

    // Compute linear predictor: eta = sum(A[l,k] * theta[person,k])
    eta_l = dot_product(A[l,], theta[person,]);

    // Compute probability for this response category
    // P(response in category) = Phi(eta + dR) - Phi(eta + dL)
    // Since -Inf→-999 and +Inf→+999, Phi will naturally handle boundaries:
    // Phi(eta - 999) ≈ 0 and Phi(eta + 999) ≈ 1
    pr = Phi(eta_l + dR[l]) - Phi(eta_l + dL[l]);

    // Add to log likelihood with person weight
    target += wgt[person] * log(pr);
  }
}

generated quantities {
  // Correlation matrix of dimensions
  matrix[K, K] Omega;

  // Reconstruct correlation matrix from Cholesky factor
  Omega = L_Omega * L_Omega';

  // Log-likelihood for model comparison
  real log_lik = 0.0;

  for (l in 1:L) {
    int person = pid[l];
    real eta_l;
    real pr;

    eta_l = dot_product(A[l,], theta[person,]);
    pr = Phi(eta_l + dR[l]) - Phi(eta_l + dL[l]);

    log_lik += wgt[person] * log(pr);
  }
}
