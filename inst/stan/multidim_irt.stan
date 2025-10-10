// Multidimensional IRT Model with Person Covariates
//
// This Stan model implements a flexible multidimensional IRT model with:
// - K latent dimensions (unidimensional when K=1)
// - Response-level factor loadings (L x K matrix)
// - Person-level covariates affecting each dimension
// - Correlated dimensions via Cholesky-factored correlation matrix
// - Ordered categorical responses with threshold boundaries
// - Infinite threshold handling (coded as -999)

data {
  // Infinite threshold code
  real inf;                    // Code for infinity (typically -999)

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
  // Check for boundary thresholds and convert infinities
  vector[L] dL_finite;
  vector[L] dR_finite;
  array[L] int is_dL_inf;
  array[L] int is_dR_inf;

  for (l in 1:L) {
    // Check if thresholds are infinite (coded as inf)
    is_dL_inf[l] = (dL[l] == inf) ? 1 : 0;
    is_dR_inf[l] = (dR[l] == inf) ? 1 : 0;

    // Replace infinite values with finite placeholders for computation
    dL_finite[l] = is_dL_inf[l] ? -10.0 : dL[l];
    dR_finite[l] = is_dR_inf[l] ? 10.0 : dR[l];
  }
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

    if (is_dL_inf[l] && is_dR_inf[l]) {
      // Both boundaries infinite: probability = 1 (degenerate)
      pr = 1.0;
    } else if (is_dL_inf[l]) {
      // Left boundary is -infinity: P = Phi(eta + dR)
      pr = Phi(eta_l + dR_finite[l]);
    } else if (is_dR_inf[l]) {
      // Right boundary is +infinity: P = 1 - Phi(eta + dL)
      pr = 1.0 - Phi(eta_l + dL_finite[l]);
    } else {
      // Both boundaries finite: P = Phi(eta + dR) - Phi(eta + dL)
      pr = Phi(eta_l + dR_finite[l]) - Phi(eta_l + dL_finite[l]);
    }

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

    if (is_dL_inf[l] && is_dR_inf[l]) {
      pr = 1.0;
    } else if (is_dL_inf[l]) {
      pr = Phi(eta_l + dR_finite[l]);
    } else if (is_dR_inf[l]) {
      pr = 1.0 - Phi(eta_l + dL_finite[l]);
    } else {
      pr = Phi(eta_l + dR_finite[l]) - Phi(eta_l + dL_finite[l]);
    }

    log_lik += wgt[person] * log(pr);
  }
}
