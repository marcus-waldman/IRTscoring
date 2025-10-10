# Manual testing script for estimation functions
# This script compiles the Stan model once and reuses it across tests

library(devtools)
load_all()

cat("=== Testing Phase 3: Stan Estimation Functions ===\n\n")

# Create simple test data
cat("Creating test data...\n")
response_data <- data.frame(
  pid = c(1, 1, 2, 2, 3, 3, 3),
  iid = c(1, 2, 1, 2, 1, 2, 3)
)

loadings <- matrix(c(1, 1.5, 1, 1.5, 2, 1.2, 0.8), ncol = 1)
thresh_L <- c(-Inf, -1, -Inf, -1, -Inf, -2, 0)
thresh_R <- c(-1, Inf, -1, Inf, -2, 0, Inf)

stan_data <- prepare_irt_data(response_data, loadings, thresh_L, thresh_R)

cat("Data prepared: N =", stan_data$N, ", L =", stan_data$L, ", K =", stan_data$K, "\n\n")

# Test 1: Compile model (will be cached)
cat("=== Test 1: Compiling Stan model with cmdstanr ===\n")
model <- get_stan_model(backend = "cmdstanr")
cat("Model compiled successfully!\n")
cat("Backend:", get_stan_backend(), "\n")
cat("Model cached:", is_stan_model_compiled(), "\n\n")

# Test 2: MAP estimation
cat("=== Test 2: MAP Estimation ===\n")
fit_map <- fit_irt_map(stan_data, backend = "cmdstanr")
cat("MAP estimation completed!\n")
cat("Method:", fit_map$method, "\n")
cat("Backend:", fit_map$backend, "\n")
cat("Theta dimensions:", paste(dim(fit_map$theta), collapse = " x "), "\n")
cat("Theta estimates:\n")
print(fit_map$theta)
cat("Beta0:", fit_map$beta0, "\n")
cat("Tau:", fit_map$tau, "\n")
cat("Omega:\n")
print(fit_map$Omega)
cat("Log-likelihood:", fit_map$log_lik, "\n\n")

# Test 3: Check that model is reused (should be fast)
cat("=== Test 3: Reusing compiled model ===\n")
start_time <- Sys.time()
fit_map2 <- fit_irt_map(stan_data, backend = "cmdstanr")
elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
cat("Second MAP estimation completed in", round(elapsed, 3), "seconds\n")
cat("(Should be fast since model is cached)\n\n")

# Test 4: Multidimensional model
cat("=== Test 4: Multidimensional Model ===\n")
response_data_2d <- data.frame(
  pid = c(1, 1, 2, 2),
  iid = c(1, 2, 1, 2)
)

loadings_2d <- matrix(c(1, 0.5, 1, 0.5,    # Dim 1
                        0.5, 1, 0.5, 1),    # Dim 2
                      ncol = 2)
thresh_L_2d <- c(-Inf, -1, -Inf, -1)
thresh_R_2d <- c(-1, Inf, -1, Inf)

stan_data_2d <- prepare_irt_data(response_data_2d, loadings_2d,
                                  thresh_L_2d, thresh_R_2d)

fit_2d <- fit_irt_map(stan_data_2d, backend = "cmdstanr")
cat("Multidimensional MAP completed!\n")
cat("Theta dimensions:", paste(dim(fit_2d$theta), collapse = " x "), "\n")
cat("Theta estimates:\n")
print(fit_2d$theta)
cat("Correlation matrix (Omega):\n")
print(fit_2d$Omega)
cat("\n")

# Test 5: VB estimation (fast alternative)
cat("=== Test 5: Variational Bayes Estimation ===\n")
fit_vb <- fit_irt_vb(stan_data, backend = "cmdstanr",
                     algorithm = "meanfield", draws = 500)
cat("VB estimation completed!\n")
cat("Method:", fit_vb$method, "\n")
cat("Theta estimates:\n")
print(fit_vb$theta)
cat("\nComparing MAP vs VB theta estimates:\n")
cat("Mean absolute difference:", mean(abs(fit_map$theta - fit_vb$theta)), "\n\n")

# Test 6: MCMC with very few iterations (for testing only)
cat("=== Test 6: MCMC Estimation (minimal iterations) ===\n")
fit_mcmc <- fit_irt_mcmc(stan_data, backend = "cmdstanr",
                         chains = 2, iter = 200, warmup = 100)
cat("MCMC completed!\n")
cat("Method:", fit_mcmc$method, "\n")
cat("Theta estimates:\n")
print(fit_mcmc$theta)
cat("\n")

cat("=== All tests completed successfully! ===\n")
