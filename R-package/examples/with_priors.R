# Example: HiSpa Analysis with Cluster Initialization

library(HiSpaR)

# Load contact matrix
contact_matrix <- as.matrix(read.table("contacts_chr1.txt"))

# Run analysis with cluster-based initialization
# This runs MCMC on individual clusters first, then assembles them
result <- hispa_analyze(
  contact_matrix = contact_matrix,
  output_dir = "output_cluster_init",
  mcmc_iterations = 6000,
  mcmc_burn_in = 1000,
  num_clusters = 5,  # Specify number of clusters
  use_cluster_init = TRUE,  # Use cluster-based initialization
  cluster_init_iterations = 1000,  # MCMC iterations per cluster
  cluster_initial_sd = 0.1,  # SD for cluster MCMC
  save_samples = TRUE,  # Save MCMC trace
  sample_interval = 50,
  verbose = TRUE
)

# View results
print(result)
summary(result)

# Extract results
final_positions <- result$final_positions
initial_positions <- result$initial_positions
ll_trace <- result$log_likelihood_trace

# Plot convergence
plot(ll_trace, type = "l",
     xlab = "MCMC Iteration",
     ylab = "Log-Likelihood",
     main = "MCMC Convergence with Cluster Initialization")
abline(v = 1000, col = "red", lty = 2)

# Compare with random initialization
result_random <- hispa_analyze(
  contact_matrix = contact_matrix,
  output_dir = "output_random_init",
  mcmc_iterations = 6000,
  mcmc_burn_in = 1000,
  num_clusters = 5,
  use_cluster_init = FALSE,  # Random initialization
  verbose = FALSE
)

# Compare convergence
par(mfrow = c(1, 2))
plot(result$log_likelihood_trace, type = "l",
     main = "Cluster Initialization",
     xlab = "Iteration", ylab = "Log-Likelihood")
plot(result_random$log_likelihood_trace, type = "l",
     main = "Random Initialization",
     xlab = "Iteration", ylab = "Log-Likelihood")

cat("\nFinal log-likelihood (cluster init):", tail(result$log_likelihood_trace, 1), "\n")
cat("Final log-likelihood (random init):", tail(result_random$log_likelihood_trace, 1), "\n")
