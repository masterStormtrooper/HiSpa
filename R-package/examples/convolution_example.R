# Example: Comparing Different Settings

library(HiSpaR)

# Load contact matrix
contact_matrix <- as.matrix(read.table("contact_matrix.txt"))

# Test 1: Random initialization
cat("Running with random initialization...\n")
result_random <- hispa_analyze(
  contact_matrix = contact_matrix,
  output_dir = "output_random",
  mcmc_iterations = 3000,
  mcmc_burn_in = 500,
  use_cluster_init = FALSE,
  verbose = FALSE
)

# Test 2: Cluster initialization
cat("Running with cluster initialization...\n")
result_cluster <- hispa_analyze(
  contact_matrix = contact_matrix,
  output_dir = "output_cluster",
  mcmc_iterations = 3000,
  mcmc_burn_in = 500,
  use_cluster_init = TRUE,
  cluster_init_iterations = 1000,
  verbose = FALSE
)

# Test 3: Different number of clusters
cat("Running with more clusters...\n")
result_more_clusters <- hispa_analyze(
  contact_matrix = contact_matrix,
  output_dir = "output_more_clusters",
  mcmc_iterations = 3000,
  mcmc_burn_in = 500,
  num_clusters = 10,
  use_cluster_init = TRUE,
  cluster_init_iterations = 500,
  verbose = FALSE
)

# Compare results
cat("\n=== Comparison ===\n")
cat(sprintf("Random init - Final LL: %.2f\n", 
            tail(result_random$log_likelihood_trace, 1)))
cat(sprintf("Cluster init - Final LL: %.2f\n", 
            tail(result_cluster$log_likelihood_trace, 1)))
cat(sprintf("More clusters - Final LL: %.2f\n", 
            tail(result_more_clusters$log_likelihood_trace, 1)))

# Plot convergence comparison
par(mfrow = c(2, 2))

plot(result_random$log_likelihood_trace, type = "l",
     main = "Random Initialization",
     xlab = "Iteration", ylab = "Log-Likelihood")

plot(result_cluster$log_likelihood_trace, type = "l",
     main = "Cluster Initialization",
     xlab = "Iteration", ylab = "Log-Likelihood")

plot(result_more_clusters$log_likelihood_trace, type = "l",
     main = "More Clusters",
     xlab = "Iteration", ylab = "Log-Likelihood")

# Compare all on same plot
plot(result_random$log_likelihood_trace, type = "l", col = "red",
     main = "Convergence Comparison",
     xlab = "Iteration", ylab = "Log-Likelihood",
     ylim = range(c(result_random$log_likelihood_trace,
                    result_cluster$log_likelihood_trace,
                    result_more_clusters$log_likelihood_trace)))
lines(result_cluster$log_likelihood_trace, col = "blue")
lines(result_more_clusters$log_likelihood_trace, col = "green")
legend("bottomright", 
       legend = c("Random", "Cluster", "More Clusters"),
       col = c("red", "blue", "green"), lty = 1)
