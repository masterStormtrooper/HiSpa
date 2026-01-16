# Example: Basic HiSpa Analysis
# This example demonstrates how to run a basic HiSpa analysis

library(HiSpaR)

# Save to file (required by HiSpaR)
input_file <- "/Users/lyc/Proj/HiSpa/R-package/data/su1_contact_mat.txt"
# Run HiSpa analysis with default settings (random initialization)
hispa_analyze(
  input_file = input_file,
  output_dir = "example_output",
  mcmc_iterations = 2000,
  mcmc_burn_in = 10,
#   num_clusters = 0,  # Auto-detect (sqrt(n))
  use_cluster_init = TRUE,  # Use random initialization
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
     main = "MCMC Convergence")
abline(v = 1000, col = "red", lty = 2)  # Burn-in line
legend("bottomright", legend = "Burn-in", col = "red", lty = 2)

# Basic 3D visualization (requires rgl package)
if (requireNamespace("rgl", quietly = TRUE)) {
  rgl::plot3d(final_positions, 
              col = rainbow(nrow(final_positions)),
              size = 5,
              xlab = "X", ylab = "Y", zlab = "Z",
              main = "Inferred 3D Chromatin Structure")
  
  # Connect sequential loci
  for (i in 1:(nrow(final_positions)-1)) {
    rgl::lines3d(final_positions[c(i, i+1), ], col = "gray", lwd = 2)
  }
}

# Compare initial vs final positions
par(mfrow = c(1, 2))
plot(initial_positions[, 1], initial_positions[, 2], 
     col = rainbow(n), pch = 16,
     xlab = "X", ylab = "Y", main = "Initial Positions (XY plane)")
plot(final_positions[, 1], final_positions[, 2], 
     col = rainbow(n), pch = 16,
     xlab = "X", ylab = "Y", main = "Final Positions (XY plane)")

# Clean up temp file
unlink(input_file)

cat("\nOutput files saved in:", result$output_dir, "\n")
cat("  - final_positions.txt\n")
cat("  - initial_positions.txt\n")
cat("  - log_likelihood_trace.txt\n")
cat("  - block_timings.txt\n")
