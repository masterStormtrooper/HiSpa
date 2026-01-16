#' Run HiSpa MCMC Analysis
#'
#' @description
#' Performs hierarchical Bayesian inference of 3D chromatin structure from
#' Hi-C contact matrix using MCMC sampling. This function follows the exact
#' workflow from the HiSpa C++ implementation.
#'
#' @param input_file Character string, path to contact matrix file.
#' @param output_dir Character string specifying the output directory path.
#' @param mcmc_iterations Integer, number of MCMC iterations (default: 6000).
#' @param num_clusters Integer, number of clusters for hierarchical analysis. 
#'   If 0 (default), automatically determined as sqrt(n).
#' @param mcmc_burn_in Integer, number of burn-in iterations to discard (default: 0).
#' @param mcmc_initial_sd Numeric, initial standard deviation for MCMC proposals (default: 0.1).
#' @param mcmc_sd_floor Numeric, minimum allowed standard deviation (default: 0.0001).
#' @param mcmc_sd_ceil Numeric, maximum allowed standard deviation (default: 0.3).
#' @param use_cluster_init Logical, use cluster-based initialization instead of 
#'   random initialization (default: FALSE).
#' @param cluster_init_iterations Integer, number of iterations for cluster 
#'   initialization MCMC (default: 1000).
#' @param cluster_initial_sd Numeric, initial standard deviation for cluster 
#'   initialization MCMC (default: 0.1).
#' @param save_samples Logical, whether to save MCMC trace samples (default: FALSE).
#' @param sample_interval Integer, save samples every n iterations (default: 50).
#' @param verbose Logical, enable verbose output (default: TRUE).
#'
#' @return Invisibly returns the output directory path. All analysis results
#'   are saved as text files in the output directory (see Details).
#'
#' @details
#' This function implements the complete HiSpa workflow:
#' \enumerate{
#'   \item Load contact matrix from file
#'   \item Assign loci to clusters (k-means)
#'   \item Build cluster relationships
#'   \item Initialize structure (random or cluster-based)
#'   \item Assemble global structure
#'   \item Run main MCMC algorithm
#'   \item Save results to output directory
#' }
#' 
#' All results are automatically saved to the output directory:
#' \itemize{
#'   \item \strong{final_positions.txt} - Final inferred 3D coordinates (n x 3 matrix)
#'   \item \strong{initial_positions.txt} - Initial positions before MCMC (n x 3 matrix)
#'   \item \strong{log_likelihood_trace.txt} - MCMC log-likelihood values (convergence diagnostic)
#'   \item \strong{block_timings.txt} - Computation time for each MCMC block
#'   \item \strong{mcmc_log.txt} - Detailed analysis log with parameter values and settings
#' }
#' 
#' Read results using standard R functions:
#' \code{final_pos <- read.table("output_dir/final_positions.txt")}
#' \code{ll_trace <- scan("output_dir/log_likelihood_trace.txt")}
#'
#' @examples
#' \dontrun{
#' # Save contact matrix to file
#' write.table(contact_mat, "contact_matrix.txt", 
#'             row.names = FALSE, col.names = FALSE)
#' 
#' # Run analysis - all results saved to output directory
#' hispa_analyze(
#'   input_file = "contact_matrix.txt",
#'   output_dir = "output",
#'   mcmc_iterations = 6000,
#'   mcmc_burn_in = 1000,
#'   num_clusters = 0,  # auto-detect
#'   use_cluster_init = TRUE,
#'   cluster_init_iterations = 1000,
#'   verbose = TRUE
#' )
#' 
#' # Read results from output directory
#' final_pos <- as.matrix(read.table("output/final_positions.txt"))
#' initial_pos <- as.matrix(read.table("output/initial_positions.txt"))
#' ll_trace <- scan("output/log_likelihood_trace.txt")
#' 
#' # Plot convergence
#' plot(ll_trace, type = "l",
#'      xlab = "Iteration", ylab = "Log-Likelihood")
#' 
#' # 3D visualization
#' library(rgl)
#' plot3d(final_pos, col = rainbow(nrow(final_pos)), size = 5)
#' }
#'
#' @export
hispa_analyze <- function(
    input_file,
    output_dir,
    mcmc_iterations = 6000L,
    num_clusters = 0L,
    mcmc_burn_in = 0L,
    mcmc_initial_sd = 0.1,
    mcmc_sd_floor = 0.0001,
    mcmc_sd_ceil = 0.3,
    use_cluster_init = FALSE,
    cluster_init_iterations = 1000L,
    cluster_initial_sd = 0.1,
    save_samples = FALSE,
    sample_interval = 50L,
    verbose = TRUE
) {
  # Input validation
  if (!is.character(input_file) || length(input_file) != 1) {
    stop("input_file must be a single character string")
  }
  
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }
  
  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be a single character string")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Call C++ function (returns output directory path)
  output_path <- .hispa_analyze_cpp(
                  input_file,
                  output_dir,
                  as.integer(mcmc_iterations),
                  as.integer(num_clusters),
                  as.integer(mcmc_burn_in),
                  as.numeric(mcmc_initial_sd),
                  as.numeric(mcmc_sd_floor),
                  as.numeric(mcmc_sd_ceil),
                  as.logical(use_cluster_init),
                  as.integer(cluster_init_iterations),
                  as.numeric(cluster_initial_sd),
                  as.logical(save_samples),
                  as.integer(sample_interval),
                  as.logical(verbose))
  
  # Return output directory invisibly
  invisible(output_path)
}
