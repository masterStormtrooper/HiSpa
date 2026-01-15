// =========================================================================
// hispa_main.cpp - Command-line interface for HiSpa chromosome analysis using CLI11
// =========================================================================
#include <iostream>
#include <string>
#include <filesystem>
#include <chrono>
#include "CLI11.hpp"
#include "chromosome.h"

int main(int argc, char* argv[]) {
    // --- CLI Setup ---
    CLI::App app{"HiSpa - Hi-C Spatial Analysis Tool"};
    app.set_version_flag("-v,--version", "1.0.0");
    
    // Required arguments
    std::string input_file;
    std::string output_dir;
    
    // Optional parameters with defaults
    std::string prior_positions_file = "";
    std::string prior_contacts_file = "";
    std::vector<int> prior_distances = {1};  // Default: only adjacent distances
    bool use_prior_positions = false;
    bool sample_from_prior = false;
    int mcmc_iterations = 6000;
    int mcmc_burn_in = 0;
    double mcmc_initial_sd = 0.1;
    double mcmc_sd_floor = 0.0001;
    double mcmc_sd_ceil = 0.3;
    int cluster_distance_threshold = 2;
    double cluster_quantile_threshold = -1;
    bool use_cluster_initialization = false;
    int cluster_init_iterations = 1000;
    double cluster_initial_sd = 0.1;
    bool save_samples = false;
    int sample_interval = 50;
    bool verbose = false;
    int num_clusters = 0;
    bool use_quantile_threshold = false;
    bool skip_zero_contact_loci = true;  // Default: skip loci with no contacts
    bool use_convoluted_sampling = false;  // Use convoluted matrices for MCMC sampling
    int convolution_half_k = 3;  // Half window size for convolution (default: 3)
    
    // Add command-line options
    app.add_option("-i,--input", input_file, "Input contact matrix file path")
        ->required()
        ->check(CLI::ExistingFile);
    
    app.add_option("-o,--output", output_dir, "Output directory path")
        ->required();
    
    app.add_option("--prior-positions", prior_positions_file, 
        "Position matrix file to fit gamma prior for adjacent loci distances")
        ->check(CLI::ExistingFile);
    
    app.add_option("--prior-distances", prior_distances, 
        "Genomic distances for fitting priors (e.g., '1,2,5,10' for adjacent, skip-1, skip-4, skip-9)")
        ->delimiter(',');
    
    app.add_option("--prior-contacts", prior_contacts_file, 
        "Contact matrix file from prior data (loaded before clustering)")
        ->check(CLI::ExistingFile);
    
    app.add_flag("--use-prior-positions", use_prior_positions, 
        "Use prior positions as initial structure (skips cluster initialization)");
    
    app.add_flag("--sample-from-prior", sample_from_prior, 
        "Sample new locus positions from prior positions (requires --prior-positions)");
    
    app.add_option("--mcmc-iterations", mcmc_iterations, 
        "Number of MCMC iterations for main algorithm (default: 6000)")
        ->check(CLI::PositiveNumber);

    app.add_option("--num-clusters", num_clusters, 
        "Number of clusters for clustering (default: sqrt(n))")
        ->check(CLI::PositiveNumber);
    
    app.add_option("--mcmc-burn-in", mcmc_burn_in, 
        "Number of burn-in iterations to discard (default: 0)")
        ->check(CLI::NonNegativeNumber);
    
    app.add_option("--mcmc-initial-sd", mcmc_initial_sd, 
        "Initial standard deviation for MCMC proposals (default: 0.3)")
        ->check(CLI::PositiveNumber);
    
    app.add_option("--mcmc-sd-floor", mcmc_sd_floor, 
        "Minimum allowed standard deviation (default: 0.0001)")
        ->check(CLI::PositiveNumber);
    
    app.add_option("--mcmc-sd-ceil", mcmc_sd_ceil, 
        "Maximum allowed standard deviation (default: 2.0)")
        ->check(CLI::PositiveNumber);
    
    app.add_option("--cluster-distance-threshold", cluster_distance_threshold, 
        "Distance threshold for cluster relationships (default: 2)")
        ->check(CLI::PositiveNumber);
    
    app.add_option("--cluster-quantile-threshold", cluster_quantile_threshold, 
        "Quantile threshold for cluster relationships")
        ->check(CLI::PositiveNumber);
    
    app.add_flag("--use-cluster-init", use_cluster_initialization, 
        "Use cluster-based initialization instead of random initialization");
    
    app.add_option("--cluster-init-iterations", cluster_init_iterations, 
        "Number of iterations for cluster initialization MCMC (default: 1000)")
        ->check(CLI::PositiveNumber);
    
    app.add_option("--cluster-initial-sd", cluster_initial_sd, 
        "Initial standard deviation for cluster initialization MCMC (default: 0.1)")
        ->check(CLI::PositiveNumber);
    
    app.add_flag("--save-samples", save_samples, 
        "Save MCMC samples to file during run");
    
    app.add_option("--sample-interval", sample_interval, 
        "Save every kth sample after burn-in (default: 50)")
        ->check(CLI::PositiveNumber);
    
    app.add_flag("--skipzerocontact{true},--no-skipzerocontact{false}", skip_zero_contact_loci,
        "Skip loci with zero contacts during MCMC sampling (default: true)");
    
    app.add_flag("--use-convoluted", use_convoluted_sampling,
        "Use convoluted contact/distance matrices for MCMC sampling");
    
    app.add_option("--convolution-half-k", convolution_half_k,
        "Half window size for convolution (default: 3)")
        ->check(CLI::PositiveNumber);
    
    app.add_flag("--verbose", verbose, "Enable verbose output");
    
    // Parse command line
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }
    
    // --- Validation and Setup ---
    if (verbose) {
        std::cout << "HiSpa - Hi-C Spatial Analysis Tool v1.0.0\n";
        std::cout << "==========================================\n";
        std::cout << "Input file: " << input_file << "\n";
        std::cout << "Output directory: " << output_dir << "\n";
        std::cout << "MCMC iterations: " << mcmc_iterations << "\n";
        std::cout << "MCMC burn-in: " << mcmc_burn_in << "\n";
        std::cout << "MCMC initial SD: " << mcmc_initial_sd << "\n";
        std::cout << "MCMC SD floor: " << mcmc_sd_floor << "\n";
        std::cout << "MCMC SD ceiling: " << mcmc_sd_ceil << "\n";
        std::cout << "Use cluster initialization: " << (use_cluster_initialization ? "Yes" : "No") << "\n";
        if (use_cluster_initialization) {
            std::cout << "Cluster init iterations: " << cluster_init_iterations << "\n";
            std::cout << "Cluster initial SD: " << cluster_initial_sd << "\n";
            std::cout << "Cluster distance threshold: " << cluster_distance_threshold << "\n";
        }
        std::cout << "Save samples: " << (save_samples ? "Yes" : "No") << "\n";
        if (save_samples) {
            std::cout << "Sample interval: " << sample_interval << "\n";
        }
        std::cout << "Use convoluted sampling: " << (use_convoluted_sampling ? "Yes" : "No") << "\n";
        if (use_convoluted_sampling) {
            std::cout << "Convolution half window size: " << convolution_half_k << "\n";
        }
        std::cout << "==========================================\n";
    }
    
    // Create output directory if it doesn't exist
    try {
        if (!std::filesystem::exists(output_dir)) {
            std::filesystem::create_directories(output_dir);
            if (verbose) {
                std::cout << "Created output directory: " << output_dir << "\n";
            }
        }
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error: Could not create output directory '" << output_dir << "': " << e.what() << std::endl;
        return 1;
    }
    
    if (!verbose) {
        std::cout << "HiSpa - Hi-C Spatial Analysis\n";
        std::cout << "==============================\n";
        std::cout << "Input: " << input_file << "\n";
        std::cout << "Output: " << output_dir << "\n\n";
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        // Initialize Chromosome object
        Chromosome my_chromosome(output_dir);
        
        // Load contact matrix data
        if (verbose) {
            std::cout << "Loading contact matrix from: " << input_file << "\n";
        }
        if (!my_chromosome.load_data_from_file(input_file)) {
            std::cerr << "Error: Could not load contact matrix from '" << input_file << "'" << std::endl;
            return 1;
        }
        if (verbose) {
            std::cout << "Successfully loaded contact matrix data.\n";
        }
        
        // Set MCMC options
        my_chromosome.set_skip_zero_contact_loci(skip_zero_contact_loci);
        my_chromosome.set_sample_from_prior(sample_from_prior);
        
        // Compute convoluted matrices if requested
        if (use_convoluted_sampling) {
            if (verbose) {
                std::cout << "Computing convoluted contact matrix (half_k=" << convolution_half_k << ")...\n";
            }
            my_chromosome.compute_convoluted_contact_matrix(convolution_half_k);
            std::cout << "Convoluted contact matrix computed and stored.\n";
        }
        
        // Fit gamma prior if position file is provided
        if (!prior_positions_file.empty()) {
            if (verbose) {
                std::cout << "Fitting gamma priors from position file: " << prior_positions_file << "\n";
                std::cout << "Prior distances: ";
                for (size_t i = 0; i < prior_distances.size(); ++i) {
                    std::cout << prior_distances[i];
                    if (i < prior_distances.size() - 1) std::cout << ", ";
                }
                std::cout << "\n";
            }
            
            if (my_chromosome.fit_distance_priors_from_file(prior_positions_file, prior_distances)) {
                const auto& distance_priors = my_chromosome.get_distance_priors();
                std::cout << "Distance priors fitted for " << distance_priors.active_distances.size() << " distance(s):\n";
                for (int dist : distance_priors.active_distances) {
                    const auto& prior = distance_priors.get_prior(dist);
                    std::cout << "  Distance " << dist << ": shape = " << prior.shape << ", rate = " << prior.rate 
                              << " (mean = " << prior.mean() << ", sd = " << std::sqrt(prior.variance()) << ")\n";
                }
                
                // Validate prior positions if they will be used for initialization
                if (use_prior_positions) {
                    const arma::mat& prior_positions = my_chromosome.get_prior_position_matrix();
                    if (prior_positions.n_rows != my_chromosome.get_contact_matrix().n_rows) {
                        std::cerr << "Error: Prior position matrix has " << prior_positions.n_rows 
                                  << " rows, but contact matrix has " << my_chromosome.get_contact_matrix().n_rows 
                                  << " rows\n";
                        return 1;
                    }
                    if (prior_positions.n_cols != 3) {
                        std::cerr << "Error: Prior position matrix has " << prior_positions.n_cols 
                                  << " columns, but expected 3 (x, y, z)\n";
                        return 1;
                    }
                    if (verbose) {
                        std::cout << "Prior positions validated for initialization (" << prior_positions.n_rows << " x " 
                                  << prior_positions.n_cols << ")\n";
                    }
                }
            } else {
                std::cerr << "Warning: Failed to fit distance priors from file\n";
                if (use_prior_positions) {
                    std::cerr << "Error: Cannot use prior positions if distance prior fitting failed\n";
                    return 1;
                }
            }
        }
        
        // Validate use_prior_positions flag
        if (use_prior_positions && prior_positions_file.empty()) {
            std::cerr << "Error: --use-prior-positions requires --prior-positions to be specified\n";
            return 1;
        }
        
        // Validate sample_from_prior flag
        if (sample_from_prior && prior_positions_file.empty()) {
            std::cerr << "Error: --sample-from-prior requires --prior-positions to be specified\n";
            return 1;
        }
        
        // Load prior contact matrix if provided
        if (!prior_contacts_file.empty()) {
            if (verbose) {
                std::cout << "Loading prior contact matrix from: " << prior_contacts_file << "\n";
            }
            if (!my_chromosome.load_prior_contact_matrix_from_file(prior_contacts_file)) {
                std::cerr << "Error: Failed to load prior contact matrix from file\n";
                return 1;
            }
        }
        
        // --- 1. Pre-processing steps ---
        if (verbose) {
            std::cout << "\n=== Pre-processing ===\n";
            std::cout << "Assigning clusters...\n";
        }
        if (num_clusters > 0) {
            my_chromosome.assign_clusters(num_clusters);
        } else {
            my_chromosome.assign_clusters();
        }
        
        
        
        if (verbose) {
            std::cout << "Building cluster relationships (distance threshold: " << cluster_distance_threshold << ")...\n";
        }
        if (cluster_quantile_threshold > 0) {
            my_chromosome.build_cluster_relationships(cluster_quantile_threshold);
        } else {
            my_chromosome.build_cluster_relationships_by_distance(cluster_distance_threshold);
        }
        
        // --- 2. Initialize structure ---
        if (verbose) {
            std::cout << "\n=== Structure Initialization ===\n";
        }
        
        if (use_prior_positions) {
            if (verbose) {
                std::cout << "Using prior positions as initial structure (skipping cluster initialization)...\n";
            }
            // Skip cluster initialization and use prior positions directly
            // Note: We still need to call initialize_positions() to allocate the position matrix,
            // then we'll overwrite it with prior positions after assemble_global_structure
            my_chromosome.initialize_positions();
        } else if (use_cluster_initialization) {
            if (verbose) {
                std::cout << "Initializing structure from individual cluster MCMC...\n";
                std::cout << "Cluster MCMC iterations: " << cluster_init_iterations << "\n";
                std::cout << "Cluster initial SD: " << cluster_initial_sd << "\n";
            }
            my_chromosome.initialize_structure_from_clusters(cluster_init_iterations, 0, cluster_initial_sd);
        } else {
            if (verbose) {
                std::cout << "Using random position initialization...\n";
            }
            my_chromosome.initialize_positions();
        }
        
        // --- 3. Assemble the global structure ---
        if (verbose) {
            std::cout << "Assembling global structure...\n";
        }
        my_chromosome.assemble_global_structure();
        
        // Override with prior positions if requested
        if (use_prior_positions) {
            if (verbose) {
                std::cout << "Overriding assembled structure with prior positions...\n";
            }
            const arma::mat& prior_positions = my_chromosome.get_prior_position_matrix();
            my_chromosome.set_position_matrix(prior_positions);
            std::cout << "Initial structure set from prior positions.\n";
        }
        
        // Save initial positions
        arma::mat initial_positions = my_chromosome.get_position_matrix();
        std::string initial_pos_file = output_dir + "/initial_positions.txt";
        initial_positions.save(initial_pos_file, arma::raw_ascii);
        if (verbose) {
            std::cout << "Initial positions saved to: " << initial_pos_file << "\n";
        }
        
        // --- 4. Run the main MCMC on the assembled structure ---
        if (verbose) {
            std::cout << "\n=== Main MCMC Algorithm ===\n";
            std::cout << "Running MCMC with " << mcmc_iterations << " iterations...\n";
            std::cout << "Burn-in: " << mcmc_burn_in << "\n";
            std::cout << "Initial SD: " << mcmc_initial_sd << "\n";
            std::cout << "SD floor: " << mcmc_sd_floor << "\n";
            std::cout << "SD ceiling: " << mcmc_sd_ceil << "\n";
            std::cout << "Using convoluted sampling: " << (use_convoluted_sampling ? "Yes" : "No") << "\n";
            if (save_samples) {
                std::cout << "Saving samples: every " << sample_interval << " iterations after burn-in\n";
            }
        } else {
            std::cout << "Running MCMC (" << mcmc_iterations << " iterations";
            if (save_samples) {
                std::cout << ", saving samples";
            }
            if (use_convoluted_sampling) {
                std::cout << ", using convoluted matrices";
            }
            std::cout << ")...\n";
        }
        
        if (use_convoluted_sampling) {
            my_chromosome.run_mcmc_convoluted(mcmc_iterations, mcmc_burn_in, mcmc_initial_sd, mcmc_sd_floor, mcmc_sd_ceil, save_samples, sample_interval, convolution_half_k);
        } else {
            my_chromosome.run_mcmc(mcmc_iterations, mcmc_burn_in, mcmc_initial_sd, mcmc_sd_floor, mcmc_sd_ceil, save_samples, sample_interval);
        }
        
        // --- 5. Save final results ---
        if (verbose) {
            std::cout << "\n=== Saving Results ===\n";
        }
        std::string chr_name = my_chromosome.get_name();
        
        // Save final positions
        arma::mat best_positions = my_chromosome.get_best_position_matrix();
        std::string final_pos_file = chr_name + "/final_positions.txt";
        best_positions.save(final_pos_file, arma::raw_ascii);
        std::cout << "Final 3D coordinates saved to: " << final_pos_file << std::endl;
        
        // Save log-likelihood trace
        arma::vec ll_trace = my_chromosome.get_mcmc_trace_log_likelihood();
        std::string ll_trace_file = chr_name + "/log_likelihood_trace.txt";
        ll_trace.save(ll_trace_file, arma::raw_ascii);
        if (verbose) {
            std::cout << "MCMC log-likelihood trace saved to: " << ll_trace_file << std::endl;
        }
        
        // Save block timings
        const auto& time_trace_vec = my_chromosome.get_mcmc_trace_block_durations();
        arma::vec time_trace_arma(time_trace_vec);
        std::string timing_file = chr_name + "/block_timings.txt";
        time_trace_arma.save(timing_file, arma::raw_ascii);
        if (verbose) {
            std::cout << "MCMC block timings saved to: " << timing_file << std::endl;
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::cout << "\n=== Analysis Complete ===\n";
        std::cout << "Total execution time: " << duration.count() << " seconds\n";
        std::cout << "Results saved in: " << output_dir << "\n";
        if (verbose) {
            std::cout << "Program finished successfully.\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error during algorithm execution: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
