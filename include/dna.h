// =========================================================================
// dna.h - Header file for the DNA class
// =========================================================================
#ifndef DNA_H
#define DNA_H

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <armadillo>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include <cmath>        // For std::exp, std::log, std::lgamma, std::isfinite
#include <limits>       // For std::numeric_limits

/*
 * @struct ClusterAdjacency
 * @brief Stores the row and column indices of segments belonging to a cluster, its neighbors,
 * and its non-neighboring "stranger" clusters for DNA microscopy data.
 */
struct ClusterAdjacency {
    arma::uvec self_row_indices;
    arma::uvec self_col_indices;
    arma::uvec neighbor_row_indices;
    arma::uvec neighbor_col_indices;
    arma::uvec stranger_row_indices;
    arma::uvec stranger_col_indices;
};

/*
 * @class DNA
 * @brief Represents a DNA molecule and its associated microscopy contact data.
 * The contact matrix is asymmetric and should be read as-is.
 */
class DNA {
private:
    std::string dna_name;
    arma::mat contact_matrix; // Asymmetric matrix from DNA microscopy
    arma::uvec row_cluster_labels;
    arma::uvec col_cluster_labels;
    arma::mat geometric_mean_matrix; // n_clusters x n_clusters matrix (used as backbone contact matrix)
    arma::mat backbone_contact_matrix; // Symmetric backbone contact matrix derived from geometric means
    std::vector<ClusterAdjacency> cluster_adjacencies;
    arma::mat row_position_matrix; // 2D positions for rows (n_rows x 2)
    arma::mat col_position_matrix; // 2D positions for columns (n_cols x 2)
    arma::mat pairwise_distance_matrix; // Distance matrix between rows and columns
    double beta0;
    double beta1;
    arma::mat best_row_position_matrix;
    arma::mat best_col_position_matrix;
    double max_log_likelihood;
    
    // Initialization results
    std::vector<arma::mat> initial_cluster_row_structures;
    std::vector<arma::mat> initial_cluster_col_structures;
    arma::mat backbone_row_structure;
    arma::mat backbone_col_structure;
    std::vector<double> initial_beta0s;
    std::vector<double> initial_beta1s;
    double backbone_beta0;
    double backbone_beta1;
    double dispersion_k = 0.2;  // Dispersion parameter for Negative Binomial model
    // Logger
    std::shared_ptr<spdlog::logger> logger;

    /*
     * @brief Log-likelihood calculation for DNA microscopy data using d^2 instead of log(d).
     * Lambda = exp(b0 + b1 * d^2)
     * Log-likelihood = sum(c * log(lambda) - lambda)
     */
    template<typename T1, typename T2>
    double calculate_log_likelihood(
        const arma::Base<double, T1>& distances_expr, 
        const arma::Base<double, T2>& contacts_expr, 
        double b0, double b1) const 
    {
        const T1& distances = distances_expr.get_ref();
        const T2& contacts = contacts_expr.get_ref();

        double total_log_likelihood = 0.0;
        
        const double* dist_mem = distances.memptr();
        const double* cont_mem = contacts.memptr();
        const arma::uword n_elem = distances.n_elem;

        for (arma::uword i = 0; i < n_elem; ++i) {
            const double d = dist_mem[i];
            const double c = cont_mem[i];
            
            if (d <= 0) {
                continue;
            }
            // if (c == 0) {
            //     continue;
            // }

            // const double c = cont_mem[i];
            const double d_squared = d * d;
            const double log_lam = b0 + b1 * d_squared;
            const double term = c * log_lam - std::exp(log_lam);
            
            if (std::isfinite(term)) {
                total_log_likelihood += term;
            }
        }
        return total_log_likelihood;
    }


    // template<typename T1, typename T2>
    // double calculate_log_likelihood(
    //     const arma::Base<double, T1>& distances_expr, 
    //     const arma::Base<double, T2>& contacts_expr, 
    //     double b0, double b1) const 
    // {
    //     const T1& distances = distances_expr.get_ref();
    //     const T2& contacts = contacts_expr.get_ref();

    //     // The dispersion parameter 'k' must be positive.
    //     // If not, return -infinity as this is invalid for maximization.
    //     if (dispersion_k <= 0) {
    //         return -std::numeric_limits<double>::infinity();
    //     }

    //     double total_log_likelihood = 0.0;
        
    //     const double* dist_mem = distances.memptr();
    //     const double* cont_mem = contacts.memptr();
    //     const arma::uword n_elem = distances.n_elem;

    //     // Pre-calculate log(k) since it's constant in the loop
    //     const double log_k = std::log(dispersion_k);
    //     const double lgamma_k = std::lgamma(dispersion_k);

    //     for (arma::uword i = 0; i < n_elem; ++i) {
    //         const double d = dist_mem[i];
    //         const double c = cont_mem[i];
            
    //         // Mean (µ) calculation, same as your Poisson's lambda
    //         const double d_squared = d * d;
    //         const double log_mu = b0 + b1 * d_squared;
    //         const double mu = std::exp(log_mu);

    //         // Negative Binomial log-likelihood for a single observation (c, mu, k)
    //         // We drop std::lgamma(c + 1) [which is log(c!)] as it's
    //         // a constant with respect to the parameters (b0, b1, k)
    //         // and doesn't affect optimization.
    //         const double term = std::lgamma(c + dispersion_k) - lgamma_k + 
    //                             dispersion_k * log_k + 
    //                             c * log_mu - 
    //                             (c + dispersion_k) * std::log(mu + dispersion_k);
            
    //         if (std::isfinite(term)) {
    //             total_log_likelihood += term;
    //         }
    //     }
    //     return total_log_likelihood;
    // }

    arma::mat calculate_pairwise_distances(const arma::mat& pos1, const arma::mat& pos2) const;

public:
    DNA(const std::string& name);
    ~DNA();

    // Data loading functions
    bool load_data_from_file(const std::string& filename);
    bool load_row_cluster_labels(const std::string& filename);
    bool load_col_cluster_labels(const std::string& filename);

    // Cluster relationship building
    void build_cluster_relationship(int top_k = 3);

    // Position initialization and MCMC sampling
    void initialize_positions();
    // Initialize positions from provided matrices (rows x 2) and set parameters/log-likelihood similarly
    void initialize_positions(const arma::mat& row_pos, const arma::mat& col_pos);
    void run_mcmc(int iterations, int burn_in, double initial_sd = 0.1, double sd_floor = 0.001, double sd_ceiling = 0.3, bool save_samples = false, int sample_interval = 5);

    double calculate_log_likelihood(const arma::mat& distances, const arma::mat& contacts) const {
        return calculate_log_likelihood(distances, contacts, this->beta0, this->beta1);
    }
    arma::mat rotate_positions(const arma::mat& positions, const arma::mat& rotation) const;
    arma::mat get_clockwise_rotation(double angle_degrees);
    // arma::mat optimize_via_rotation(const arma::mat& positions, int num_angles = 36);
    // Structure initialization and assembly
    void initialize_structure_from_clusters(int sub_iterations, int sub_burn_in, double initial_sd = 0.1);
    void assemble_global_structure();
    void save_assembled_structures(const std::string& output_dir = "");
    
    // Getters
    const std::string& get_name() const;
    const arma::mat& get_contact_matrix() const;
    const arma::uvec& get_row_cluster_labels() const;
    const arma::uvec& get_col_cluster_labels() const;
    const arma::mat& get_geometric_mean_matrix() const;
    const arma::mat& get_backbone_contact_matrix() const;
    const std::vector<ClusterAdjacency>& get_cluster_adjacencies() const;
    const arma::mat& get_row_position_matrix() const;
    const arma::mat& get_col_position_matrix() const;
    const arma::mat& get_best_row_position_matrix() const;
    const arma::mat& get_best_col_position_matrix() const;
};

#endif // DNA_H
