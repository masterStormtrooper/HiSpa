// =========================================================================
// dna.cpp - Implementation file for the DNA class
// =========================================================================
#include "dna.h"

DNA::DNA(const std::string& name) : dna_name(name), beta0(0.0), beta1(0.0), max_log_likelihood(-arma::datum::inf) {
    std::cout << "DNA object '" << dna_name << "' created." << std::endl;
    try {
        std::filesystem::create_directory(dna_name);
        std::string log_path = dna_name + "/mcmc_log.txt";
        std::string logger_name = "logger_" + name;
        std::replace(logger_name.begin(), logger_name.end(), '/', '_');
        logger = spdlog::basic_logger_mt(logger_name, log_path, true); 
    } catch (const spdlog::spdlog_ex &ex) {
        std::cerr << "Log initialization failed: " << ex.what() << std::endl;
    }
}

DNA::~DNA() {
    std::cout << "DNA object '" << dna_name << "' destroyed." << std::endl;
}

bool DNA::load_data_from_file(const std::string& filename) {
    std::cout << "Loading DNA microscopy data from '" << filename << "'..." << std::endl;
    bool success = contact_matrix.load(filename, arma::raw_ascii);

    if (success) {
        std::cout << "Successfully loaded asymmetric contact matrix with size: "
                  << contact_matrix.n_rows << " x " << contact_matrix.n_cols
                  << std::endl;
        logger->info("Loaded DNA microscopy data: {} rows x {} columns", 
                    contact_matrix.n_rows, contact_matrix.n_cols);
    } else {
        std::cerr << "Error: Failed to load DNA microscopy data from '" << filename << "'." << std::endl;
    }
    return success;
}

bool DNA::load_row_cluster_labels(const std::string& filename) {
    std::cout << "Loading row cluster labels from '" << filename << "'..." << std::endl;
    bool success = row_cluster_labels.load(filename, arma::raw_ascii);

    if (success) {
        std::cout << "Successfully loaded " << row_cluster_labels.n_elem << " row cluster labels" << std::endl;
        
        // Validate that the number of labels matches the number of rows in contact matrix
        if (!contact_matrix.is_empty() && row_cluster_labels.n_elem != contact_matrix.n_rows) {
            std::cerr << "Warning: Number of row cluster labels (" << row_cluster_labels.n_elem 
                      << ") does not match contact matrix rows (" << contact_matrix.n_rows << ")" << std::endl;
        }
        
        logger->info("Loaded {} row cluster labels", row_cluster_labels.n_elem);
    } else {
        std::cerr << "Error: Failed to load row cluster labels from '" << filename << "'." << std::endl;
    }
    return success;
}

bool DNA::load_col_cluster_labels(const std::string& filename) {
    std::cout << "Loading column cluster labels from '" << filename << "'..." << std::endl;
    bool success = col_cluster_labels.load(filename, arma::raw_ascii);

    if (success) {
        std::cout << "Successfully loaded " << col_cluster_labels.n_elem << " column cluster labels" << std::endl;
        
        // Validate that the number of labels matches the number of columns in contact matrix
        if (!contact_matrix.is_empty() && col_cluster_labels.n_elem != contact_matrix.n_cols) {
            std::cerr << "Warning: Number of column cluster labels (" << col_cluster_labels.n_elem 
                      << ") does not match contact matrix columns (" << contact_matrix.n_cols << ")" << std::endl;
        }
        
        logger->info("Loaded {} column cluster labels", col_cluster_labels.n_elem);
    } else {
        std::cerr << "Error: Failed to load column cluster labels from '" << filename << "'." << std::endl;
    }
    return success;
}

void DNA::build_cluster_relationship(int top_k) {
    if (row_cluster_labels.is_empty() || col_cluster_labels.is_empty()) {
        std::cerr << "Error: Cannot build relationships. Load row and column cluster labels first." << std::endl;
        return;
    }
    
    if (contact_matrix.is_empty()) {
        std::cerr << "Error: Cannot build relationships. Load contact matrix first." << std::endl;
        return;
    }

    // Get unique cluster labels from both row and column labels
    arma::uvec unique_row_labels = arma::unique(row_cluster_labels);
    arma::uvec unique_col_labels = arma::unique(col_cluster_labels);
    
    // For DNA microscopy, we'll use row clusters as the main clusters
    arma::uword num_row_clusters = unique_row_labels.n_elem;
    arma::uword num_col_clusters = unique_col_labels.n_elem;
    std::cout << "Building relationships for " << num_row_clusters << " row clusters and " << num_col_clusters << " column clusters using modified geometric means..." << std::endl;
    
    // Create indices by cluster for rows
    std::vector<arma::uvec> row_indices_by_cluster(num_row_clusters);
    for (arma::uword i = 0; i < row_cluster_labels.n_elem; ++i) {
        arma::uword cluster_id = row_cluster_labels(i);
        if (row_indices_by_cluster[cluster_id].is_empty()) {
            row_indices_by_cluster[cluster_id] = arma::uvec{i};
        } else {
            row_indices_by_cluster[cluster_id] = arma::join_cols(row_indices_by_cluster[cluster_id], arma::uvec{i});
        }
    }
    
    // Create indices by cluster for columns
    std::vector<arma::uvec> col_indices_by_cluster(num_col_clusters);
    for (arma::uword i = 0; i < col_cluster_labels.n_elem; ++i) {
        arma::uword cluster_id = col_cluster_labels(i);
        if (col_indices_by_cluster[cluster_id].is_empty()) {
            col_indices_by_cluster[cluster_id] = arma::uvec{i};
        } else {
            col_indices_by_cluster[cluster_id] = arma::join_cols(col_indices_by_cluster[cluster_id], arma::uvec{i});
        }
    }
    
    // Calculate modified geometric mean matrix (num_row_clusters x num_row_clusters)
    geometric_mean_matrix.set_size(num_row_clusters, num_row_clusters);
    
    for (arma::uword c1 = 0; c1 < num_row_clusters; ++c1) {
        for (arma::uword c2 = 0; c2 < num_row_clusters; ++c2) {
            // Get the submatrix between cluster c1 (rows) and cluster c2 (columns)
            arma::mat sub_contact_matrix = contact_matrix.submat(
                row_indices_by_cluster[c1], 
                col_indices_by_cluster[c2]
            );
            
            // Calculate modified geometric mean: add 1 to all entries, calculate geometric mean, then subtract 1
            arma::mat modified_matrix = sub_contact_matrix + 1.0; // Add 1 to all entries
            arma::vec flat_modified = arma::vectorise(modified_matrix);
            
            // Calculate geometric mean
            double geometric_mean = 0.0;
            if (flat_modified.n_elem > 0) {
                double sum_log = 0.0;
                arma::uword valid_count = 0;
                for (arma::uword k = 0; k < flat_modified.n_elem; ++k) {
                    if (flat_modified(k) > 0) {
                        sum_log += std::log(flat_modified(k));
                        valid_count++;
                    }
                }
                if (valid_count > 0) {
                    geometric_mean = std::exp(sum_log / valid_count) - 1.0; // Subtract 1 after geometric mean
                }
            }
            
            geometric_mean_matrix(c1, c2) = geometric_mean;
        }
    }
    
    // Save geometric mean matrix as backbone contact matrix for later use
    backbone_contact_matrix = geometric_mean_matrix;
    
    // Build cluster adjacencies based on top-k geometric means
    cluster_adjacencies.resize(num_row_clusters);
    
    for (arma::uword c = 0; c < num_row_clusters; ++c) {
        cluster_adjacencies[c].self_row_indices = row_indices_by_cluster[c];
        cluster_adjacencies[c].self_col_indices = col_indices_by_cluster[c];
        
        // Get the row of geometric means for this cluster
        arma::rowvec geometric_means = geometric_mean_matrix.row(c);
        
        // Create pairs of (geometric_mean, cluster_index) excluding self
        std::vector<std::pair<double, arma::uword>> cluster_means;
        for (arma::uword other_c = 0; other_c < num_row_clusters; ++other_c) {
            if (other_c != c) {
                cluster_means.push_back(std::make_pair(geometric_means(other_c), other_c));
            }
        }
        
        // Sort by geometric mean (descending)
        std::sort(cluster_means.begin(), cluster_means.end(), 
                  [](const std::pair<double, arma::uword>& a, const std::pair<double, arma::uword>& b) {
                      return a.first > b.first;
                  });
        
        // Select top-k clusters as neighbors
        std::vector<arma::uword> neighbor_row_indices, neighbor_col_indices;
        std::vector<arma::uword> stranger_row_indices, stranger_col_indices;
        
        for (arma::uword i = 0; i < cluster_means.size(); ++i) {
            arma::uword cluster_id = cluster_means[i].second;
            const arma::uvec& cluster_row_indices = row_indices_by_cluster[cluster_id];
            const arma::uvec& cluster_col_indices = col_indices_by_cluster[cluster_id];
            
            if (i < top_k) {
                // This is a neighbor cluster
                neighbor_row_indices.insert(neighbor_row_indices.end(), 
                                          cluster_row_indices.begin(), 
                                          cluster_row_indices.end());
                neighbor_col_indices.insert(neighbor_col_indices.end(), 
                                          cluster_col_indices.begin(), 
                                          cluster_col_indices.end());
            } else {
                // This is a stranger cluster
                stranger_row_indices.insert(stranger_row_indices.end(), 
                                          cluster_row_indices.begin(), 
                                          cluster_row_indices.end());
                stranger_col_indices.insert(stranger_col_indices.end(), 
                                          cluster_col_indices.begin(), 
                                          cluster_col_indices.end());
            }
        }
        
        cluster_adjacencies[c].neighbor_row_indices = arma::conv_to<arma::uvec>::from(neighbor_row_indices);
        cluster_adjacencies[c].neighbor_col_indices = arma::conv_to<arma::uvec>::from(neighbor_col_indices);
        cluster_adjacencies[c].stranger_row_indices = arma::conv_to<arma::uvec>::from(stranger_row_indices);
        cluster_adjacencies[c].stranger_col_indices = arma::conv_to<arma::uvec>::from(stranger_col_indices);
    }
    
    // Log the relationships
    logger->info("--- Cluster Adjacency Report (Geometric Mean-based, top_k={}) ---", top_k);
    for (arma::uword c = 0; c < num_row_clusters; ++c) {
        const auto& adj = cluster_adjacencies[c];
        std::string neighbors_str = "None";
        if (!adj.neighbor_row_indices.is_empty()) {
            arma::uvec neighbor_labels = arma::unique(row_cluster_labels(adj.neighbor_row_indices));
            std::stringstream ss;
            neighbor_labels.t().print(ss);
            neighbors_str = ss.str();
            neighbors_str.erase(std::remove(neighbors_str.begin(), neighbors_str.end(), '\n'), neighbors_str.end());
        }
        std::string strangers_str = "None";
        if (!adj.stranger_row_indices.is_empty()) {
            arma::uvec stranger_labels = arma::unique(row_cluster_labels(adj.stranger_row_indices));
            std::stringstream ss;
            stranger_labels.t().print(ss);
            strangers_str = ss.str();
            strangers_str.erase(std::remove(strangers_str.begin(), strangers_str.end(), '\n'), strangers_str.end());
        }
        logger->info("Cluster {}: Self rows={}, cols={}; Neighbors [{}], Strangers [{}]", 
                    c, adj.self_row_indices.n_elem, adj.self_col_indices.n_elem, neighbors_str, strangers_str);
    }
    
    std::cout << "Cluster relationships built successfully. Geometric mean matrix size: " 
              << geometric_mean_matrix.n_rows << " x " << geometric_mean_matrix.n_cols << std::endl;
}

// Getters
const std::string& DNA::get_name() const { 
    return dna_name; 
}

const arma::mat& DNA::get_contact_matrix() const { 
    return contact_matrix; 
}

const arma::uvec& DNA::get_row_cluster_labels() const { 
    return row_cluster_labels; 
}

const arma::uvec& DNA::get_col_cluster_labels() const { 
    return col_cluster_labels; 
}

const arma::mat& DNA::get_geometric_mean_matrix() const { 
    return geometric_mean_matrix; 
}

const arma::mat& DNA::get_backbone_contact_matrix() const { 
    return backbone_contact_matrix; 
}

const std::vector<ClusterAdjacency>& DNA::get_cluster_adjacencies() const { 
    return cluster_adjacencies; 
}
arma::mat DNA::rotate_positions(const arma::mat& positions, const arma::mat& rotation) const {
    // Validate input
    if (positions.is_empty()) {
        logger->warn("rotate_positions called with empty positions matrix.");
        return positions;
    }
    if (rotation.is_empty()) {
        logger->warn("rotate_positions called with empty rotation matrix. Returning original positions.");
        return positions;
    }

    arma::uword n = positions.n_rows;
    arma::uword d = positions.n_cols;

    if (rotation.n_rows != d || rotation.n_cols != d) {
        logger->warn("rotate_positions: rotation matrix dimension ({}, {}) does not match positions dimensionality ({}). Returning original positions.", rotation.n_rows, rotation.n_cols, d);
        return positions;
    }

    // Compute centroid (row vector)
    arma::rowvec centroid = arma::mean(positions, 0);

    // Convert to relative positions (n x d), apply rotation, and add centroid back
    arma::mat relative = positions.each_row() - centroid;
    arma::mat rotated_relative = relative * rotation.t();
    arma::mat result = rotated_relative;
    result.each_row() += centroid;

    return result;
}


arma::mat DNA::get_clockwise_rotation(double degrees) {
    // 1. Convert the input angle from degrees to radians
    // C++ math functions (sin, cos) expect radians.
    double radians = degrees * arma::datum::pi / 180.0;

    // 2. Calculate the cosine and sine of the angle
    double c = std::cos(radians);
    double s = std::sin(radians);
    arma::mat R = { { c,  s },
                    { -s, c } };

    return R;
}

// arma::mat DNA::optimize_via_rotation(const arma::mat& positions, int num_angles) {
//     if (positions.is_empty()) {
//         logger->warn("optimize_via_rotation called with empty positions matrix.");
//         return positions;
//     }

//     double best_loglikelihood = -arma::datum::inf;
//     arma::mat best_positions = positions;

//     for (int i = 0; i < num_angles; ++i) {
//         double angle = (360.0 / num_angles) * i;
//         arma::mat rotation = get_clockwise_rotation(angle);
//         arma::mat rotated_positions = rotate_positions(positions, rotation);
//         double variance = arma::var(arma::vectorise(rotated_positions));

//         if (variance < best_variance) {
//             best_variance = variance;
//             best_positions = rotated_positions;
//         }
//     }

//     logger->info("optimize_via_rotation: Best variance after optimization: {:.4f}", best_variance);
//     return best_positions;
// }


arma::mat DNA::calculate_pairwise_distances(const arma::mat& pos1, const arma::mat& pos2) const {
    // Efficient pairwise distance calculation using matrix operations
    // pos1: n1 x 2 matrix (2D positions)
    // pos2: n2 x 2 matrix (2D positions)
    // Returns: n1 x n2 distance matrix
    
    arma::mat p1_2 = arma::sum(arma::square(pos1), 1); // Column vector of squared norms
    arma::mat p2_2 = arma::sum(arma::square(pos2), 1); // Column vector of squared norms
    arma::mat p1p2 = -2 * pos1 * pos2.t();              // Cross product term
    p1p2.each_col() += p1_2;                             // Add p1_2 to each column
    p1p2.each_row() += p2_2.t();                         // Add p2_2 to each row
    return arma::sqrt(p1p2.clamp(0, p1p2.max()));        // Clamp to avoid negative values due to numerical error
}

void DNA::initialize_positions() {
    if (contact_matrix.is_empty()) {
        std::cerr << "Error: Load contact matrix before initializing positions." << std::endl;
        return;
    }
    
    arma::uword n_rows = contact_matrix.n_rows;
    arma::uword n_cols = contact_matrix.n_cols;
    
    std::cout << "Initializing " << n_rows << " row positions and " << n_cols << " column positions randomly..." << std::endl;
    
    // Initialize row positions in 2D - random Gaussian distribution
    row_position_matrix.set_size(n_rows, 2);
    row_position_matrix = arma::randn<arma::mat>(n_rows, 2); // Standard Gaussian N(0,1)
    
    // Initialize column positions in 2D - random Gaussian distribution
    col_position_matrix.set_size(n_cols, 2);
    col_position_matrix = arma::randn<arma::mat>(n_cols, 2); // Standard Gaussian N(0,1)
    
    std::cout << "Row positions - mean: [" << arma::mean(row_position_matrix.col(0)) << ", " 
              << arma::mean(row_position_matrix.col(1)) << "], std: [" 
              << arma::stddev(row_position_matrix.col(0)) << ", " 
              << arma::stddev(row_position_matrix.col(1)) << "]" << std::endl;
    std::cout << "Col positions - mean: [" << arma::mean(col_position_matrix.col(0)) << ", " 
              << arma::mean(col_position_matrix.col(1)) << "], std: [" 
              << arma::stddev(col_position_matrix.col(0)) << ", " 
              << arma::stddev(col_position_matrix.col(1)) << "]" << std::endl;
    
    // Calculate initial pairwise distance matrix
    std::cout << "Calculating initial pairwise distance matrix..." << std::endl;
    pairwise_distance_matrix = calculate_pairwise_distances(row_position_matrix, col_position_matrix);
    
    // Initialize beta parameters
    beta0 = 4;
    beta1 = -1;
    
    // Initialize best state
    max_log_likelihood = calculate_log_likelihood(pairwise_distance_matrix, contact_matrix);
    best_row_position_matrix = row_position_matrix;
    best_col_position_matrix = col_position_matrix;
    
    std::cout << "Initial log-likelihood: " << max_log_likelihood << std::endl;
    logger->info("Positions initialized randomly: {} rows, {} cols, initial LL: {:.4f}", n_rows, n_cols, max_log_likelihood);
}

// Initialize positions using provided matrices and set up parameters/log-likelihood
void DNA::initialize_positions(const arma::mat& row_pos, const arma::mat& col_pos) {
    if (contact_matrix.is_empty()) {
        std::cerr << "Error: Load contact matrix before initializing positions." << std::endl;
        return;
    }

    arma::uword n_rows = contact_matrix.n_rows;
    arma::uword n_cols = contact_matrix.n_cols;

    if (row_pos.n_rows != n_rows || row_pos.n_cols != 2) {
        std::cerr << "Error: Provided row_pos must be " << n_rows << " x 2." << std::endl;
        return;
    }
    if (col_pos.n_rows != n_cols || col_pos.n_cols != 2) {
        std::cerr << "Error: Provided col_pos must be " << n_cols << " x 2." << std::endl;
        return;
    }

    std::cout << "Initializing positions from provided matrices for " << n_rows << " rows and " << n_cols << " cols..." << std::endl;

    // Assign provided positions
    row_position_matrix = row_pos;
    col_position_matrix = col_pos;

    // Calculate pairwise distances
    pairwise_distance_matrix = calculate_pairwise_distances(row_position_matrix, col_position_matrix);

    // Initialize beta parameters to same defaults as random init
    beta0 = 4.0;
    beta1 = -1.0;

    // Initialize best state and log-likelihood
    max_log_likelihood = calculate_log_likelihood(pairwise_distance_matrix, contact_matrix);
    best_row_position_matrix = row_position_matrix;
    best_col_position_matrix = col_position_matrix;

    std::cout << "Initialized log-likelihood: " << max_log_likelihood << std::endl;
    logger->info("Positions initialized from provided matrices: {} rows, {} cols, initial LL: {:.4f}", n_rows, n_cols, max_log_likelihood);
}

void DNA::initialize_structure_from_clusters(int sub_iterations, int sub_burn_in, double initial_sd) {
    if (cluster_adjacencies.empty()) {
        std::cerr << "Error: Build cluster relationships before initializing from clusters." << std::endl;
        return;
    }
    
    arma::uword num_clusters = cluster_adjacencies.size();
    initial_cluster_row_structures.resize(num_clusters);
    initial_cluster_col_structures.resize(num_clusters);
    initial_beta0s.resize(num_clusters);
    initial_beta1s.resize(num_clusters);
    
    // Create initialization directory
    std::string init_dir = dna_name + "/initialization";
    std::filesystem::create_directories(init_dir);
    
    std::cout << "\n--- Initializing Structures for Each Cluster and Backbone (in parallel) ---" << std::endl;
    auto init_start_time = std::chrono::steady_clock::now();
    
    #pragma omp parallel for
    for (arma::uword c = 0; c <= num_clusters; ++c) {
        if (c < num_clusters) {
            // --- Process Individual Clusters ---
            const arma::uvec& self_row_idx = cluster_adjacencies[c].self_row_indices;
            const arma::uvec& self_col_idx = cluster_adjacencies[c].self_col_indices;
            arma::mat sub_contacts = contact_matrix.submat(self_row_idx, self_col_idx);
            
            std::string cluster_name = init_dir + "/cluster_" + std::to_string(c);
            
            #pragma omp critical
            std::cout << "--- Processing cluster " << c << " (" << self_row_idx.n_elem << " rows x " 
                      << self_col_idx.n_elem << " cols) on thread " << omp_get_thread_num() << " ---" << std::endl;
            
            DNA cluster_dna(cluster_name);
            cluster_dna.contact_matrix = sub_contacts;
            cluster_dna.initialize_positions();
            cluster_dna.beta0 = 3.0;
            cluster_dna.beta1 = -2.0;
            cluster_dna.run_mcmc(sub_iterations, sub_burn_in, initial_sd);
            
            initial_cluster_row_structures[c] = cluster_dna.best_row_position_matrix;
            initial_cluster_col_structures[c] = cluster_dna.best_col_position_matrix;
            initial_beta0s[c] = cluster_dna.beta0;
            initial_beta1s[c] = cluster_dna.beta1;
            
            // Save results
            cluster_dna.best_row_position_matrix.save(cluster_name + "/final_row_positions.txt", arma::raw_ascii);
            cluster_dna.best_col_position_matrix.save(cluster_name + "/final_col_positions.txt", arma::raw_ascii);
        } else {
            // --- Process Backbone Structure ---
            std::string backbone_name = init_dir + "/backbone";
            
            #pragma omp critical
            std::cout << "--- Processing backbone (" << backbone_contact_matrix.n_rows 
                      << " clusters) on thread " << omp_get_thread_num() << " ---" << std::endl;
            
            DNA backbone_dna(backbone_name);
            backbone_dna.contact_matrix = backbone_contact_matrix;
            backbone_dna.initialize_positions();
            backbone_dna.beta0 = 1.0;
            backbone_dna.beta1 = -1.0;
            backbone_dna.run_mcmc(sub_iterations, sub_burn_in, initial_sd);
            
            #pragma omp critical
            {
                backbone_row_structure = backbone_dna.best_row_position_matrix;
                backbone_col_structure = backbone_dna.best_col_position_matrix;
                backbone_beta0 = backbone_dna.beta0;
                backbone_beta1 = backbone_dna.beta1;
                
                // Save results
                backbone_dna.best_row_position_matrix.save(backbone_name + "/final_row_positions.txt", arma::raw_ascii);
                backbone_dna.best_col_position_matrix.save(backbone_name + "/final_col_positions.txt", arma::raw_ascii);
            }
        }
    }
    
    std::cout << "--- Finished Initializing All Cluster and Backbone Structures ---" << std::endl;
    auto init_end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> init_duration = init_end_time - init_start_time;
    logger->info("Initialization completed in {:.2f} seconds.", init_duration.count());
}

void DNA::assemble_global_structure() {
    if (initial_cluster_row_structures.empty() || initial_cluster_col_structures.empty() || 
        backbone_row_structure.is_empty() || backbone_col_structure.is_empty()) {
        std::cerr << "Error: Must run initialize_structure_from_clusters() before assembling." << std::endl;
        return;
    }
    
    std::cout << "\n--- Assembling Global DNA Structure ---" << std::endl;
    auto assembly_start_time = std::chrono::steady_clock::now();
    
    // Calculate median beta values from all clusters
    arma::mat beta_mat(initial_beta0s.size(), 2);
    beta_mat.col(0) = arma::conv_to<arma::vec>::from(initial_beta0s);
    beta_mat.col(1) = arma::conv_to<arma::vec>::from(initial_beta1s);
    
    double median_beta0 = arma::median(beta_mat.col(0));
    double median_beta1 = arma::median(beta_mat.col(1));
    
    this->beta0 = median_beta0;
    this->beta1 = median_beta1;
    
    logger->info("Median beta0 for assembly: {:.4f}", median_beta0);
    logger->info("Median beta1 for assembly: {:.4f}", median_beta1);
    
    // --- Scale Backbone Structure using median beta values ---
    std::cout << "Scaling backbone structure..." << std::endl;
    
    // Calculate current distances in backbone
    arma::mat backbone_current_distances = calculate_pairwise_distances(backbone_row_structure, backbone_col_structure);
    
    // Calculate ratios for backbone scaling
    std::vector<double> backbone_ratios_vec;
    backbone_ratios_vec.reserve(backbone_contact_matrix.n_elem);
    
    for (arma::uword i = 0; i < backbone_contact_matrix.n_rows; ++i) {
        for (arma::uword j = 0; j < backbone_contact_matrix.n_cols; ++j) {
            if (backbone_contact_matrix(i, j) > 0 && backbone_current_distances(i, j) > 0) {
                // Using d^2 model: d = sqrt((log(contact) - beta0) / beta1)
                double d_squared_val = (std::log(backbone_contact_matrix(i, j)) - median_beta0) / median_beta1;
                if (d_squared_val > 0) {  // Only take sqrt if positive
                    double expected_dist = std::sqrt(d_squared_val);
                    if (std::isfinite(expected_dist)) {
                        backbone_ratios_vec.push_back(expected_dist / backbone_current_distances(i, j));
                    }
                }
            }
        }
    }
    
    // Use median of ratios as the scaling factor for backbone
    arma::vec backbone_ratios(backbone_ratios_vec);
    double backbone_scale = backbone_ratios.is_empty() ? 1.0 : arma::median(backbone_ratios);
    backbone_scale = std::max(0.5, std::min(2.0, backbone_scale)); // Clamp scale factor
    
    logger->info("Scaling backbone by a factor of {:.4f} (based on {} contact pairs)", backbone_scale, backbone_ratios_vec.size());
    
    // Apply scaling to backbone structures
    arma::rowvec backbone_center = arma::mean(backbone_row_structure, 0);
    backbone_row_structure.each_row() -= backbone_center;
    backbone_col_structure.each_row() -= backbone_center;
    backbone_row_structure *= backbone_scale;
    backbone_col_structure *= backbone_scale;
    backbone_row_structure.each_row() += backbone_center;
    backbone_col_structure.each_row() += backbone_center;
    
    arma::uword n_rows = contact_matrix.n_rows;
    arma::uword n_cols = contact_matrix.n_cols;
    arma::uword num_clusters = initial_cluster_row_structures.size();
    
    row_position_matrix.set_size(n_rows, 2);
    col_position_matrix.set_size(n_cols, 2);
    
    for (arma::uword c = 0; c < num_clusters; ++c) {
        arma::mat cluster_row_pos = initial_cluster_row_structures[c];
        arma::mat cluster_col_pos = initial_cluster_col_structures[c];
        const arma::uvec& cluster_row_indices = cluster_adjacencies[c].self_row_indices;
        const arma::uvec& cluster_col_indices = cluster_adjacencies[c].self_col_indices;
        
        // --- Scaling using median beta values ---
        // Get the submatrix of contacts for this cluster
        const arma::mat& cluster_contacts = contact_matrix.submat(cluster_row_indices, cluster_col_indices);
        
        // Calculate current distances between row and column positions
        arma::mat current_distances = calculate_pairwise_distances(cluster_row_pos, cluster_col_pos);
        
        // Calculate the ratio of expected distance to inferred distance for all pairs
        std::vector<double> ratios_vec;
        ratios_vec.reserve(cluster_contacts.n_elem);
        
        for (arma::uword i = 0; i < cluster_contacts.n_rows; ++i) {
            for (arma::uword j = 0; j < cluster_contacts.n_cols; ++j) {
                if (cluster_contacts(i, j) > 0 && current_distances(i, j) > 0) {
                    // Using d^2 model: lambda = exp(beta0 + beta1 * d^2)
                    // Solving for d: d = sqrt((log(contact) - beta0) / beta1)
                    double d_squared_val = (std::log(cluster_contacts(i, j)) - median_beta0) / median_beta1;
                    if (d_squared_val > 0) {  // Only take sqrt if positive
                        double expected_dist = std::sqrt(d_squared_val);
                        if (std::isfinite(expected_dist)) {
                            ratios_vec.push_back(expected_dist / current_distances(i, j));
                        }
                    }
                }
            }
        }
        
        // Use median of ratios as the scaling factor
        arma::vec ratios(ratios_vec);
        double scale = ratios.is_empty() ? 1.0 : arma::median(ratios);
        scale = std::max(0.5, std::min(2.0, scale)); // Clamp scale factor
        
        logger->info("Scaling cluster {} by a factor of {:.4f} (based on {} contact pairs)", c, scale, ratios_vec.size());
        
        cluster_row_pos.each_row() -= arma::mean(cluster_row_pos, 0);
        cluster_col_pos.each_row() -= arma::mean(cluster_col_pos, 0);

        arma::mat final_cluster_row_pos = cluster_row_pos * scale;
        arma::mat final_cluster_col_pos = cluster_col_pos * scale;
        final_cluster_row_pos.each_row() += backbone_row_structure.row(c);
        final_cluster_col_pos.each_row() += backbone_col_structure.row(c);
        // --- Final placement ---
        row_position_matrix.rows(cluster_row_indices) = final_cluster_row_pos;
        col_position_matrix.rows(cluster_col_indices) = final_cluster_col_pos;
    }
    

    // // ---- now rotate each cluster to maximize the loglikelihood ----
    // std::cout << "Refining cluster orientations via rotation..." << std::endl;
    // #pragma omp parallel for
    // for (arma::uword c = 0; c < num_clusters; ++c) {
    //     const arma::uvec& cluster_row_indices = cluster_adjacencies[c].self_row_indices;
    //     const arma::uvec& cluster_col_indices = cluster_adjacencies[c].self_col_indices;
    //     const arma::uvec& other_row 

    //     // Extract the positions for the current cluster
    //     arma::mat cluster_row_pos = row_position_matrix.rows(cluster_row_indices);
    //     arma::mat cluster_col_pos = col_position_matrix.rows(cluster_col_indices);
    //     arma::mat best_cluster_row_pos = cluster_row_pos;
    //     arma::mat best_cluster_col_pos = cluster_col_pos;

    //      for (int i = 0; i < num_angles; ++i) {
    //         double angle = (360.0 / num_angles) * i;
    //         arma::mat rotation = get_clockwise_rotation(angle);
    //         arma::mat rotated_row_pos = rotate_positions(cluster_row_pos, rotation);
    //         arma::mat rotated_col_pos = rotate_positions(cluster_col_pos, rotation);
            
    //      }

    // }



    std::cout << "--- Global Assembly Complete ---" << std::endl;
    auto assembly_end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> assembly_duration = assembly_end_time - assembly_start_time;
    logger->info("Assembly completed in {:.2f} seconds.", assembly_duration.count());
    
    // Recalculate pairwise distance matrix with assembled structure
    pairwise_distance_matrix = calculate_pairwise_distances(row_position_matrix, col_position_matrix);
    
    // Update best state
    max_log_likelihood = calculate_log_likelihood(pairwise_distance_matrix, contact_matrix);
    best_row_position_matrix = row_position_matrix;
    best_col_position_matrix = col_position_matrix;
    
    logger->info("Assembled structure log-likelihood: {:.4f}", max_log_likelihood);
}

void DNA::run_mcmc(int iterations, int burn_in, double initial_sd, double sd_floor, double sd_ceiling, bool save_samples, int sample_interval) {
    if (row_position_matrix.is_empty() || col_position_matrix.is_empty()) {
        std::cerr << "Error: Row and column positions must be initialized before running MCMC." << std::endl;
        return;
    }
    if (pairwise_distance_matrix.is_empty()) {
        std::cerr << "Error: Pairwise distance matrix must be initialized before running MCMC." << std::endl;
        return;
    }
    
    // Create intermediate results directory if saving is enabled
    std::string intermediate_dir;
    if (save_samples) {
        intermediate_dir = dna_name + "/intermediate_results";
        std::filesystem::create_directories(intermediate_dir);
        logger->info("Saving intermediate results to: {}", intermediate_dir);
    }
    
    double current_ll = calculate_log_likelihood(pairwise_distance_matrix, contact_matrix);
    std::cout << "\n--- Starting MCMC Sampling for DNA ---" << std::endl;
    logger->info("--- Starting MCMC Sampling for DNA ---");
    logger->info("Initial beta0: {:.4f}, Initial beta1: {:.4f}, Initial loglikelihood: {:.4f}", this->beta0, this->beta1, current_ll);
    logger->info("Iterations: {}, Burn-in: {}", iterations, burn_in);
    logger->info("Using d^2 likelihood model (lambda = exp(b0 + b1*d^2))");
    
    // --- Adaptive MCMC setup ---
    double sd_b0 = initial_sd, sd_b1 = initial_sd, sd_row = initial_sd, sd_col = initial_sd;
    int accepted_b0 = 0, accepted_b1 = 0, accepted_row = 0, accepted_col = 0;
    int total_accepted_b0 = 0, total_accepted_b1 = 0, total_accepted_row = 0, total_accepted_col = 0;
    
    arma::uword n_rows = row_position_matrix.n_rows;
    arma::uword n_cols = col_position_matrix.n_rows;
    
    auto block_start_time = std::chrono::steady_clock::now();
    auto total_start_time = std::chrono::steady_clock::now();
    double delta_ll = 0;
    
    for (int i = 0; i < iterations; ++i) {
        // --- Sample beta0 ---
        double proposed_beta0 = this->beta0 + arma::randn() * sd_b0;
        double old_beta0 = this->beta0;
        this->beta0 = proposed_beta0;
        delta_ll = calculate_log_likelihood(pairwise_distance_matrix, contact_matrix) - current_ll;
        if (arma::randu() < exp(delta_ll)) {
            accepted_b0++; total_accepted_b0++; current_ll += delta_ll; 
        } else { 
            this->beta0 = old_beta0; 
        }
        delta_ll = 0;
        
        // --- Sample beta1 (must be negative for d^2 model) ---
        double proposed_beta1;
        do { 
            proposed_beta1 = this->beta1 + arma::randn() * sd_b1; 
        } while (proposed_beta1 >= 0);
        double old_beta1 = this->beta1;
        this->beta1 = proposed_beta1;
        delta_ll = calculate_log_likelihood(pairwise_distance_matrix, contact_matrix) - current_ll;
        if (arma::randu() < exp(delta_ll)) {
            accepted_b1++; total_accepted_b1++; current_ll += delta_ll;
        } else { 
            this->beta1 = old_beta1; 
        }
        delta_ll = 0;
        
        // --- Sample Row Positions ---
        for (arma::uword r = 0; r < n_rows; ++r) {
            arma::rowvec current_row_pos = row_position_matrix.row(r);
            arma::rowvec proposed_row_pos = current_row_pos + arma::randn<arma::rowvec>(2) * sd_row;
            
            // Calculate proposed distances for this row to all columns
            arma::mat proposed_dists_r = calculate_pairwise_distances(proposed_row_pos, col_position_matrix);
            arma::mat current_dists_r = pairwise_distance_matrix.row(r);
            arma::mat contacts_r = contact_matrix.row(r);
            
            delta_ll = calculate_log_likelihood(proposed_dists_r, contacts_r)  - calculate_log_likelihood(current_dists_r, contacts_r);
            
            if (arma::randu() < exp(delta_ll)) {
                row_position_matrix.row(r) = proposed_row_pos;
                pairwise_distance_matrix.row(r) = proposed_dists_r;
                accepted_row++; total_accepted_row++; current_ll += delta_ll;
            }
            delta_ll = 0;
        }
        
        // --- Sample Column Positions ---
        for (arma::uword c = 0; c < n_cols; ++c) {
            arma::rowvec current_col_pos = col_position_matrix.row(c);
            arma::rowvec proposed_col_pos = current_col_pos + arma::randn<arma::rowvec>(2) * sd_col;
            
            // Calculate proposed distances from all rows to this column
            arma::mat proposed_dists_c = calculate_pairwise_distances(row_position_matrix, proposed_col_pos);
            arma::mat current_dists_c = pairwise_distance_matrix.col(c);
            arma::mat contacts_c = contact_matrix.col(c);
            
            delta_ll = calculate_log_likelihood(proposed_dists_c, contacts_c) - calculate_log_likelihood(current_dists_c, contacts_c);
            
            if (arma::randu() < exp(delta_ll)) {
                col_position_matrix.row(c) = proposed_col_pos;
                pairwise_distance_matrix.col(c) = proposed_dists_c;
                accepted_col++; total_accepted_col++; current_ll += delta_ll;
            }
            delta_ll = 0;
        }
        
        // --- Track best state ---
        if (current_ll > max_log_likelihood) {
            max_log_likelihood = current_ll;
            best_row_position_matrix = row_position_matrix;
            best_col_position_matrix = col_position_matrix;
        }
        
        // --- Adapt proposal SDs and log block progress ---
        if ((i + 1) % 50 == 0) {
            auto block_end_time = std::chrono::steady_clock::now();
            std::chrono::duration<double> block_duration = block_end_time - block_start_time;
            
            double rate_b0 = (double)accepted_b0 / 50.0;
            double rate_b1 = (double)accepted_b1 / 50.0;
            double rate_row = (double)accepted_row / (50.0 * n_rows);
            double rate_col = (double)accepted_col / (50.0 * n_cols);
            
            logger->info("-------------------- Block Summary (Iter {}) --------------------", i + 1);
            logger->info("Time for last 50 iterations: {:.2f} seconds.", block_duration.count());
            logger->info("Acceptance Rates -> beta0: {:.2f}%, beta1: {:.2f}%, row_pos: {:.2f}%, col_pos: {:.2f}%", 
                        rate_b0 * 100.0, rate_b1 * 100.0, rate_row * 100.0, rate_col * 100.0);
            logger->info("Proposal SDs -> beta0: {:.4f}, beta1: {:.4f}, row_pos: {:.4f}, col_pos: {:.4f}", 
                        sd_b0, sd_b1, sd_row, sd_col);
            logger->info("Current State -> max LL: {:.4f}, beta0: {:.4f}, beta1: {:.4f}", 
                        max_log_likelihood, this->beta0, this->beta1);
            
            // Save intermediate results if requested
            if (save_samples) {
                std::string iter_str = std::to_string(i + 1);
                best_row_position_matrix.save(intermediate_dir + "/row_positions_iter" + iter_str + ".txt", arma::raw_ascii);
                best_col_position_matrix.save(intermediate_dir + "/col_positions_iter" + iter_str + ".txt", arma::raw_ascii);
                
                // Save parameters
                std::ofstream param_file(intermediate_dir + "/parameters_iter" + iter_str + ".txt");
                param_file << "iteration " << (i + 1) << std::endl;
                param_file << "max_log_likelihood " << max_log_likelihood << std::endl;
                param_file << "current_log_likelihood " << current_ll << std::endl;
                param_file << "beta0 " << this->beta0 << std::endl;
                param_file << "beta1 " << this->beta1 << std::endl;
                param_file.close();
            }
            
            // Adapt proposal standard deviations
            if (rate_b0 > 0.3) sd_b0 = std::min(sd_b0*2.0, sd_ceiling); 
            else if (rate_b0 < 0.2) sd_b0 = std::max(sd_b0/2.0, sd_floor);
            
            if (rate_b1 > 0.3) sd_b1 = std::min(sd_b1*2.0, sd_ceiling); 
            else if (rate_b1 < 0.2) sd_b1 = std::max(sd_b1/2.0, sd_floor);
            
            if (rate_row > 0.3) sd_row = std::min(sd_row*2.0, sd_ceiling); 
            else if (rate_row < 0.2) sd_row = std::max(sd_row/2.0, sd_floor);
            
            if (rate_col > 0.3) sd_col = std::min(sd_col*2.0, sd_ceiling); 
            else if (rate_col < 0.2) sd_col = std::max(sd_col/2.0, sd_floor);
            
            accepted_b0 = 0; accepted_b1 = 0; accepted_row = 0; accepted_col = 0;
            block_start_time = std::chrono::steady_clock::now();
        }
        
        if ((i + 1) % 1000 == 0) {
            std::cout << "Iter " << i+1 << "/" << iterations << ". Current LL: " << current_ll << ". Best LL: " << max_log_likelihood << std::endl;
        }
    }
    
    std::cout << "--- MCMC Finished ---" << std::endl;
    std::cout << "Final max log-likelihood found: " << max_log_likelihood << std::endl;
    std::cout << "Beta0 acceptance rate: " << (double)total_accepted_b0 / iterations * 100.0 << "%" << std::endl;
    std::cout << "Beta1 acceptance rate: " << (double)total_accepted_b1 / iterations * 100.0 << "%" << std::endl;
    std::cout << "Row position acceptance rate: " << (double)total_accepted_row / (iterations * n_rows) * 100.0 << "%" << std::endl;
    std::cout << "Column position acceptance rate: " << (double)total_accepted_col / (iterations * n_cols) * 100.0 << "%" << std::endl;
    
    // Save final results
    std::string final_row_file = dna_name + "/final_row_positions.txt";
    std::string final_col_file = dna_name + "/final_col_positions.txt";
    std::string final_param_file = dna_name + "/final_parameters.txt";
    
    best_row_position_matrix.save(final_row_file, arma::raw_ascii);
    best_col_position_matrix.save(final_col_file, arma::raw_ascii);
    
    std::ofstream param_file(final_param_file);
    param_file << "# Final MCMC Results" << std::endl;
    param_file << "iterations " << iterations << std::endl;
    param_file << "burn_in " << burn_in << std::endl;
    param_file << "max_log_likelihood " << max_log_likelihood << std::endl;
    param_file << "beta0 " << this->beta0 << std::endl;
    param_file << "beta1 " << this->beta1 << std::endl;
    param_file << "n_rows " << n_rows << std::endl;
    param_file << "n_cols " << n_cols << std::endl;
    param_file << "beta0_acceptance_rate " << (double)total_accepted_b0 / iterations << std::endl;
    param_file << "beta1_acceptance_rate " << (double)total_accepted_b1 / iterations << std::endl;
    param_file << "row_pos_acceptance_rate " << (double)total_accepted_row / (iterations * n_rows) << std::endl;
    param_file << "col_pos_acceptance_rate " << (double)total_accepted_col / (iterations * n_cols) << std::endl;
    param_file.close();
    
    std::cout << "Final results saved to:" << std::endl;
    std::cout << "  Row positions: " << final_row_file << std::endl;
    std::cout << "  Col positions: " << final_col_file << std::endl;
    std::cout << "  Parameters: " << final_param_file << std::endl;
    
    auto total_end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> total_duration = total_end_time - total_start_time;
    logger->info("Total MCMC duration: {:.2f} seconds.", total_duration.count());
}

// Save assembled structures to files
void DNA::save_assembled_structures(const std::string& output_dir) {
    std::string base_dir = output_dir.empty() ? dna_name : output_dir;
    
    // Create output directory if it doesn't exist
    std::string mkdir_cmd = "mkdir -p " + base_dir;
    int ret = std::system(mkdir_cmd.c_str());
    (void)ret;
    
    // Save assembled row positions
    std::string row_file = base_dir + "/assembled_row_positions.txt";
    bool row_saved = row_position_matrix.save(row_file, arma::raw_ascii);
    
    // Save assembled column positions
    std::string col_file = base_dir + "/assembled_col_positions.txt";
    bool col_saved = col_position_matrix.save(col_file, arma::raw_ascii);
    
    // Save parameters
    std::string param_file_path = base_dir + "/assembled_parameters.txt";
    std::ofstream param_file(param_file_path);
    param_file << "# Assembled Structure Parameters" << std::endl;
    param_file << "beta0 " << this->beta0 << std::endl;
    param_file << "beta1 " << this->beta1 << std::endl;
    param_file << "n_rows " << row_position_matrix.n_rows << std::endl;
    param_file << "n_cols " << col_position_matrix.n_rows << std::endl;
    param_file << "max_log_likelihood " << max_log_likelihood << std::endl;
    param_file.close();
    
    if (row_saved && col_saved) {
        logger->info("Assembled structures saved to:");
        logger->info("  Row positions: {}", row_file);
        logger->info("  Col positions: {}", col_file);
        logger->info("  Parameters: {}", param_file_path);
        
        std::cout << "\n✓ Assembled structures saved:" << std::endl;
        std::cout << "  - " << row_file << std::endl;
        std::cout << "  - " << col_file << std::endl;
        std::cout << "  - " << param_file_path << std::endl;
    } else {
        logger->error("Failed to save assembled structures");
    }
}

// Getters
const arma::mat& DNA::get_row_position_matrix() const {
    return row_position_matrix;
}

const arma::mat& DNA::get_col_position_matrix() const {
    return col_position_matrix;
}

const arma::mat& DNA::get_best_row_position_matrix() const {
    return best_row_position_matrix;
}

const arma::mat& DNA::get_best_col_position_matrix() const {
    return best_col_position_matrix;
}
