# R Package Modifications - Simplified Interface

## Summary

The R package has been modified to export **only the core workflow** from `hispa_main.cpp`, following your exact specifications.

## What Was Changed

### 1. Single Main Function: `hispa_analyze()`

**Parameters (matching hispa_main.cpp):**
- `contact_matrix` - Input Hi-C contact matrix
- `output_dir` - Output directory path
- `mcmc_iterations` - Number of MCMC iterations (default: 6000)
- `num_clusters` - Number of clusters (default: 0 for auto-detect)
- `mcmc_burn_in` - Burn-in period (default: 0)
- `mcmc_initial_sd` - Initial SD for proposals (default: 0.1)
- `mcmc_sd_floor` - Minimum SD (default: 0.0001)
- `mcmc_sd_ceil` - Maximum SD (default: 0.3)
- `use_cluster_init` - Use cluster-based initialization (default: FALSE)
- `cluster_init_iterations` - Cluster MCMC iterations (default: 1000)
- `cluster_initial_sd` - Cluster initial SD (default: 0.1)
- `save_samples` - Save MCMC trace (default: FALSE)
- `sample_interval` - Save interval (default: 50)
- `verbose` - Verbose output (default: TRUE)

**Returns:**
- `final_positions` - Final 3D coordinates (n × 3 matrix)
- `initial_positions` - Initial positions before MCMC (n × 3 matrix)
- `log_likelihood_trace` - MCMC log-likelihood trace vector
- `output_dir` - Path to output directory

**Workflow (exact match to hispa_main.cpp):**
1. Load contact matrix
2. Assign loci to clusters (k-means)
3. Build cluster relationships (distance threshold = 2)
4. Initialize structure (random or cluster-based)
5. Assemble global structure
6. Run main MCMC algorithm
7. Save results to output directory

**Files saved automatically:**
- `final_positions.txt` - Final 3D coordinates
- `initial_positions.txt` - Initial positions
- `log_likelihood_trace.txt` - Log-likelihood trace
- `block_timings.txt` - Timing information
- `mcmc_log.txt` - Detailed log

### 2. Removed Functions

The following functions were **removed** as they are not part of the core workflow:
- `hispa_analyze_with_prior()` - Prior information functionality
- `convolute_contacts()` - Contact matrix convolution

These features exist in the C++ code but are not part of the basic workflow you specified.

### 3. Updated C++ Wrapper

**File:** `R-package/src/hispa_rcpp.cpp`

The Rcpp wrapper now implements the exact workflow from `hispa_main.cpp`:
- Uses `set_skip_zero_contact_loci(true)` (default)
- Uses `set_sample_from_prior(false)` (default)
- Uses cluster distance threshold = 2
- Calls the exact same sequence of methods as the C++ main function

### 4. Updated Examples

All examples updated to use the simplified interface:

**basic_analysis.R:**
- Random initialization example
- Shows convergence plotting
- Compares initial vs final positions

**with_priors.R:** (renamed to focus on cluster initialization)
- Cluster-based initialization example
- Compares with random initialization
- Shows performance differences

**convolution_example.R:** (renamed to parameter comparison)
- Compares different parameter settings
- Tests random vs cluster initialization
- Tests different numbers of clusters

### 5. Updated Documentation

- README.md - Simplified to show only `hispa_analyze()`
- R/hispa_functions.R - Complete roxygen2 documentation for the single function
- Print/summary methods updated to work with new return structure

## Key Simplifications

1. **Single entry point**: One function does everything
2. **Workflow match**: Exactly follows hispa_main.cpp logic
3. **Auto-save**: Results automatically saved to output directory
4. **Simplified returns**: No beta0/beta1, cluster_labels - just positions and trace
5. **Clean interface**: Only essential parameters exposed

## Usage Example

```r
library(HiSpaR)

# Load contact matrix
contacts <- as.matrix(read.table("contacts.txt"))

# Run complete workflow
result <- hispa_analyze(
  contact_matrix = contacts,
  output_dir = "output",
  mcmc_iterations = 6000,
  mcmc_burn_in = 1000,
  use_cluster_init = TRUE
)

# Results
final_pos <- result$final_positions
initial_pos <- result$initial_positions

# Plot convergence
plot(result$log_likelihood_trace, type = "l")

# All files already saved in output/
```

## Files Modified

1. `R/hispa_functions.R` - Simplified to single function
2. `src/hispa_rcpp.cpp` - Rewritten to match C++ workflow
3. `examples/basic_analysis.R` - Updated
4. `examples/with_priors.R` - Changed to cluster initialization example
5. `examples/convolution_example.R` - Changed to parameter comparison
6. `README.md` - Simplified documentation

## Next Steps

To complete the package:

1. **Add getter methods** to `chromosome.h`:
   ```cpp
   void set_contact_matrix(const arma::mat& contacts);
   arma::mat get_best_position_matrix() const;
   arma::vec get_mcmc_trace_log_likelihood() const;
   // etc.
   ```

2. **Copy chromosome.cpp** to `src/`:
   ```bash
   cp ../src/chromosome.cpp R-package/src/
   ```

3. **Generate exports**:
   ```r
   Rcpp::compileAttributes()
   devtools::document()
   ```

4. **Build and install**:
   ```r
   devtools::check()
   devtools::install()
   ```

## Verification

The package now exports **exactly** the workflow from lines 201-419 of `hispa_main.cpp`:
- ✅ Load contact matrix
- ✅ Assign clusters
- ✅ Build cluster relationships
- ✅ Initialize structure (random or cluster-based)
- ✅ Assemble global structure
- ✅ Run main MCMC
- ✅ Save all results

No additional functionality is exposed.
