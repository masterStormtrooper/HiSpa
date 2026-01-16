# HiSpaR Package Skeleton - Complete

## ✅ Package Structure Created

Your R package skeleton has been successfully created in `/Users/lyc/Proj/HiSpa/R-package/`

### Directory Structure

```
R-package/
├── DESCRIPTION              ✅ Package metadata and dependencies
├── NAMESPACE               ✅ Function exports
├── LICENSE                 ✅ MIT License
├── .Rbuildignore          ✅ Build exclusions
├── .gitignore             ✅ Git exclusions
├── README.md              ✅ Main documentation
├── SETUP_GUIDE.md         ✅ Detailed setup instructions
├── QUICK_REFERENCE.md     ✅ Quick reference guide
│
├── R/                      ✅ R source code
│   ├── HiSpaR-package.R       - Package documentation
│   ├── hispa_functions.R      - Main R wrapper functions
│   └── RcppExports.R          - Rcpp exports (auto-generated)
│
├── src/                    ✅ C++ source code
│   ├── Makevars               - Unix/macOS compilation settings
│   ├── Makevars.win           - Windows compilation settings
│   └── RcppExports_skeleton.cpp - Template for C++ wrappers
│
├── man/                    📝 Documentation (will be auto-generated)
│
├── inst/                   📦 Installed files
│   └── include/               - Header files (optional)
│
├── examples/               ✅ Example scripts
│   ├── basic_analysis.R       - Simple usage
│   ├── with_priors.R          - Prior information
│   └── convolution_example.R  - Convolution
│
├── tests/                  ✅ Unit tests
│   ├── testthat.R             - Test runner
│   └── testthat/
│       └── test-hispa.R       - Test cases
│
├── vignettes/              ✅ Long-form documentation
│   └── getting-started.Rmd    - Tutorial vignette
│
└── data/                   📊 Example datasets (optional)
```

## 📋 Main Functions Provided

### 1. `hispa_analyze()`
Main function for running HiSpa MCMC analysis.

**Parameters:**
- `contact_matrix` - Hi-C contact matrix (symmetric)
- `output_dir` - Output directory
- `mcmc_iterations` - Number of MCMC iterations
- `mcmc_burn_in` - Burn-in period
- `num_clusters` - Number of clusters (auto if 0)
- `verbose` - Enable detailed output

**Returns:** List with `position_matrix`, `beta0`, `beta1`, `log_likelihood`, `cluster_labels`

### 2. `hispa_analyze_with_prior()`
Analysis with gamma priors from prior position data.

**Parameters:**
- `contact_matrix` - Hi-C contact matrix
- `prior_positions` - Prior 3D positions (n × 3)
- `prior_distances` - Genomic distances for priors
- Other MCMC parameters

**Returns:** Same as above plus `priors` list

### 3. `convolute_contacts()`
Convolve contact matrix for smoothing.

**Parameters:**
- `contact_matrix` - Input matrix
- `half_k` - Half window size (default: 3)

**Returns:** Convoluted matrix

## 🔧 Next Steps to Complete Integration

### Step 1: Link C++ Implementation

You need to connect your existing C++ code to the R package. You have two main tasks:

#### A. Add Getter/Setter Methods to `chromosome.h`

Add these public methods to your `Chromosome` class:

```cpp
public:
    // Setters
    void set_contact_matrix(const arma::mat& contacts) {
        contact_matrix = contacts;
    }
    
    void set_skip_zero_contact_loci(bool skip) {
        skip_zero_contact_loci = skip;
    }
    
    void set_sample_from_prior(bool sample) {
        sample_from_prior = sample;
    }
    
    // Getters
    arma::mat get_position_matrix() const {
        return position_matrix;
    }
    
    double get_beta0() const {
        return beta0;
    }
    
    double get_beta1() const {
        return beta1;
    }
    
    double get_max_log_likelihood() const {
        return max_log_likelihood;
    }
    
    arma::uvec get_cluster_labels() const {
        return cluster_labels;
    }
```

#### B. Create Main Rcpp Wrapper File

Create `R-package/src/hispa_rcpp.cpp`:

```cpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include "../../include/chromosome.h"

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
Rcpp::List hispa_analyze(
    const arma::mat& contact_matrix,
    const std::string& output_dir,
    int mcmc_iterations = 6000,
    int mcmc_burn_in = 1000,
    int num_clusters = 0,
    double mcmc_initial_sd = 0.1,
    double mcmc_sd_floor = 0.0001,
    double mcmc_sd_ceil = 0.3,
    int cluster_distance_threshold = 2,
    bool save_samples = false,
    int sample_interval = 50,
    bool skip_zero_contact_loci = true,
    bool use_convoluted_sampling = false,
    int convolution_half_k = 3,
    bool verbose = true
) {
    // Create Chromosome object
    Chromosome chrom(output_dir);
    
    // Set contact matrix
    chrom.set_contact_matrix(contact_matrix);
    chrom.set_skip_zero_contact_loci(skip_zero_contact_loci);
    
    // Initialize positions
    chrom.initialize_positions();
    
    // Clustering
    if (num_clusters <= 0) {
        num_clusters = static_cast<int>(std::sqrt(contact_matrix.n_rows));
    }
    chrom.cluster_loci(num_clusters, cluster_distance_threshold);
    
    // Estimate initial beta
    chrom.estimate_beta();
    
    // Convolute if requested
    if (use_convoluted_sampling) {
        chrom.convolute_contacts(convolution_half_k);
    }
    
    // Run MCMC
    chrom.run_mcmc(
        mcmc_iterations,
        mcmc_burn_in,
        mcmc_initial_sd,
        mcmc_sd_floor,
        mcmc_sd_ceil,
        save_samples,
        sample_interval
    );
    
    // Get results
    return Rcpp::List::create(
        Named("position_matrix") = chrom.get_position_matrix(),
        Named("beta0") = chrom.get_beta0(),
        Named("beta1") = chrom.get_beta1(),
        Named("log_likelihood") = chrom.get_max_log_likelihood(),
        Named("cluster_labels") = chrom.get_cluster_labels(),
        Named("contact_matrix") = contact_matrix
    );
}

//' @export
// [[Rcpp::export]]
arma::mat convolute_contacts_R(const arma::mat& contact_matrix, int half_k = 3) {
    return convolute_contacts(contact_matrix, half_k);
}

// Add other wrapper functions as needed
```

#### C. Copy C++ Implementation

```bash
cd /Users/lyc/Proj/HiSpa/R-package/src
cp ../../src/chromosome.cpp .
```

### Step 2: Generate Rcpp Exports

```r
setwd("/Users/lyc/Proj/HiSpa/R-package")
Rcpp::compileAttributes()
devtools::document()
```

### Step 3: Build and Install

```r
# Check for issues
devtools::check()

# Install
devtools::install()

# Load and test
library(HiSpaR)
```

### Step 4: Test Installation

```r
library(HiSpaR)

# Create test data
test_mat <- matrix(rpois(400, 10), 20, 20)
test_mat <- (test_mat + t(test_mat)) / 2

# Run analysis
result <- hispa_analyze(
  contact_matrix = test_mat,
  output_dir = tempdir(),
  mcmc_iterations = 100,
  verbose = TRUE
)

print(result)
```

## 📚 Documentation Files

1. **README.md** - Overview, installation, quick start
2. **SETUP_GUIDE.md** - Detailed step-by-step setup instructions
3. **QUICK_REFERENCE.md** - Quick reference for common tasks
4. **examples/** - Working example scripts
5. **vignettes/** - Tutorial documentation

## 🔍 Important Notes

### Compilation Flags

The `Makevars` file is configured for:
- **C++17** standard
- **OpenMP** parallelization
- **Armadillo** linear algebra library
- Platform-specific optimizations (macOS/Linux)

### R Function Interface

All main functions include:
- Input validation
- Type conversion (R → C++)
- Error handling
- Documentation strings
- Return value formatting

### Class Structure

- R functions (`R/hispa_functions.R`) call C++ via Rcpp
- C++ wrappers (`src/hispa_rcpp.cpp`) call your `Chromosome` class
- Your existing C++ code remains largely unchanged

## 🐛 Troubleshooting

### If compilation fails:

1. **Check Armadillo:**
   ```bash
   pkg-config --modversion armadillo
   ```

2. **Check OpenMP:**
   ```bash
   echo | clang++ -fopenmp -xc++ -E - 2>&1 | grep -i openmp
   ```

3. **Update Makevars:**
   Edit `~/.R/Makevars` if needed

4. **Verify compiler:**
   ```bash
   clang++ --version  # Should support C++17
   ```

### Common Fixes:

**macOS OpenMP:** Add to `~/.R/Makevars`:
```make
PKG_CXXFLAGS = -Xclang -fopenmp
PKG_LIBS = -lomp
```

**Missing headers:** Copy to `inst/include/`:
```bash
cp -r ../include/* inst/include/
```

## 📞 Support

- **Setup Guide:** See `SETUP_GUIDE.md` for detailed instructions
- **Quick Ref:** See `QUICK_REFERENCE.md` for common tasks
- **Examples:** See `examples/` directory
- **Issues:** https://github.com/masterStormtrooper/HiSpa/issues

## ✨ What You Get

Once installed, users can:

```r
library(HiSpaR)

# Load Hi-C data
contacts <- as.matrix(read.table("contacts.txt"))

# Run analysis with one function call
result <- hispa_analyze(contacts, "output")

# Get 3D structure
positions <- result$position_matrix

# Visualize
plot3d(positions, col = rainbow(nrow(positions)))
```

All the power of your C++ implementation, accessible from R! 🚀
