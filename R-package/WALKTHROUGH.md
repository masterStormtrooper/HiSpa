# Complete Walkthrough: Converting HiSpa C++ to R Package

This is your comprehensive guide to converting your HiSpa C++ program into an R package using Rcpp.

## 🎯 What We're Building

We're creating **HiSpaR**, an R package that wraps your existing C++ HiSpa code, allowing R users to:
- Load Hi-C contact matrices directly in R
- Run HiSpa analysis with a single function call
- Get 3D chromatin structures back as R matrices
- Use all HiSpa features from R environment

## ✅ What's Already Done

I've created a complete R package skeleton with:

### 📦 Core Package Files
- `DESCRIPTION` - Package metadata
- `NAMESPACE` - Function exports  
- `LICENSE` - MIT license
- Build configuration files

### 💻 R Interface
- `R/hispa_functions.R` - Three main R functions:
  - `hispa_analyze()` - Main analysis function
  - `hispa_analyze_with_prior()` - Analysis with priors
  - `convolute_contacts()` - Matrix convolution
- Full parameter documentation
- Input validation
- Pretty print methods

### 🔧 Compilation Setup
- `src/Makevars` - Unix/macOS compilation flags
- `src/Makevars.win` - Windows compilation flags
- C++17, OpenMP, and Armadillo configured

### 📚 Documentation
- `README.md` - User-facing documentation
- `SETUP_GUIDE.md` - Detailed setup instructions
- `QUICK_REFERENCE.md` - Quick reference
- `STRUCTURE.md` - Package structure overview
- `vignettes/getting-started.Rmd` - Tutorial

### 🎨 Examples
- `examples/basic_analysis.R` - Simple usage
- `examples/with_priors.R` - Prior information
- `examples/convolution_example.R` - Convolution

### ✓ Tests
- `tests/testthat/test-hispa.R` - Unit tests

## 🔨 What You Need to Do

### Step 1: Add Methods to Your Chromosome Class

Your `Chromosome` class needs getter/setter methods for R integration.

**Edit: `include/chromosome.h`**

Add these public methods (if they don't already exist):

```cpp
public:
    // Setters for R interface
    void set_contact_matrix(const arma::mat& contacts) {
        contact_matrix = contacts;
    }
    
    void set_skip_zero_contact_loci(bool skip) {
        skip_zero_contact_loci = skip;
    }
    
    void set_sample_from_prior(bool sample) {
        sample_from_prior = sample;
    }
    
    void set_prior_positions(const arma::mat& prior_pos) {
        prior_position_matrix = prior_pos;
    }
    
    void set_position_matrix(const arma::mat& positions) {
        position_matrix = positions;
    }
    
    // Getters for R interface
    arma::mat get_position_matrix() const {
        return position_matrix;
    }
    
    arma::mat get_contact_matrix() const {
        return contact_matrix;
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

### Step 2: Create Rcpp Wrapper File

This is the bridge between R and your C++ code.

**Create: `R-package/src/hispa_rcpp.cpp`**

```cpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <omp.h>

// Include your HiSpa headers - adjust path as needed
#include "../../include/chromosome.h"

using namespace Rcpp;

//' Run HiSpa MCMC Analysis
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
    try {
        // Create Chromosome object
        Chromosome chrom(output_dir);
        
        // Set contact matrix
        chrom.set_contact_matrix(contact_matrix);
        chrom.set_skip_zero_contact_loci(skip_zero_contact_loci);
        
        if (verbose) {
            Rcpp::Rcout << "Initializing positions...\n";
        }
        chrom.initialize_positions();
        
        // Clustering
        if (num_clusters <= 0) {
            num_clusters = static_cast<int>(std::sqrt(contact_matrix.n_rows));
        }
        
        if (verbose) {
            Rcpp::Rcout << "Clustering into " << num_clusters << " clusters...\n";
        }
        chrom.cluster_loci(num_clusters, cluster_distance_threshold);
        
        // Estimate initial beta parameters
        if (verbose) {
            Rcpp::Rcout << "Estimating initial parameters...\n";
        }
        chrom.estimate_beta();
        
        // Convolute if requested
        if (use_convoluted_sampling) {
            if (verbose) {
                Rcpp::Rcout << "Convoluting contact matrix...\n";
            }
            chrom.convolute_contacts(convolution_half_k);
        }
        
        // Run MCMC
        if (verbose) {
            Rcpp::Rcout << "Running MCMC (" << mcmc_iterations << " iterations)...\n";
        }
        
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
        arma::mat positions = chrom.get_position_matrix();
        double beta0 = chrom.get_beta0();
        double beta1 = chrom.get_beta1();
        double log_likelihood = chrom.get_max_log_likelihood();
        arma::uvec cluster_labels = chrom.get_cluster_labels();
        
        // Save results
        chrom.save_position_matrix_to_file(output_dir + "/position_matrix.txt");
        
        if (verbose) {
            Rcpp::Rcout << "Complete! Log-likelihood: " << log_likelihood << "\n";
        }
        
        // Return as R list
        return Rcpp::List::create(
            Named("position_matrix") = positions,
            Named("beta0") = beta0,
            Named("beta1") = beta1,
            Named("log_likelihood") = log_likelihood,
            Named("cluster_labels") = cluster_labels,
            Named("contact_matrix") = contact_matrix
        );
        
    } catch (const std::exception& e) {
        Rcpp::stop("Error in hispa_analyze: " + std::string(e.what()));
    }
}


//' Run HiSpa Analysis with Prior Information
//' @export
// [[Rcpp::export]]
Rcpp::List hispa_analyze_with_prior(
    const arma::mat& contact_matrix,
    const arma::mat& prior_positions,
    const std::vector<int>& prior_distances,
    const std::string& output_dir,
    int mcmc_iterations = 6000,
    int mcmc_burn_in = 1000,
    bool use_prior_positions = false,
    bool sample_from_prior = false,
    bool verbose = true
) {
    try {
        Chromosome chrom(output_dir);
        
        chrom.set_contact_matrix(contact_matrix);
        
        if (verbose) {
            Rcpp::Rcout << "Fitting gamma priors from prior positions...\n";
        }
        chrom.set_prior_positions(prior_positions);
        chrom.fit_gamma_priors(prior_distances);
        
        // Initialize positions
        if (use_prior_positions) {
            if (verbose) {
                Rcpp::Rcout << "Using prior positions as initial structure...\n";
            }
            chrom.set_position_matrix(prior_positions);
        } else {
            chrom.initialize_positions();
        }
        
        chrom.set_sample_from_prior(sample_from_prior);
        
        // Clustering
        int num_clusters = static_cast<int>(std::sqrt(contact_matrix.n_rows));
        if (verbose) {
            Rcpp::Rcout << "Clustering into " << num_clusters << " clusters...\n";
        }
        chrom.cluster_loci(num_clusters, 2);
        
        // Estimate beta
        if (verbose) {
            Rcpp::Rcout << "Estimating parameters...\n";
        }
        chrom.estimate_beta();
        
        // Run MCMC
        if (verbose) {
            Rcpp::Rcout << "Running MCMC with priors...\n";
        }
        chrom.run_mcmc(mcmc_iterations, mcmc_burn_in, 0.1, 0.0001, 0.3, false, 50);
        
        // Get results
        arma::mat positions = chrom.get_position_matrix();
        double beta0 = chrom.get_beta0();
        double beta1 = chrom.get_beta1();
        double log_likelihood = chrom.get_max_log_likelihood();
        arma::uvec cluster_labels = chrom.get_cluster_labels();
        
        chrom.save_position_matrix_to_file(output_dir + "/position_matrix.txt");
        
        if (verbose) {
            Rcpp::Rcout << "Complete!\n";
        }
        
        return Rcpp::List::create(
            Named("position_matrix") = positions,
            Named("beta0") = beta0,
            Named("beta1") = beta1,
            Named("log_likelihood") = log_likelihood,
            Named("cluster_labels") = cluster_labels
        );
        
    } catch (const std::exception& e) {
        Rcpp::stop("Error in hispa_analyze_with_prior: " + std::string(e.what()));
    }
}


//' Convolute Hi-C Contact Matrix
//' @export
// [[Rcpp::export]]
arma::mat convolute_contacts_R(const arma::mat& contact_matrix, int half_k = 3) {
    return convolute_contacts(contact_matrix, half_k);
}
```

### Step 3: Copy C++ Implementation

```bash
cd /Users/lyc/Proj/HiSpa/R-package/src
cp ../../src/chromosome.cpp .
```

**Optional:** Copy headers to `inst/include/` for cleaner distribution:

```bash
cd /Users/lyc/Proj/HiSpa/R-package
mkdir -p inst/include
cp -r ../include/* inst/include/
```

Then update `src/Makevars` to use `../inst/include` instead of `../../include`.

### Step 4: Generate Rcpp Exports

Open R in the package directory:

```bash
cd /Users/lyc/Proj/HiSpa/R-package
R
```

In R:

```r
# Install/load required packages
install.packages(c("Rcpp", "RcppArmadillo", "devtools", "roxygen2"))

library(devtools)

# Generate Rcpp exports from your [[Rcpp::export]] annotations
Rcpp::compileAttributes()

# Generate documentation from roxygen2 comments
document()
```

This will:
- Create/update `src/RcppExports.cpp`
- Create/update `R/RcppExports.R`
- Update `NAMESPACE`
- Generate `.Rd` files in `man/`

### Step 5: Build and Check

```r
# Check for errors and warnings
check()

# If check passes, install
install()
```

Fix any errors that appear. Common issues:
- Missing getter/setter methods
- Include path problems
- OpenMP linking issues (see troubleshooting)

### Step 6: Test the Package

```r
library(HiSpaR)

# Help documentation
?hispa_analyze

# Create test data
set.seed(123)
n <- 30
contact_mat <- matrix(rpois(n*n, lambda = 10), n, n)
contact_mat <- (contact_mat + t(contact_mat)) / 2

# Run analysis
result <- hispa_analyze(
  contact_matrix = contact_mat,
  output_dir = tempdir(),
  mcmc_iterations = 500,  # Small for testing
  mcmc_burn_in = 100,
  verbose = TRUE
)

# Check results
print(result)
summary(result)

# Verify output
str(result)
dim(result$position_matrix)  # Should be 30 x 3
```

## 📊 Complete Example Usage

Once installed, here's how users will use your package:

```r
library(HiSpaR)

# Load real Hi-C data
contact_matrix <- as.matrix(read.table("chr1_contacts.txt"))

# Run HiSpa analysis
result <- hispa_analyze(
  contact_matrix = contact_matrix,
  output_dir = "chr1_output",
  mcmc_iterations = 6000,
  mcmc_burn_in = 1000,
  num_clusters = 0,  # Auto-detect
  verbose = TRUE
)

# View results
print(result)
# HiSpa Analysis Result
# =====================
# Number of loci: 500
# Beta0: 5.234
# Beta1: -0.876
# Log-likelihood: -12345.67

# Extract 3D positions
positions <- result$position_matrix
write.table(positions, "chr1_3d_structure.txt")

# Visualize (requires rgl)
library(rgl)
plot3d(positions, 
       col = rainbow(nrow(positions)),
       size = 5,
       xlab = "X", ylab = "Y", zlab = "Z")

# Connect sequential loci
for (i in 1:(nrow(positions)-1)) {
  lines3d(positions[c(i, i+1), ], col = "gray")
}
```

## 🔧 Troubleshooting

### Issue: "chromosome.h: No such file or directory"

**Solution:** Check include paths in `src/Makevars`:
```make
PKG_CPPFLAGS = -I../inst/include -I../../include
```

### Issue: OpenMP errors on macOS

**Solution:** Create/edit `~/.R/Makevars`:
```make
PKG_CXXFLAGS = -Xclang -fopenmp
PKG_LIBS = -L/usr/local/opt/libomp/lib -lomp
CPPFLAGS = -I/usr/local/opt/libomp/include
```

### Issue: "undefined symbol" errors

**Solution:** Ensure `chromosome.cpp` is being compiled. Check that it's in `src/` directory.

### Issue: Armadillo not found

**Solution:**
```bash
# Install Armadillo
brew install armadillo  # macOS

# Verify
pkg-config --modversion armadillo
```

### Issue: spdlog errors in R context

**Solution:** The logger initialization might need adjustment. You may need to modify the `Chromosome` constructor to handle R environment better, or disable file logging when called from R.

## 🎓 Development Workflow

Once set up, iterate with:

```r
# 1. Edit C++ or R code

# 2. Recompile (if changed C++ exports)
Rcpp::compileAttributes()

# 3. Reload package without reinstalling
devtools::load_all()

# 4. Test
result <- hispa_analyze(test_matrix, tempdir(), mcmc_iterations=100)

# 5. When ready, full check and install
devtools::check()
devtools::install()
```

## 📚 Documentation Files Guide

- **README.md** - Start here for overview and installation
- **SETUP_GUIDE.md** - Detailed step-by-step instructions (this file)
- **QUICK_REFERENCE.md** - Common tasks quick reference
- **STRUCTURE.md** - Complete package structure
- **PACKAGE_SUMMARY.md** - Summary and checklist

## ✨ What Users Get

Your users will be able to:

```r
# Install
devtools::install_github("masterStormtrooper/HiSpa", subdir="R-package")

# Use
library(HiSpaR)
result <- hispa_analyze(contacts, "output")
plot3d(result$position_matrix)
```

All the power of your C++ HiSpa implementation, seamlessly integrated into R! 🚀

## 🎯 Summary Checklist

- ✅ Package skeleton created
- 🔲 Add getter/setter methods to `Chromosome` class
- 🔲 Create `src/hispa_rcpp.cpp` wrapper
- 🔲 Copy `chromosome.cpp` to `src/`
- 🔲 Run `Rcpp::compileAttributes()`
- 🔲 Run `devtools::document()`
- 🔲 Run `devtools::check()`
- 🔲 Run `devtools::install()`
- 🔲 Test with example data
- 🔲 Celebrate! 🎉

**You're ready to go!** Follow steps 1-6 above to complete the integration.
