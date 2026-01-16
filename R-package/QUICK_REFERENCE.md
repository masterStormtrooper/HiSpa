# Quick Reference: Converting HiSpa C++ to R Package

## Files Created

✅ **Package Structure**
- `DESCRIPTION` - Package metadata
- `NAMESPACE` - Export declarations
- `LICENSE` - MIT license
- `.Rbuildignore` - Build exclusions
- `.gitignore` - Git exclusions

✅ **R Interface** (`R/`)
- `HiSpaR-package.R` - Package documentation
- `hispa_functions.R` - R wrapper functions
- `RcppExports.R` - Auto-generated Rcpp exports

✅ **C++ Integration** (`src/`)
- `Makevars` - Unix/macOS compilation settings
- `Makevars.win` - Windows compilation settings
- `RcppExports_skeleton.cpp` - Template for C++ wrappers

✅ **Examples** (`examples/`)
- `basic_analysis.R` - Simple usage example
- `with_priors.R` - Prior information example
- `convolution_example.R` - Convolution example

✅ **Documentation**
- `README.md` - Main package documentation
- `SETUP_GUIDE.md` - Step-by-step setup instructions

✅ **Tests** (`tests/`)
- `testthat.R` - Test runner
- `testthat/test-hispa.R` - Unit tests

✅ **Vignette** (`vignettes/`)
- `getting-started.Rmd` - Tutorial vignette

## Next Steps

### 1. Copy/Link C++ Source Files

```bash
cd R-package/src

# Option A: Copy C++ implementation
cp ../../src/chromosome.cpp .

# Option B: Headers are already linked via Makevars
# The -I../../include flag makes your headers accessible
```

### 2. Create Main Rcpp Wrapper

Create `R-package/src/hispa_rcpp.cpp`:

```cpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include "../../include/chromosome.h"

// Add wrapper functions here
// See RcppExports.cpp for template
```

Key functions to wrap:
- `Chromosome::run_mcmc()`
- `Chromosome::get_position_matrix()`
- `Chromosome::get_beta0()` / `get_beta1()`
- `convolute_contacts()`

### 3. Add Getter Methods to Chromosome Class

Add to `include/chromosome.h` if not present:

```cpp
public:
    // Getters for R interface
    void set_contact_matrix(const arma::mat& contacts);
    arma::uvec get_cluster_labels() const;
    double get_max_log_likelihood() const;
    arma::mat get_position_matrix() const;
    double get_beta0() const;
    double get_beta1() const;
```

### 4. Generate Rcpp Exports

```r
setwd("/path/to/HiSpa/R-package")
Rcpp::compileAttributes()
devtools::document()
```

### 5. Build and Install

```r
# Check for issues
devtools::check()

# Install
devtools::install()

# Test
library(HiSpaR)
?hispa_analyze
```

## Key Modifications Needed

### In `chromosome.h`:
1. Add getter/setter methods for R interface
2. Ensure all needed methods are public
3. Consider adding simpler constructors

### In `chromosome.cpp`:
1. May need to remove `main()` if compiling together
2. Ensure spdlog initialization works in R context
3. Handle file I/O paths appropriately

### Create `hispa_rcpp.cpp`:
```cpp
//' @export
// [[Rcpp::export]]
Rcpp::List hispa_analyze(
    const arma::mat& contact_matrix,
    const std::string& output_dir,
    int mcmc_iterations = 6000,
    // ... other parameters
) {
    Chromosome chrom(output_dir);
    chrom.set_contact_matrix(contact_matrix);
    chrom.initialize_positions();
    // ... run analysis
    
    return Rcpp::List::create(
        Rcpp::Named("position_matrix") = chrom.get_position_matrix(),
        Rcpp::Named("beta0") = chrom.get_beta0(),
        // ...
    );
}
```

## Testing

```r
library(HiSpaR)

# Simple test
mat <- matrix(rpois(400, 10), 20, 20)
mat <- (mat + t(mat)) / 2

result <- hispa_analyze(
  contact_matrix = mat,
  output_dir = tempdir(),
  mcmc_iterations = 100,
  verbose = TRUE
)

print(result)
```

## Troubleshooting

### Compilation Errors
- Check compiler supports C++17
- Verify Armadillo installation: `pkg-config --modversion armadillo`
- Check OpenMP: Add flags to `~/.R/Makevars`

### Linking Errors
- Ensure all C++ source files are listed in Makevars
- Add explicit library links in PKG_LIBS

### Runtime Errors
- Check file paths in C++ code
- Verify spdlog works in R context
- Test with small matrices first

## Resources

- [Rcpp Gallery](https://gallery.rcpp.org/)
- [RcppArmadillo Examples](https://github.com/RcppCore/RcppArmadillo)
- [R Packages Book](https://r-pkgs.org/)
- [Writing R Extensions](https://cran.r-project.org/doc/manuals/R-exts.html)
