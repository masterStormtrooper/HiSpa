# HiSpaR Setup Guide

This guide walks you through setting up the HiSpaR R package from your existing HiSpa C++ code.

## Step-by-Step Setup Process

### Step 1: Verify Prerequisites

#### 1.1 Check R Installation
```r
# In R console
R.version.string
# Should be R version 4.0.0 or higher
```

#### 1.2 Check Compiler
```bash
# macOS/Linux
clang++ --version  # or g++ --version
# Should support C++17

# Check for OpenMP
echo '#include <omp.h>' | clang++ -fopenmp -x c++ -E - >/dev/null 2>&1 && echo "OpenMP available" || echo "OpenMP not found"
```

#### 1.3 Install Armadillo
```bash
# macOS
brew install armadillo

# Verify installation
pkg-config --modversion armadillo
pkg-config --cflags armadillo
```

### Step 2: Install R Dependencies

```r
# Install required R packages
install.packages(c("Rcpp", "RcppArmadillo", "devtools", "roxygen2"))

# Verify installation
library(Rcpp)
library(RcppArmadillo)
packageVersion("Rcpp")
packageVersion("RcppArmadillo")
```

### Step 3: Prepare C++ Source Code

The R package needs access to your C++ source files. You have two options:

#### Option A: Copy Headers to inst/include (Recommended for distribution)
```bash
cd /path/to/HiSpa/R-package
mkdir -p inst/include
cp ../include/chromosome.h inst/include/
cp -r ../include/spdlog inst/include/
cp ../include/CLI11.hpp inst/include/
```

#### Option B: Link to Existing Headers (Development)
The `src/Makevars` already points to `../../include`, so your existing headers will be used.

### Step 4: Create Rcpp Wrapper

You need to create a C++ file that wraps your existing `Chromosome` class for R:

```bash
cd /path/to/HiSpa/R-package/src
```

Create `hispa_rcpp.cpp` with Rcpp exports:

```cpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include "../../include/chromosome.h"

// Your wrapper functions here (see RcppExports.cpp template)
```

### Step 5: Add Missing C++ Methods (if needed)

Your `Chromosome` class may need additional getter/setter methods for R integration. Add to `chromosome.h`:

```cpp
// In chromosome.h, add public methods:
void set_contact_matrix(const arma::mat& contacts) {
    contact_matrix = contacts;
}

arma::uvec get_cluster_labels() const {
    return cluster_labels;
}

double get_max_log_likelihood() const {
    return max_log_likelihood;
}
```

Then implement in `chromosome.cpp` if needed.

### Step 6: Generate Rcpp Exports

```r
# In R, from the R-package directory
setwd("/path/to/HiSpa/R-package")

# This generates RcppExports.R and updates NAMESPACE
Rcpp::compileAttributes()
```

### Step 7: Generate Documentation

```r
# Generate man/ documentation from roxygen2 comments
devtools::document()

# This creates .Rd files in man/ directory
```

### Step 8: Build and Check Package

```r
# Check for errors
devtools::check()

# Build the package
devtools::build()

# Install locally
devtools::install()
```

### Step 9: Test Installation

```r
# Load the package
library(HiSpaR)

# Test with simple data
test_matrix <- matrix(rpois(100, 5), 10, 10)
test_matrix <- (test_matrix + t(test_matrix)) / 2

# Run analysis
result <- hispa_analyze(
  contact_matrix = test_matrix,
  output_dir = tempdir(),
  mcmc_iterations = 100,  # Small for testing
  verbose = TRUE
)

print(result)
```

## Common Issues and Solutions

### Issue 1: "chromosome.h not found"

**Solution:** Check `src/Makevars` includes correct path:
```make
PKG_CPPFLAGS = -I../inst/include -I../../include
```

### Issue 2: "undefined symbol: _ZN10Chromosome..."

**Problem:** Missing implementation files.

**Solution:** Add `chromosome.cpp` to compilation:

1. Copy source file:
```bash
cp ../src/chromosome.cpp R-package/src/
```

2. Update `Makevars`:
```make
SOURCES = hispa_rcpp.cpp chromosome.cpp
OBJECTS = $(SOURCES:.cpp=.o)
```

### Issue 3: OpenMP not found on macOS

**Solution:** Edit `~/.R/Makevars`:
```make
PKG_CXXFLAGS = -Xclang -fopenmp
PKG_LIBS = -L/usr/local/opt/libomp/lib -lomp
CPPFLAGS = -I/usr/local/opt/libomp/include
```

### Issue 4: Armadillo linking errors

**Solution:** Explicitly link Armadillo in `Makevars`:
```make
PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -larmadillo
```

## Directory Structure Overview

```
R-package/
├── DESCRIPTION          # Package metadata
├── NAMESPACE           # Exported functions (auto-generated)
├── LICENSE             # MIT license
├── README.md           # Package documentation
├── .Rbuildignore       # Files to exclude from build
├── .gitignore          # Git ignore patterns
├── R/
│   ├── HiSpaR-package.R      # Package documentation
│   ├── hispa_functions.R     # R wrapper functions
│   └── RcppExports.R         # Auto-generated Rcpp exports
├── src/
│   ├── Makevars              # Unix compilation flags
│   ├── Makevars.win          # Windows compilation flags
│   ├── hispa_rcpp.cpp        # Rcpp wrapper code
│   ├── chromosome.cpp        # Your C++ implementation (copy)
│   └── RcppExports.cpp       # Auto-generated Rcpp code
├── inst/
│   └── include/
│       ├── chromosome.h      # Your headers (optional copy)
│       └── spdlog/          # Dependencies
├── man/                      # Documentation (auto-generated)
├── examples/                 # Example R scripts
├── tests/                    # Unit tests (optional)
└── vignettes/               # Long-form documentation (optional)
```

## Next Steps

1. **Add Unit Tests**: Create tests in `tests/testthat/`
2. **Write Vignettes**: Add tutorials in `vignettes/`
3. **Create Sample Data**: Add example datasets in `data/`
4. **CI/CD Setup**: Configure GitHub Actions for automated testing
5. **CRAN Submission**: Prepare for CRAN if desired

## Workflow for Development

```r
# 1. Make changes to R or C++ code

# 2. Recompile attributes (if changed C++ exports)
Rcpp::compileAttributes()

# 3. Regenerate documentation
devtools::document()

# 4. Reload package
devtools::load_all()

# 5. Test changes
# ... run your tests ...

# 6. Check package
devtools::check()

# 7. Install updated version
devtools::install()
```

## Resources

- [Rcpp Documentation](http://www.rcpp.org/)
- [RcppArmadillo](https://github.com/RcppCore/RcppArmadillo)
- [R Packages Book](https://r-pkgs.org/)
- [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html)
