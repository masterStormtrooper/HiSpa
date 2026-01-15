# HiSpa: Hierarchical Inference of Spatial Positions from Hi-C Data

HiSpa is a Hierarchical Bayesian model for inferring 3D chromatin structures from Hi-C contact matrices.

## Dependencies

### Required
- **C++17 compatible compiler**: 
  - macOS: `clang++` (Apple Clang 12.0+)
  - Linux: `g++` (GCC 7.0+)
- **Armadillo C++ Library** (9.0+): Linear algebra and matrix operations
  ```bash
  # macOS (Homebrew)
  brew install armadillo
  
  # Ubuntu/Debian
  sudo apt-get install libarmadillo-dev
  
  # From source
  wget http://sourceforge.net/projects/arma/files/armadillo-12.0.0.tar.xz
  tar -xvf armadillo-12.0.0.tar.xz
  cd armadillo-12.0.0
  cmake .
  make
  sudo make install
  ```

- **OpenMP**: Parallel processing support
  ```bash
  # macOS (Homebrew)
  brew install libomp
  
  # Usually included with g++ on Linux
  ```

### Optional
- **spdlog** (bundled): Fast C++ logging library (included in `include/spdlog/`)
- **CLI11** (bundled): Command line parser (included as `include/CLI11.hpp`)

## Installation

### Clone the Repository
```bash
git clone https://github.com/masterStormtrooper/HiSpa.git
cd HiSpa
```

### Compile
```bash
make
```

The compiled binary will be created at `build/hispa`.

### Clean Build
```bash
make clean
make
```

## Usage

### Basic Usage
```bash
build/hispa -i <contact_matrix.txt> -o <output_directory>
```

### Example
```bash
build/hispa \
  -i data/example_contact_matrix.txt \
  -o outputs/example_run \
  --num-clusters 1 \
  --mcmc-iterations 4000 \
  --mcmc-burn-in 1000
```

### Command Line Options

#### Input/Output
- `-i, --input <file>`: Input Hi-C contact matrix file (required)
- `-o, --output <dir>`: Output directory (required)

#### MCMC Parameters
- `--mcmc-iterations <n>`: Number of MCMC iterations (default: 6000)
- `--mcmc-burn-in <n>`: Burn-in period (default: 1000)
- `--mcmc-initial-sd <val>`: Initial proposal standard deviation (default: 0.1)
- `--mcmc-sd-floor <val>`: Minimum SD for adaptive MCMC (default: 0.001)
- `--mcmc-sd-ceil <val>`: Maximum SD for adaptive MCMC (default: 0.3)

#### Clustering
- `--num-clusters <n>`: Number of clusters (default: auto-detect)
- `--cluster-distance-threshold <val>`: number of adjacent clusters that are neighbours(default: 2.0)
- `--use-cluster-init`: Use cluster-based initialization

#### Advanced Options
- `--save-samples`: Save MCMC trace samples
- `--sample-interval <n>`: Save samples every n iterations (default: 5)
- `--verbose`: Enable verbose output

#### Help
- `-h, --help`: Display help message

### Input Format

Contact matrix should be a space/tab-delimited text file with `n x n` numeric values:
```
0 10 5 2
10 0 8 3
5 8 0 6
2 3 6 0
```

## Output Files

HiSpa generates the following output files in the specified output directory:

- `final_positions.txt`: Final inferred 3D positions (n × 3 matrix)
- `best_pos_mat.txt`: Best position matrix from MCMC chain
- `initial_positions.txt`: Initial positions before MCMC
- `parameters.txt`: Final model parameters (beta0, beta1, log-likelihood)
- `hispa.log`: Detailed execution log
- `samples/mcmc_samples_*.txt`: MCMC trace samples (if `--save-samples` enabled)
- `intermediate_results/`: Checkpoints every 50 iterations

```

## License

This project is licensed under the GNU General Public License v3.0 - see below for details.

```
HiSpa - Hierarchical Inference of Spatial Positions from Hi-C Data
Copyright (C) 2026

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
```
