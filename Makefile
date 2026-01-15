# --- Compiler and Flags ---
# Define the C++ compiler
CXX := /opt/homebrew/opt/llvm/bin/clang++

# Define compiler flags, including standards, warnings, and include paths
CXXFLAGS := -std=c++17 \
            -fcolor-diagnostics -fansi-escape-codes \
            -fopenmp -g -O3 \
            -I/opt/homebrew/include \
            -I/opt/homebrew/opt/libomp/include \
            -Iinclude \
            -DARMA_DONT_USE_WRAPPER

# Define linker flags, including library paths and the libraries to link
LDFLAGS := -L/opt/homebrew/lib \
           -L/opt/homebrew/opt/libomp/lib \
           -lfmt -lm -lgsl -lgslcblas \
           -framework Accelerate

# --- Project Structure ---
# Define the output directory and the final executable name
BUILD_DIR := build
TARGET := $(BUILD_DIR)/hispa

# List of source files with their paths
SOURCES := src/hispa_main.cpp src/chromosome.cpp


# --- Rules ---
# Phony targets are not actual files; they are just names for commands.
.PHONY: all clean test test-circle test-cluster-init test-real-dna

# The default rule, executed when you just run `make`
all: $(TARGET)

# Rule to compile and link the final executable in a single step
# It depends directly on the source files.
$(TARGET): $(SOURCES)
	@# Create the build directory if it doesn't exist
	@mkdir -p $(BUILD_DIR)
	@echo "  COMPILING & LINKING... $@"
	@$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Rule to clean up all generated files
clean:
	@echo "  CLEANING $(BUILD_DIR)"
	@rm -rf $(BUILD_DIR)

