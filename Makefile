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
SOURCES := src/hispa_main.cpp src/chromosome.cpp src/dna.cpp

# Test files
TEST_SOURCES := tests/test_dna.cpp src/dna.cpp
TEST_TARGET := $(BUILD_DIR)/test_dna
TEST_CIRCLE_SOURCES := tests/test_dna_circle.cpp src/dna.cpp
TEST_CIRCLE_TARGET := $(BUILD_DIR)/test_dna_circle
TEST_CLUSTER_INIT_SOURCES := tests/test_dna_cluster_init.cpp src/dna.cpp
TEST_CLUSTER_INIT_TARGET := $(BUILD_DIR)/test_dna_cluster_init
TEST_REAL_DNA_TARGET := $(BUILD_DIR)/test_real_dna
TEST_REAL_DNA_SOURCES := tests/test_real_dna.cpp src/dna.cpp


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

# Rule to build and run the DNA test
test: $(TEST_TARGET)
	@echo "  RUNNING DNA TEST..."
	@$(TEST_TARGET)

# Rule to compile the DNA test executable
$(TEST_TARGET): $(TEST_SOURCES)
	@# Create the build directory if it doesn't exist
	@mkdir -p $(BUILD_DIR)
	@echo "  COMPILING DNA TEST... $@"
	@$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Rule to build and run the DNA circle test
test-circle: $(TEST_CIRCLE_TARGET)
	@echo "  RUNNING DNA CIRCLE TEST..."
	@$(TEST_CIRCLE_TARGET)

# Rule to compile the DNA circle test executable
$(TEST_CIRCLE_TARGET): $(TEST_CIRCLE_SOURCES)
	@# Create the build directory if it doesn't exist
	@mkdir -p $(BUILD_DIR)
	@echo "  COMPILING DNA CIRCLE TEST... $@"
	@$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Rule to build and run the DNA cluster initialization test
test-cluster-init: $(TEST_CLUSTER_INIT_TARGET)
	@echo "  RUNNING DNA CLUSTER INITIALIZATION TEST..."
	@$(TEST_CLUSTER_INIT_TARGET)

# Rule to build and run the DNA cluster initialization test
test-real-dna: $(TEST_REAL_DNA_TARGET)
	@echo "  RUNNING DNA CLUSTER INITIALIZATION TEST..."
	@$(TEST_REAL_DNA_TARGET)

# Rule to compile the DNA cluster initialization test executable
$(TEST_CLUSTER_INIT_TARGET): $(TEST_CLUSTER_INIT_SOURCES)
	@# Create the build directory if it doesn't exist
	@mkdir -p $(BUILD_DIR)
	@echo "  COMPILING DNA CLUSTER INIT TEST... $@"
	@$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Rule to compile the DNA cluster initialization test executable
$(TEST_REAL_DNA_TARGET): $(TEST_REAL_DNA_SOURCES)
	@# Create the build directory if it doesn't exist
	@mkdir -p $(BUILD_DIR)
	@echo "  COMPILING DNA REAL DNA TEST... $@"
	@$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)
# Rule to clean up all generated files
clean:
	@echo "  CLEANING $(BUILD_DIR)"
	@rm -rf $(BUILD_DIR)

