# Compiler and flags
CC = nvc
CFLAGS = -O3 -acc -Minfo=accel \
         -I/leonardo/prod/spack/06/install/0.22/linux-rhel8-icelake/nvhpc-24.5/hdf5-1.14.3-2soj72l4cjp6mx7chrqu3h2qxzmrugnc/include \
         -I/leonardo/prod/spack/06/install/0.22/linux-rhel8-icelake/gcc-8.5.0/nvhpc-24.5-torlmnyzcexnrs6pq4cccabv7ehkv3xy/Linux_x86_64/24.5/comm_libs/12.4/hpcx/hpcx-2.19/ompi/include

LDFLAGS = -acc \
          -L/leonardo/prod/spack/06/install/0.22/linux-rhel8-icelake/nvhpc-24.5/hdf5-1.14.3-2soj72l4cjp6mx7chrqu3h2qxzmrugnc/lib \
          -L/leonardo/prod/spack/06/install/0.22/linux-rhel8-icelake/gcc-8.5.0/nvhpc-24.5-torlmnyzcexnrs6pq4cccabv7ehkv3xy/Linux_x86_64/24.5/compilers/lib \
          -L/leonardo/prod/spack/06/install/0.22/linux-rhel8-icelake/gcc-8.5.0/nvhpc-24.5-torlmnyzcexnrs6pq4cccabv7ehkv3xy/Linux_x86_64/24.5/cuda/lib64 \
          -lnvomp -lhdf5 -lm

# Directories
SRC_DIR := source
OBJ_DIR := object
BIN := vlasov_poisson

# Source and object files
SOURCES := $(wildcard $(SRC_DIR)/*.c)
OBJECTS := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SOURCES))

# Default target
all: $(BIN)

# Link the binary
$(BIN): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Compile .c to .o into object/
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Create object directory if it doesn't exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Clean all
clean:
	rm -f $(BIN) $(OBJ_DIR)/*.o *.csv *.bin *.log

# Debug build
debug: CFLAGS = -g -O0 -acc -Minfo=accel
debug: clean all

# Install dependencies (example)
install-deps:
	sudo apt-get update
	sudo apt-get install nvidia-hpc-sdk

# Test
test: $(BIN)
	./$(BIN) -n 256 -T 1.0 -freq 5 -tree

# Performance test
perf: $(BIN)
	./$(BIN) -n 4096 -T 5.0 -freq 10 -tree

# Help
help:
	@echo "Available targets:"
	@echo "  all       - Build the simulation (default)"
	@echo "  clean     - Remove build files and output"
	@echo "  test      - Run small test simulation"
	@echo "  perf      - Run performance test"
	@echo "  debug     - Build debug version"
	@echo "  help      - Show this help"

.PHONY: all clean debug test perf help install-deps

