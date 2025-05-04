# Top-level Makefile to manage all versions
.PHONY: all v1 v2 v3 run-v1 run-v2 run-v3 clean help

all: v1

# Build targets
v1:
	@echo "Building sequential version..."
	$(MAKE) -C V1

v2:
	@echo "V2 (OpenMP) is a placeholder, not yet implemented."

v3:
	@echo "V3 (OpenCL) is a placeholder, not yet implemented."

# Run targets (build + execute)
run-v1: v1
	@echo "Running sequential version..."
	$(MAKE) -C V1 run

run-v2: v2
	@echo "V2 (OpenMP) is a placeholder, not yet implemented."

run-v3: v3
	@echo "V3 (OpenCL) is a placeholder, not yet implemented."

# Clean all
clean:
	@echo "Cleaning all versions..."
	$(MAKE) -C V1 clean
	@if [ -d V2 ]; then $(MAKE) -C V2 clean; fi
	@if [ -d V3 ]; then $(MAKE) -C V3 clean; fi

# Help message
help:
	@echo "Available targets:"
	@echo "  all       - Build all versions (currently only v1)"
	@echo "  v1        - Build sequential version"
	@echo "  v2        - Placeholder for OpenMP version"
	@echo "  v3        - Placeholder for OpenCL version"
	@echo "  run-v1    - Build and run sequential version"
	@echo "  run-v2    - Placeholder for running OpenMP version"
	@echo "  run-v3    - Placeholder for running OpenCL version"
	@echo "  clean     - Remove all build artifacts"
	@echo "  help      - Show this help message"