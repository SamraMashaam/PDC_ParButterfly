# Top-level Makefile for PDC_ParButterfly

.PHONY: all v1 v2 v3 run-v1 run-v2 cluster-v2 clean help

all: v1 v2

v1:
	@echo "Building sequential version..."
	$(MAKE) -C V1

v2:
	@echo "Building MPI parallel version locally..."
	$(MAKE) -C V2 run-v2

v3:
	@echo "V3 (OpenMP) is a placeholder, not yet implemented."

run-v1:
	@echo "Running sequential version..."
	$(MAKE) -C V1 run

run-v2:
	@echo "Running MPI parallel version locally..."
	$(MAKE) -C V2 run-v2

cluster-v2:
	@echo "Running MPI parallel version on cluster..."
	$(MAKE) -C V2 cluster-v2

clean:
	@echo "Cleaning all versions..."
	$(MAKE) -C V1 clean
	@if [ -d V2 ]; then $(MAKE) -C V2 clean; fi
	@if [ -d V3 ]; then $(MAKE) -C V3 clean; fi

help:
	@echo "Available targets:"
	@echo "  all       - Build all versions (v1, v2)"
	@echo "  v1        - Build sequential version"
	@echo "  v2        - Build and run MPI version locally"
	@echo "  v3        - Placeholder for OpenMP version"
	@echo "  run-v1    - Run sequential version"
	@echo "  run-v2    - Run MPI version locally"
	@echo "  cluster-v2 - Run MPI version on cluster"
	@echo "  clean     - Remove all build artifacts"
	@echo "  help      - Show this help message"