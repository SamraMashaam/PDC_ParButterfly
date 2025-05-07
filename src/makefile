# Top-level Makefile for PDC_ParButterfly

.PHONY: all v1 v2 v3 run-v1 run-v2 run-v3 run-v2-cluster run-v3-cluster clean help

all: v1 v2 v3

v1:
	@echo "Building sequential version..."
	$(MAKE) -C V1 $(TARGET)

v2:
	@echo "Building MPI parallel version..."
	$(MAKE) -C V2 $(TARGET)

v3:
	@echo "Building OpenMP parallel version..."
	$(MAKE) -C V3 $(TARGET)

run-v1:
	@echo "Running sequential version..."
	$(MAKE) -C V1 run

run-v2:
	@echo "Running MPI parallel version locally..."
	$(MAKE) -C V2 run-v2

run-v3:
	@echo "Running OpenMP parallel version..."
	$(MAKE) -C V3 run-v3

run-v2-cluster:
	@echo "Running MPI parallel version on cluster..."
	$(MAKE) -C V2 cluster-v2

run-v3-cluster:
	@echo "Running OpenMP parallel version on cluster..."
	$(MAKE) -C V3 cluster-v3

clean:
	@echo "Cleaning all versions..."
	$(MAKE) -C V1 clean
	@if [ -d V2 ]; then $(MAKE) -C V2 clean; fi
	@if [ -d V3 ]; then $(MAKE) -C V3 clean; fi

help:
	@echo "Available targets:"
	@echo "  all          - Build all versions (v1, v2, v3)"
	@echo "  v1           - Build sequential version"
	@echo "  v2           - Build MPI version"
	@echo "  v3           - Build OpenMP version"
	@echo "  run-v1       - Run sequential version"
	@echo "  run-v2       - Run MPI version locally"
	@echo "  run-v3       - Run OpenMP version"
	@echo "  run-v2-cluster - Run MPI version on cluster"
	@echo "  run-v3-cluster - Run OpenMP version on cluster"
	@echo "  clean        - Remove all build artifacts"
	@echo "  help         - Show this help message"