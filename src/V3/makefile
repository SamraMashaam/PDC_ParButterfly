# Makefile for V3: MPI+OpenMP parallel butterfly counting

CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -Wall -fopenmp
DEBUGFLAGS = -std=c++17 -g -Wall -fopenmp -fsanitize=address
TARGET = butterfly_v3
LOG_DIR = logs
TEST_INPUT = test_input.txt
NUM_PROCS = 2  # Default number of MPI processes

.PHONY: all run-v3 test debug cluster-v3 clean

all: $(TARGET)

$(TARGET): v3.cpp
	$(CXX) $(CXXFLAGS) v3.cpp -o $(TARGET)

run-v3: $(TARGET)
	@mkdir -p $(LOG_DIR)
	@chmod -R 755 $(LOG_DIR)
	mpirun -np $(NUM_PROCS) ./$(TARGET)

test: $(TARGET)
	@mkdir -p $(LOG_DIR)
	@chmod -R 755 $(LOG_DIR)
	@echo -e "3\n2\n3\n3\n0\n1\n" > $(TEST_INPUT)  # uSize=3, vSize=2, numEdges=3, no sparsification, exact counting
	mpirun -np $(NUM_PROCS) ./$(TARGET) < $(TEST_INPUT)
	@rm -f $(TEST_INPUT)
	@cat $(LOG_DIR)/output_*.txt

debug: v3.cpp
	$(CXX) $(DEBUGFLAGS) v3.cpp -o $(TARGET)
	@mkdir -p $(LOG_DIR)
	@chmod -R 755 $(LOG_DIR)
	@echo -e "3\n2\n3\n3\n0\n1\n" > $(TEST_INPUT)
	mpirun -np $(NUM_PROCS) ./$(TARGET) < $(TEST_INPUT)
	@rm -f $(TEST_INPUT)
	@cat $(LOG_DIR)/output_*.txt

cluster-v3:
	@echo "Submitting cluster job (modify run_v3.sh for your cluster)"
	./run_v3.sh

clean:
	rm -f $(TARGET) $(TEST_INPUT) *.o
	rm -rf $(LOG_DIR)