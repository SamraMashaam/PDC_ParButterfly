# Makefile for V2: MPI parallel butterfly counting

CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -Wall
NUM_PROCESSES = 8
SPARSIFICATION_P = 1.0
DATASET = ../Data/rMatGraph_J_5_100.txt
LOG_DIR = logs
TARGET = butterfly_v2

.PHONY: all run-v2 cluster-v2 clean

all: $(TARGET)

$(TARGET): v2.cpp
	$(CXX) $(CXXFLAGS) v2.cpp -o $(TARGET)

run-v2: $(TARGET)
	@mkdir -p $(LOG_DIR)
	@chmod -R 777 $(LOG_DIR)
	mpirun -np $(NUM_PROCESSES) --oversubscribe ./$(TARGET)

cluster-v2:
	./run_v2.sh $(NUM_PROCESSES) $(SPARSIFICATION_P)

clean:
	rm -f $(TARGET) *.o
	rm -rf $(LOG_DIR)