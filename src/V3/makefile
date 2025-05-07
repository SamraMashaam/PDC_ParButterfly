# Makefile for V3: MPI+OpenMP parallel butterfly counting

CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -Wall -fopenmp
DEBUGFLAGS = -std=c++17 -g -Wall -fopenmp -fsanitize=address
TARGET = butterfly_v3
LOG_DIR = logs
TEST_INPUT = test_input.txt

.PHONY: all run-v3 test debug cluster-v3 clean

all: $(TARGET)

$(TARGET): v3.cpp
	$(CXX) $(CXXFLAGS) v3.cpp -o $(TARGET)

run-v3: $(TARGET)
	@mkdir -p $(LOG_DIR)
	@chmod -R 777 $(LOG_DIR)
	./$(TARGET)

test: $(TARGET)
	@mkdir -p $(LOG_DIR)
	@chmod -R 777 $(LOG_DIR)
	@echo -e "3\n2\n3\n6\n0\n1\n" > $(TEST_INPUT)
	@./$(TARGET) < $(TEST_INPUT)
	@rm -f $(TEST_INPUT)
	@cat $(LOG_DIR)/output_*.txt

debug: v3.cpp
	$(CXX) $(DEBUGFLAGS) v3.cpp -o $(TARGET)
	@mkdir -p $(LOG_DIR)
	@chmod -R 777 $(LOG_DIR)
	@echo -e "3\n2\n3\n6\n0\n1\n" > $(TEST_INPUT)
	@./$(TARGET) < $(TEST_INPUT)
	@rm -f $(TEST_INPUT)
	@cat $(LOG_DIR)/output_*.txt

cluster-v3:
	./run_v3.sh

clean:
	rm -f $(TARGET) $(TEST_INPUT) *.o
	rm -rf $(LOG_DIR)