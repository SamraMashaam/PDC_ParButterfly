# Makefile for parallel version using MPI and METIS
CC = mpic++  # Use MPI compiler
CFLAGS = -std=c++17 -O3 -Wall -fopenmp 
METIS_LIB = -lmetis  # Link METIS library
MPI_LIB = -lmpi  # Link MPI library
OMP_FLAGS = -fopenmp
LDFLAGS = $(OMP_FLAGS) $(MPI_LIB) $(METIS_LIB)
TARGET = parallel
SRC = parallel.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o $@ $< $(LDFLAGS)

run: $(TARGET)
	mpirun -np 4 ./$(TARGET)  

clean:
	rm -f $(TARGET) *.o

.PHONY: all run clean
