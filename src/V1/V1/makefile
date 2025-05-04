# Makefile for parallel version using MPI and METIS
CC = mpic++  # Use MPI compiler
CFLAGS = -std=c++17 -O3 -Wall
METIS_LIB = -lmetis  # Link METIS library
MPI_LIB = -lmpi  # Link MPI library
TARGET = v1
SRC = v1.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(MPI_LIB) $(METIS_LIB)

run: $(TARGET)
	mpirun -np 4 ./$(TARGET)  

clean:
	rm -f $(TARGET) *.o

.PHONY: all run clean
