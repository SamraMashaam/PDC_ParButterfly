# Makefile for sequential version with METIS
CC = g++
CFLAGS = -std=c++17 -O3 -Wall
LDFLAGS = -lmetis
TARGET = v1
SRC = v1.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.o

.PHONY: all run clean
