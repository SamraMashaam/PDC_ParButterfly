# Makefile for sequential version
CC = g++
CFLAGS = -std=c++17 -O3 -Wall
TARGET = serial
SRC = serial.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $<

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.o

.PHONY: all run clean
