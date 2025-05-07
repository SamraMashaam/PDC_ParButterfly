# Makefile for sequential version

CC = g++
CFLAGS = -std=c++17 -O3 -Wall
TARGET = v1
SRC = v1.cpp
LOG_DIR = logs

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $<

run: $(TARGET)
	@mkdir -p $(LOG_DIR)
	@chmod -R 777 $(LOG_DIR)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.o
	rm -rf $(LOG_DIR)