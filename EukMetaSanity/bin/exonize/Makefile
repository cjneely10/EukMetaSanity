CC = g++
CFLAGS = -g -O3 -Wall -Wvla -Wshadow -Wunreachable-code -Wconversion -Wno-sign-conversion -Wno-sign-compare -Wno-write-strings -fdiagnostics-color=always

.PHONY: all

all:
	$(CC) $(CFLAGS) gene_structure.h exonize.cpp -o exonize

