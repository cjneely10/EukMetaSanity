INSPATH = /home/kmarvos/lib
CC = g++
CPP_FLAGS = -O3 -g --std=c++17 -Wall -I$(INSPATH)
COMPILE = $(CC) $(CPP_FLAGS)
BUILD_DIR = build
BIN_DIR = bin
SRC_DIR = src
LIB_DIR = lib
TEST_DIR = test

all: pre $(BUILD_DIR)/.dirstamp $(BIN_DIR)/.dirstamp

pre:
	cd ../ && make && cd -

run: all run-parser-test run-gene-test

run-parser-test: all $(BIN_DIR)/parser-test
	@echo Running Parser test
	$(BIN_DIR)/parser-test /home/kmarvos/Data/TestData/MergeRegions/hyp-prots.faa /home/kmarvos/Data/TestData/MergeRegions/spo.gff3

run-gene-test: all $(BIN_DIR)/gene-test
	@echo Running Genome test
	$(BIN_DIR)/gene-test

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

## Directory creation

$(BUILD_DIR)/.dirstamp:
	mkdir -p $(BUILD_DIR)
	touch $@

$(BIN_DIR)/.dirstamp:
	mkdir -p $(BIN_DIR)
	touch $@

## Parser test

$(BIN_DIR)/parser-test: $(SRC_DIR)/parser-test.cpp $(LIB_DIR)/parser-test.h
	@echo Building $@
	$(COMPILE) -o $@ $<  ../$(BUILD_DIR)/parser.o ../$(BUILD_DIR)/gene.o

$(BIN_DIR)/gene-test: $(SRC_DIR)/gene-test.cpp $(LIB_DIR)/gene-test.h
	@echo Building $@
	$(COMPILE) -o $@ $<  ../$(BUILD_DIR)/gene.o

