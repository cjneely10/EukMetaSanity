INSPATH = ..
APPATH = $(INSPATH)/MergeRegions/lib
CC = g++
CPP_FLAGS = -O3 -g --std=c++17 -Wall -I$(APPATH)
COMPILE = $(CC) $(CPP_FLAGS)
BUILD_DIR = build
BIN_DIR = bin
SRC_DIR = src
LIB_DIR = lib

all: $(BUILD_DIR)/.dirstamp $(BIN_DIR)/.dirstamp \
		$(BIN_DIR)/merge_regions
	@echo All libraries built!

run-merge-regions: all
	$(BIN_DIR)/merge_regions -h

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

########################################################################################################

## Directory creation

$(BUILD_DIR)/.dirstamp:
	mkdir -p $(BUILD_DIR)
	touch $@

$(BIN_DIR)/.dirstamp:
	mkdir -p $(BIN_DIR)
	touch $@

########################################################################################################

## MergeRegions

$(BIN_DIR)/merge_regions: $(SRC_DIR)/main.cpp $(BUILD_DIR)/gene.o $(BUILD_DIR)/parser.o
	@echo Building $@
	$(COMPILE) -I$(INSPATH) -o $@ $^

$(BUILD_DIR)/parser.o: $(SRC_DIR)/parser.cpp $(LIB_DIR)/parser.h $(BUILD_DIR)/gene.o
	@echo Building $@
	$(COMPILE) -c $< -o $@

$(BUILD_DIR)/gene.o: $(SRC_DIR)/gene.cpp $(LIB_DIR)/gene.h $(LIB_DIR)/region.h $(LIB_DIR)/record.h
	@echo Building $@
	$(COMPILE) -c $< -o $@
