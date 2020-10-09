INSPATH = /home/kmarvos/lib
CC = g++
CPP_FLAGS = -O3 -g --std=c++17 -Wall -I$(INSPATH)
COMPILE = $(CC) $(CPP_FLAGS)
BUILD_DIR = build
BIN_DIR = bin
SRC_DIR = src
LIB_DIR = lib
TEST_DIR = test

all: $(BUILD_DIR)/.dirstamp $(BIN_DIR)/.dirstamp \
		$(BUILD_DIR)/ntest.o \
		$(BUILD_DIR)/hashtable.o \
		$(BUILD_DIR)/argparse.o
	@echo All libraries built!

tests: all $(BIN_DIR)/hash_test $(BIN_DIR)/argparse_test
	@echo All tests built!

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
	rm -f $(INSPATH)/hashtable.h $(INSPATH)/ntest.h $(INSPATH)/argparse.h

run-tests: tests
	@echo Running all tests
	$(BIN_DIR)/hash_test
	$(BIN_DIR)/argparse_test
	@echo All tests complete!

run-hash-test: $(BIN_DIR)/.dirstamp $(BIN_DIR)/hash_test
	@echo Running Hash Test!
	$(BIN_DIR)/hash_test
	@echo Hash Test complete!

run-argparse-test: $(BIN_DIR)/.dirstamp $(BIN_DIR)/argparse_test
	@echo Running ArgParse Test!
	$(BIN_DIR)/argparse_test --store-true meow vroom
	@echo ArgParse test complete!


########################################################################################################

## Directory creation

$(BUILD_DIR)/.dirstamp:
	mkdir -p $(BUILD_DIR)
	touch $@

$(BIN_DIR)/.dirstamp:
	mkdir -p $(BIN_DIR)
	touch $@

########################################################################################################

## NTest

$(BUILD_DIR)/ntest.o: NTest/$(LIB_DIR)/ntest.h
	@echo Building $@
	$(COMPILE) -c $^ -o $@
	ln -fs $^ $(INSPATH)/

########################################################################################################

## HashTable

$(BUILD_DIR)/hashtable.o: HashTable/$(LIB_DIR)/hashtable.h
	@echo Building $@
	$(COMPILE) -c $^ -o $@
	ln -fs $^ $(INSPATH)/

# Tests

$(BUILD_DIR)/tester.o: HashTable/$(LIB_DIR)/tester.h $(BUILD_DIR)/ntest.o 
	@echo Building $@
	$(COMPILE) -c $< -o $@

$(BIN_DIR)/hashtable-test.o: HashTable/$(LIB_DIR)/hashtable-test.h \
		$(BUILD_DIR)/hashtable.o $(BUILD_DIR)/ntest.o
	@echo Building $@
	$(COMPILE) -c $< -o $@

$(BIN_DIR)/hash_test: HashTable/$(SRC_DIR)/hash_test.cpp $(BUILD_DIR)/hashtable.o $(BUILD_DIR)/tester.o
	@echo Building $@
	$(COMPILE) -o $@ $<

########################################################################################################

## ArgParse

$(BUILD_DIR)/argparse.o: ArgParse/$(LIB_DIR)/argparse.h
	@echo Building $@
	$(COMPILE) -c $< -o $@
	ln -fs $^ $(INSPATH)/

# Tests

$(BUILD_DIR)/aptest.o: ArgParse/$(LIB_DIR)/test.h $(BUILD_DIR)/argparse.o
	@echo Building $@
	$(COMPILE) -c $< -o $@

$(BIN_DIR)/argparse_test: ArgParse/$(SRC_DIR)/argparse-test.cpp $(BUILD_DIR)/argparse.o \
		$(BUILD_DIR)/aptest.o
	@echo Building $@
	$(COMPILE) -o $@ $<
