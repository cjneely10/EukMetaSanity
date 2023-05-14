CONDA = conda
YAPIM = yapim

BIN = bin
SRC = EukMetaSanity/src
PIPELINE_EXT = -pipeline

RUN = run
REPORT = report
REFINE = refine
DEPENDENCIES = dependencies


build: $(RUN) $(REPORT) $(REFINE)

clean: $(BIN)
	rm -rf $(BIN)/$(RUN)$(PIPELINE_EXT) $(BIN)/$(REPORT)$(PIPELINE_EXT) $(BIN)/$(REFINE)$(PIPELINE_EXT)

$(RUN): $(BIN)
	echo n | yapim create -t $(SRC)/$(RUN) -d $(SRC)/$(DEPENDENCIES) -o $(BIN)/$(RUN)$(PIPELINE_EXT)

$(REFINE): $(BIN)
	echo n | yapim create -t $(SRC)/$(REFINE) -d $(SRC)/$(DEPENDENCIES) -o $(BIN)/$(REFINE)$(PIPELINE_EXT)

$(REPORT): $(BIN)
	echo n | yapim create -t $(SRC)/$(REPORT) -d $(SRC)/$(DEPENDENCIES) -o $(BIN)/$(REPORT)$(PIPELINE_EXT)

$(BIN):
	mkdir $(BIN)

.phony: build clean
