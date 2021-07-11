CONDA = conda
YAPIM = yapim
BIN = bin
SRC = EukMetaSanity/src
RUN = run
CONFIG_EXT = -config.yaml
PIPELINE_EXT = -pipeline
REPORT = report
REFINE = refine
DEPENDENCIES = dependencies


build: $(RUN) $(REPORT) $(REFINE)

clean: $(BIN)
	rm -rf $(BIN)/$(RUN)$(PIPELINE_EXT) $(BIN)/$(REPORT)$(PIPELINE_EXT) $(BIN)/$(REFINE)$(PIPELINE_EXT)

$(RUN): $(BIN)
	echo n | yapim create -d $(SRC)/$(DEPENDENCIES) -o $(BIN)/$(RUN)$(PIPELINE_EXT) -t $(SRC)/$(RUN)
	ln $(BIN)/$(RUN)$(PIPELINE_EXT)/$(RUN)/$(RUN)$(CONFIG_EXT) $(BIN)
	rm $(BIN)/$(RUN)$(PIPELINE_EXT)/$(RUN)$(CONFIG_EXT)

$(REFINE): $(BIN)
	echo n | yapim create -d $(SRC)/$(DEPENDENCIES) -o $(BIN)/$(REFINE)$(PIPELINE_EXT) -t $(SRC)/$(REFINE)
	ln $(BIN)/$(REFINE)$(PIPELINE_EXT)/$(REFINE)/$(REFINE)$(CONFIG_EXT) $(BIN)
	rm $(BIN)/$(REFINE)$(PIPELINE_EXT)/$(REFINE)$(CONFIG_EXT)

$(REPORT): $(BIN)
	echo n | yapim create -d $(SRC)/$(DEPENDENCIES) -o $(BIN)/$(REPORT)$(PIPELINE_EXT) -t $(SRC)/$(REPORT)
	ln $(BIN)/$(REPORT)$(PIPELINE_EXT)/$(REPORT)/$(REPORT)$(CONFIG_EXT) $(BIN)
	rm $(BIN)/$(REPORT)$(PIPELINE_EXT)/$(REPORT)$(CONFIG_EXT)

$(BIN):
	mkdir $(BIN)

.phony: build clean
