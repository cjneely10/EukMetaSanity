PYTHON = python
PIP = pip
BASE = EukMetaSanity
IN_SCRIPTS_PATH = scripts
IN_BIN_PATH = bin
OUT_BIN_PATH = bin

.PHONY : build install all

build:
	$(PYTHON) -m compileall . > /dev/null
	mkdir -p $(OUT_BIN_PATH)
	ln -srf $(BASE)/$(IN_SCRIPTS_PATH)/*.py $(OUT_BIN_PATH)/
	ln -srf $(BASE)/$(IN_SCRIPTS_PATH)/*.sh $(OUT_BIN_PATH)/

install:
	$(PIP) install -r requirements.txt

clean:
	rm -r $(OUT_BIN_PATH)
	pyclean .

all: build install