PYTHON = python3
PIP = pip3
BASE = EukMetaSanity
IN_SCRIPTS_PATH = scripts
IN_BIN_PATH = bin
OUT_BIN_PATH = bin

.PHONY : build install all

build:
	$(PYTHON) -m compileall . > /dev/null
	mkdir -p $(OUT_BIN_PATH)
	ln -srf $(BASE)/$(IN_SCRIPTS_PATH)/*.py $(OUT_BIN_PATH)/

install:
	$(PIP) install -r requirements.txt

all: build install