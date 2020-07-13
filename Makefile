PYTHON = python3
PIP = pip3
IN_SCRIPTS_PATH = EukMetaSanity/scripts
OUT_SCRIPTS_PATH = bin

.PHONY : build install all

build:
	$(PYTHON) -m compileall . > /dev/null
	mkdir -p $(OUT_SCRIPTS_PATH)
	ln -srf $(IN_SCRIPTS_PATH)/*.py $(OUT_SCRIPTS_PATH)/

install:
	$(PIP) install -r requirements.txt

all: build install