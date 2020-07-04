PYTHON = python3
PIP = pip3

.PHONY : build install all

build:
	$(PYTHON) -m compileall .

install:
	$(PIP) install -r requirements.txt

all:
	build
	install