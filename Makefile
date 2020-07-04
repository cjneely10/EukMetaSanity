PYTHON = python3
PIP = pip3

.PHONY : build

build:
	$(PYTHON) -m compileall .
	$(PIP) install -r requirements.txt