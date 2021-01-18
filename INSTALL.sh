#!/bin/bash
set -e
# Get conda location information
CONDA=`which conda`
CONDA_DIRNAME=`dirname $CONDA`
MINICONDA=`dirname $CONDA_DIRNAME`
SOURCE=$MINICONDA/etc/profile.d/conda.sh

# Create each environment and install EukMetaSanity within it
conda env create -f EukMetaSanity/tasks/official/run/environment.yml
source $SOURCE
conda activate EukMS_run
python -m pip install .
conda deactivate
conda env create -f EukMetaSanity/tasks/official/report/environment.yml
source $SOURCE
conda activate EukMS_report
python -m pip install .
conda deactivate
conda env create -f EukMetaSanity/tasks/official/refine/environment.yml
source $SOURCE
conda activate EukMS_refine
python -m pip install .
conda deactivate

# Create bin directory and install non-conda dependencies
mkdir -p bin
cd bin || return
wget https://mmseqs.com/metaeuk/metaeuk-linux-sse41.tar.gz
tar "xzf" metaeuk-linux-sse41.tar.gz && rm metaeuk-linux-sse41.tar.gz
mv metaeuk/bin/metaeuk ./_metaeuk
rm -r metaeuk/ && mv _metaeuk metaeuk && cd - || return

# Add links to scripts and api to bin
ln -srf EukMetaSanity/tasks/scripts/*.py bin/
ln -srf EukMetaSanity/api/*.py bin/
