#!/bin/bash
# Get conda location information
CONDA=`which conda`
CONDA_DIRNAME=`dirname $CONDA`
MINICONDA=`dirname $CONDA_DIRNAME`
SOURCE=$MINICONDA/etc/profile.d/conda.sh

# Install mamba if not already present
conda install mamba -n base -c conda-forge

# Create each environment and install EukMetaSanity within it
mamba env create -f bin/run-pipeline/run/environment.yml
source $SOURCE
conda activate EukMS_run
python -m pip install .
conda deactivate
mamba env create -f bin/report-pipeline/report/environment.yml
source $SOURCE
conda activate EukMS_report
python -m pip install .
conda deactivate
mamba env create -f bin/refine-pipeline/refine/environment.yml
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
