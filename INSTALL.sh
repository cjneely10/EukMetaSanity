#!/bin/bash
conda env create -f environment.yml
conda activate EukMS
pip install -r requirements.txt
python -m compileall . > /dev/null
mkdir -p bin
cd bin || return
wget https://mmseqs.com/metaeuk/metaeuk-linux-sse41.tar.gz
tar "xzf" metaeuk-linux-sse41.tar.gz && rm metaeuk-linux-sse41.tar.gz
mv metaeuk/bin/metaeuk ./_metaeuk
rm -r metaeuk/ && mv _metaeuk metaeuk
ln -srf EukMetaSanity/scripts/*.py bin/
ln -srf EukMetaSanity/scripts/*.sh bin/