#!/bin/bash
conda env create -f environment.yml
mkdir -p bin
cd bin || return
wget https://mmseqs.com/metaeuk/metaeuk-linux-sse41.tar.gz
tar "xzf" metaeuk-linux-sse41.tar.gz && rm metaeuk-linux-sse41.tar.gz
mv metaeuk/bin/metaeuk ./_metaeuk
rm -r metaeuk/ && mv _metaeuk metaeuk && cd - || return
ln -srf EukMetaSanity/tasks/scripts/*.py bin/
ln -srf EukMetaSanity/api/*.py bin/
