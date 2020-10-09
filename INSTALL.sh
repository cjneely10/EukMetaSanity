#!/bin/bash
conda env create -f environment.yml
mkdir -p bin
cd bin || return
wget https://mmseqs.com/metaeuk/metaeuk-linux-sse41.tar.gz
tar "xzf" metaeuk-linux-sse41.tar.gz && rm metaeuk-linux-sse41.tar.gz
mv metaeuk/bin/metaeuk ./_metaeuk
rm -r metaeuk/ && mv _metaeuk metaeuk && cd - || return
ln -srf EukMetaSanity/scripts/*.py bin/
ln -srf EukMetaSanity/scripts/*.sh bin/
CURR_DIR=$(pwd)
cd EukMetaSanity/scripts/MergeRegions || return
LIB_DIR=$(pwd)
make argparse INSPATH="$LIB_DIR"
cd MergeRegions || return
make INSPATH="$LIB_DIR"
cd "$CURR_DIR" || return
ln -srf EukMetaSanity/scripts/MergeRegions/MergeRegions/bin/merge_regions bin/