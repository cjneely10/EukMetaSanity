#!/bin/bash
git clone https://github.com/cjneely10/EukMetaSanity.git
cd EukMetaSanity || return
conda env create -f environment.yml
conda activate EukMS
make all
cd bin || return
wget https://mmseqs.com/metaeuk/metaeuk-linux-sse41.tar.gz
tar "xzf" metaeuk-linux-sse41.tar.gz && rm metaeuk-linux-sse41.tar.gz
mv metaeuk/bin/metaeuk ./
rm -r metaeuk/