#!/bin/bash
cd ~/BioProjects/EukMetaSanity/Docker/ || return
docker build -t eukmetasanity:v0.1.0 .
docker tag eukmetasanity:v0.1.0 cjneely10/eukmetasanity:v0.1.0
docker push cjneely10/eukmetasanity:v0.1.0
