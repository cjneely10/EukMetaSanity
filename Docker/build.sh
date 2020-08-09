#!/bin/bash
cd ~/Docker/MetaSanity
git checkout v0.1.2
docker Docker -t pipedm:v0.0.6 .
docker tag pipedm:v0.0.6 cjneely10/metasanity:v0.1.2
docker push cjneely10/metasanity:v0.1.2
