#!/usr/bin/env python3
import os
import setuptools
from setuptools import setup

VERSION = '0.1.0.0'

_file = os.path.abspath(os.path.dirname(__file__))

setup(
    name='EukMetaSanity',
    version=VERSION,
    description='Structural and functional annotation of eukaryotic metagenome-assembled genomes through a highly '
                'parallel and customizable data workflow​',
    long_description=open(os.path.join(_file, 'README.md'), "r").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/cjneely10/EukMetaSanity",
    author="Christopher Neely",
    author_email="christopher.neely1200@gmail.com",
    license="GNU GPL 3",
    install_requires=open(os.path.join(_file, 'requirements.txt'), "r").readlines(),
    python_requires='>=3.6',
    packages=setuptools.find_packages(),
    include_package_data=True,
    scripts=[
        "EukMetaSanity/scripts/generate-class.py",
    ]
)
