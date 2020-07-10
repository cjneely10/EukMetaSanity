# EukMetaSanity

## Installation
```
git clone https://github.com/cjneely10/EukMetaSanity.git
cd EukMetaSanity
make all
export PATH=/path/to/EukMetaSanity/bin/:$PATH
export PYTHONPATH=/path/to/EukMetaSanity/:$PYTHONPATH
```

### Database downloads

## Usage
```
usage: EukMetaSanity.py [-h] -f FASTA_DIRECTORY -c CONFIG_FILE [-x EXTENSIONS]
                        [-o OUTPUT] [-d]
                        command

Run EukMetaSanity pipeline

positional arguments:
  command               Select from run/refine

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA_DIRECTORY, --fasta_directory FASTA_DIRECTORY
                        Directory of FASTA files to annotate
  -c CONFIG_FILE, --config_file CONFIG_FILE
                        Config file
  -x EXTENSIONS, --extensions EXTENSIONS
                        Gather files matching '/'-separated list of extensions, default .fna/.fasta/.fa
  -o OUTPUT, --output OUTPUT
                        Output directory, default out
  -d, --debug           Developer mode: display all commands on single thread, default False
```

## Citations

(All programs used in process)
(Python documentation + module documentation)
https://www.biostars.org/p/261203/
