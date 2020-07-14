# EukMetaSanity

## Installation
```
git clone https://github.com/cjneely10/EukMetaSanity.git
cd EukMetaSanity
make all
cd ..
./download-data.py databases -t <threads> -m <split-memory-limit>
```

The `download-data.py` script is provided to download all required base data. This will also generate MMseqs2
databases for the data, which can be voluminous.

```
usage: download-data.py [-h] [-b BUILD] [-x INDEX] [-o OUTPUT] [-r REWRITE]
                        [-t THREADS] [-m MAX_MEM]
                        path

Download required data and build MMseqs2 databases

positional arguments:
  path                  Download path

optional arguments:
  -h, --help            show this help message and exit
  -b BUILD, --build BUILD
                        Generate required MMseqs2 databases and linear indices, default True
  -x INDEX, --index INDEX
                        Generate search index (recommended, but takes a lot of space), default False
  -o OUTPUT, --output OUTPUT
                        Output default config files with included download paths, default True
  -r REWRITE, --rewrite REWRITE
                        Rewrite existing directory, default False
  -t THREADS, --threads THREADS
                        Number of threads to use in database generation, default 1
  -m MAX_MEM, --max_mem MAX_MEM
                        Split memory limit for database generation, default 8G
```

Add EukMetaSanity to your PATH and PYTHONPATH variables

```
export PATH=/path/to/EukMetaSanity/bin/:$PATH
export PYTHONPATH=/path/to/EukMetaSanity/:$PYTHONPATH
```

### Installing dependencies

**EukMetaSanity** consists of 3 subprograms - `run`, `refine`, `report`.

When selecting dependencies to download, please follow the instructions below:

#### Run utilities:
Install [gffread](https://github.com/gpertea/gffread)

##### Run step 1: Taxonomy identification
Install [MMseqs2](https://github.com/soedinglab/MMseqs2)

##### Run step 2: Repeats modeling
Optionally download [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) AND 
[RepeatMasker](http://www.repeatmasker.org/RMDownload.html)

##### Run step 3: *Ab initio* predictions
Choose from [Augustus](https://github.com/Gaius-Augustus/Augustus) or 
[GeneMark](http://topaz.gatech.edu/GeneMark/license_download.cgi).

*(Ensure that AUGUSTUS_CONFIG_PATH environmental variable is set prior to running)*

##### Run step 4: Initial protein evidence
Install [MetaEuk](https://github.com/soedinglab/metaeuk)

---

##### Refine step 1 (optional)
Install [MAKER3](http://www.yandell-lab.org/software/maker.html) and 
[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)

##### Refine step 2 (optional)
Install [sambamba](https://lomereiter.github.io/sambamba/), [minimap2](https://github.com/lh3/minimap2), and
[BRAKER2](https://github.com/Gaius-Augustus/BRAKER)

##### Refine step 3 (optional)
Install [GeMoMa](http://www.jstacs.de/index.php/GeMoMa)

---

##### Report step 1/2 (optional)
Install [HMMER](http://hmmer.org/)

##### Report step 3 (optional)
Install [kofamscan](ftp://ftp.genome.jp/pub/tools/kofam_scan/)

##### Report step 4 (optional)
Install [HHsuite3](https://github.com/soedinglab)

##### Report step 5 (optional)
Install [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper)

## Usage

After running `download-data.py`, config files for Run/Refine/Report will be available in the database
directory. These can be edited to fit your needs.

If the `download-data.py` script was not used, then the default config files will be available in the 
`config` directory.

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
