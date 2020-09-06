# Installation

```
./install.sh
export PATH=$(pwd)/bin/:$PATH
export PYTHONPATH=$(pwd)/:$PYTHONPATH
```

Create a link to a directory on your PATH to make **EukMetaSanity** more easily callable

```
ln -s $(pwd)/EukMetaSanity.py ~/bin/EukMetaSanity
```

**EukMetaSanity**'s conda installation is packaged with all (most) of the required dependencies.

Users who wish to use GeneMark must install it separately - see the Installing dependencies section below.

Users who wish to use kofamscan or EggNOG-mapper can also follow the instructions below.

Add EukMetaSanity to your PATH and PYTHONPATH variables

### Installing dependencies

All dependencies below are automatically installed when using the conda installation.

#### Run

Install [gffread](https://github.com/gpertea/gffread) AND [gffcompare](https://github.com/gpertea/gffcompare)

`awk`, `sed`, `grep`, `cp`, `rm`, `gunzip`, `cat` should be on your PATH

##### Run step 1: Taxonomy identification
Install [MMseqs2](https://github.com/soedinglab/MMseqs2)

##### Run step 2: Repeats modeling
Optionally download [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) AND 
[RepeatMasker](http://www.repeatmasker.org/RMDownload.html)

##### Run step 3: *Ab initio* predictions
Choose from [Augustus](https://github.com/Gaius-Augustus/Augustus) or (optionally)
[GeneMark](http://topaz.gatech.edu/GeneMark/license_download.cgi)

*(Ensure that AUGUSTUS_CONFIG_PATH environmental variable is set and writable prior to running)*

##### Run step 4: Initial protein evidence
Install [MetaEuk](https://github.com/soedinglab/metaeuk)

---
#### Report

`sqlite3` should be on your PATH in addition to the values in **Run**

##### Report step 1/2 (optional)
Install [MMseqs2](https://github.com/soedinglab/MMseqs2), and install any databases you wish to incorporate by following
the [github page](https://github.com/soedinglab/mmseqs2/wiki#downloading-databases)

Generate indices and/or linear search indices for these databases to speed up search times

##### Report step 3 (optional)
Install [kofamscan](https://www.genome.jp/tools/kofamkoala/) and its required databases

##### Report step 4 (optional)
Install [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper) and its required databases

##### Report step 5 (optional)
Install [Rfam/Infernal](https://docs.rfam.org/en/latest/genome-annotation.html) and its required databases

---
#### Fast refine

Install [sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-sort.html) AND
[BEDtools](https://github.com/arq5x/bedtools2/releases/tag/v2.29.2)

##### Fast refine step 1: Incorporate RNAseq data (optional)
Install [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml#building-from-source)

##### Fast refine step 2: Incorporate trancriptomic data (optional)
Install [GMAP](http://research-pub.gene.com/gmap/)


## Installing required databases
```
EukMetaSanity/download-data.py databases -t <threads> -m <split-memory-limit>
```

The `download-data.py` script is provided to download all required base data.

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
