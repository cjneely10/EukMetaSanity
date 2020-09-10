# Installation

### Installing dependencies

`awk`, `sed`, `grep`, `cp`, `rm`, `gunzip`, `cat`, `conda` should be on your PATH.

Ensure you have `conda`&ge;4.8.3 installed, that you have conda activated, and that you are in your `(base)` conda 
environment. Then, run the following commands:

```
git clone https://github.com/cjneely10/EukMetaSanity.git
cd EukMetaSanity && ./INSTALL.sh
conda activate EukMS
pip install -r requirements.txt
python -m compileall . > /dev/null
```

Update your `PATH` and `PYTHONPATH` variables in your .bashrc file.

```
echo export PATH=$(pwd)/bin/:'$PATH' >> ~/.bashrc
echo export PYTHONPATH=$(pwd)/:'$PYTHONPATH' >> ~/.bashrc
```

Create a link to a directory on your PATH to make **EukMetaSanity** more easily callable

```
ln -s $(pwd)/EukMetaSanity.py ~/bin/EukMetaSanity
```

**EukMetaSanity**'s conda installation is packaged with all (most) of the required dependencies.
Users who wish to use [GeneMark](http://topaz.gatech.edu/GeneMark/license_download.cgi), 
[eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper), or [kofamscan](https://www.genome.jp/tools/kofamkoala/) 
must install them separately.

### Installing required databases
The `download-data.py` script is provided to download all required base data. Run the script to download the 
required databases:

```
download-data.py -t <threads> -m <max-mem> data
```

This will download the OrthoDB and RFAM databases for use in **EukMetaSanity**. Additionally, config files will 
automatically generate for use when running **EukMetaSanity**.

General usage:

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

**Your installation is complete!**
 
If you wish to download additional databases to use in the `Report` step, use the 
`mmseqs database` command to pull them prior to running **EukMetaSanity**, and add their location to your 
`report-config.ini` file in the `[mmseqs] DATA` section.

## Non-Conda installation

If you do not wish to use Conda, the list of all required dependencies is available below.

#### Run

Install [gffread](https://github.com/gpertea/gffread) AND [gffcompare](https://github.com/gpertea/gffcompare)

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



