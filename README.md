# EukMetaSanity

## Installation
```
git clone https://github.com/cjneely10/EukMetaSanity.git
cd EukMetaSanity
make all
cd ..
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

Add EukMetaSanity to your PATH and PYTHONPATH variables

```
export PATH=/path/to/EukMetaSanity/bin/:$PATH
export PYTHONPATH=/path/to/EukMetaSanity/:$PYTHONPATH
ln -s /path/to/EukMetaSanity/EukMetaSanity.py ~/bin/EukMetaSanity
```

### Installing dependencies

**EukMetaSanity** consists of 4 subprograms - `run`, `report`, `fast_refine`, and `refine`.

When selecting dependencies to download, please follow the instructions below:

---
#### Run

Install [gffread](https://github.com/gpertea/gffread) AND [gffcompare](https://github.com/gpertea/gffcompare)

`awk`, `sed`, `grep`, `cp`, `rm`, should be on your PATH

##### Run step 1: Taxonomy identification
Install [MMseqs2](https://github.com/soedinglab/MMseqs2)

##### Run step 2: Repeats modeling
Optionally download [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) AND 
[RepeatMasker](http://www.repeatmasker.org/RMDownload.html)

##### Run step 3: *Ab initio* predictions
Choose from [Augustus](https://github.com/Gaius-Augustus/Augustus) or 
[GeneMark](http://topaz.gatech.edu/GeneMark/license_download.cgi).

*(Ensure that AUGUSTUS_CONFIG_PATH environmental variable is set and writable prior to running)*

##### Run step 4: Initial protein evidence
Install [MetaEuk](https://github.com/soedinglab/metaeuk)

---
#### Report

##### Report step 1/2 (optional)
Install [MMseqs2](https://github.com/soedinglab/MMseqs2), and install any databases you wish to incorporate by following
the [github page](https://github.com/soedinglab/mmseqs2/wiki#downloading-databases)

##### Report step 3 (optional)
Install [kofamscan](https://www.genome.jp/tools/kofamkoala/) and its required databases

##### Report step 5 (optional)
Install [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper) and its required databases

---
#### Fast refine

Install [sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-sort.html) AND
[BEDtools](https://github.com/arq5x/bedtools2/releases/tag/v2.29.2)

##### Fast refine step 1: Incorporate RNAseq data (optional)
Install [hisat2](https://ccb.jhu.edu/software/hisat2/manual.shtml#building-from-source)

##### Fast refine step 2: Incorporate trancriptomic data (optional)
Install [gmap](http://research-pub.gene.com/gmap/)

---
#### Refine

##### Refine step 1 (optional)
Install [MAKER3](http://www.yandell-lab.org/software/maker.html) AND 
[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)

##### Refine step 2 (optional)
Install [sambamba](https://lomereiter.github.io/sambamba/), [hisat2](https://github.com/lh3/minimap2), AND
[BRAKER2](https://github.com/Gaius-Augustus/BRAKER)

##### Refine step 3 (optional)
Install [GeMoMa](http://www.jstacs.de/index.php/GeMoMa)

## Usage

After running `download-data.py`, config files will be available in the database
directory. These can be edited to fit your needs.

If the `download-data.py` script was not used, then the default config files will be available in the 
`config` directory.

```
usage: EukMetaSanity.py [-h] -f FASTA_DIRECTORY -c CONFIG_FILE [-x EXTENSIONS]
                        [-o OUTPUT] [-d]
                        command

Run EukMetaSanity pipeline

positional arguments:
  command               Select from run/report/refine/fast_refine

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA_DIRECTORY, --fasta_directory FASTA_DIRECTORY
                        Directory of FASTA files to annotate, or paths_summary.tsv for report step
  -c CONFIG_FILE, --config_file CONFIG_FILE
                        Config file
  -x EXTENSIONS, --extensions EXTENSIONS
                        Gather files matching list of extensions separated by '/', default .fna/.fasta/.fa
  -o OUTPUT, --output OUTPUT
                        Output directory, default out
  -d, --debug           Developer mode: display all commands on single thread, default False
```

### Example usage


#### Run
Copy and edit the `run-config.ini` config file to fit your analysis needs.

For a directory of MAGs:
```
MAGs/
  |-- mag1.fna
  |-- mag2.fna
```

Generate initial annotation models using the command:

```
EukMetaSanity -f MAGs -c run-config.ini run
```

Add the `-x <ext>` if your extension does not match the default list, or if other files are present in the directory.

This will create a directory structure resembling:
```
out/
  |-- wdir/
  |-- run-eukmetasanity.log
  |-- results/
      |-- run-paths_summary.tsv
      |-- run/
          |-- mag1/
              ..
          |-- mag2/
              .. 
```

#### Fast refine (optional)
Copy and edit the `fast_refine-config.ini` config file to fit your analysis needs.

Generate additional models using the command:

```
EukMetaSanity -f out/run-paths_summary.tsv -c fast_refine-config.ini fast_refine
```

This will update the directory structure:
```
out/
  |-- wdir/
  |-- fast_refine-eukmetasanity.log
  |-- run-eukmetasanity.log
  |-- results/
      |-- fast_refine-paths_summary.tsv
      |-- run-paths_summary.tsv
      |-- fast_refine/
          |-- mag1/
              ..
          |-- mag2/
              ..
      |-- run/
          |-- mag1/
              ..
          |-- mag2/
              .. 
          ..
```

#### Report (optional)
Copy and edit the `report-config.ini` config file to fit your analysis needs.

Annotate gene models using the command:

```
EukMetaSanity -f out/{}-paths_summary.tsv -c report-config.ini report
```

Replacing `{}` with either `run` or `fast_refine` (if this step was completed).

This will update the directory structure:
```
out/
  |-- wdir/
  |-- report-eukmetasanity.log
  |-- fast_refine-eukmetasanity.log
  |-- run-eukmetasanity.log
  |-- results/
      |-- report-paths_summary.tsv
      |-- fast_refine-paths_summary.tsv
      |-- run-paths_summary.tsv
      |-- report/
          |-- mag1/
              ..
          |-- mag2/
              ..
      |-- fast_refine/
          |-- mag1/
              ..
          |-- mag2/
              ..
      |-- run/
          |-- mag1/
              ..
          |-- mag2/
              .. 
          ..
```

## Citations

(All programs used in process)
(Python documentation + module documentation)
https://www.biostars.org/p/261203/


