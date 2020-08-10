# EukMetaSanity

## About
Eukaryotic genome annotation is a laborious and time-intensive process. **EukMetaSanity** provides a structural and 
functional annotation of MAGs in a highly-parallel fashion, allowing for quick and in-depth analyses. The software is
customizable - users may choose from several provided options based on their analysis needs, and power users with Python
experience can easily extend the **EukMetaSanity** code base to add to or create new pipelines!

This software suite is broken up into several sub-programs

### Run
Identify putative taxonomy using the OrthoDB and annotate repeated regions of
the genome.

Generate *ab initio* and protein-driven structural predictions of coding regions of the genome.

### Fast refine
Map RNA-seq and assembled transcriptome evidence from closely related organisms (same organism or species) to the genome
to add additional evidence.

### Report
Identify RNA (noncoding, tRNA, etc.) regions of the genome.

Search KEGG, EggNOG, and any MMseqs2 database for functional annotation of putative proteins.

## Installation
See [INSTALLATION.md](https://github.com/cjneely10/EukMetaSanity/blob/master/INSTALLATION.md) for detailed installation
instructions

## Usage

#### Note for Docker users
Use the Docker CL API to call EukMetaSanity

After running `download-data.py`, config files will be available in the database
directory. These can be edited to fit your needs.

If the `download-data.py` script was not used, then the default config files will be available in this repo's 
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

Generate initial ab initio and protein-based annotation models using the command:

```
EukMetaSanity -f MAGs -c run-config.ini run
```

Add `-x <ext[,ext]>` if your extension does not match the default list, or if other files are present in the directory.

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

Integrate RNAseq and transcriptomic evidence into annotation models using the command:

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

Neely, C. J., & Tully, B. *EukMetaSanity*. Source code is available at [https://github.com/cjneely10/EukMetaSanity](). 
It is implemented in Python 3 under the GNU General Public License v3.0.

Also cite all dependencies that you used, as **EukMetaSanity** would not be possible were it not for the developers of 
these programs.


