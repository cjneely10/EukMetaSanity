[![Build Status](https://travis-ci.com/cjneely10/EukMetaSanity.svg?token=M4ut94Kepv6qucNU1mEy&branch=main)](https://travis-ci.com/cjneely10/EukMetaSanity)

# EukMetaSanity

## About
Eukaryotic genome annotation is a laborious and time-intensive process. **EukMetaSanity** provides a structural and 
functional annotation of MAGs in a highly-parallel fashion, allowing for quick and in-depth analyses. The software is
customizable - users may choose from several provided options based on their analysis needs, and power users with Python
experience can easily extend the **EukMetaSanity** code base to add to or create new pipelines!

This software suite is broken up into several sub-programs

![](assets/eukmetasanity_pipeline.png)

### Run
Identify putative taxonomy using the OrthoDB/MMETSP databases and MMseqs2 and annotate repeated regions of
the genome with either MMseqs2 or RepeatModeler/RepeatMasker. 

Generate *ab initio* structural predictions of coding regions of the genome using either Augustus or GeneMark.
Refine predictions with protein evidence using MetaEuk.

### Refine
Map RNA-seq (using HISAT2) and assembled transcriptome (using GMAP) evidence from closely related organisms (same 
organism or species) to the genome to add additional evidence using BRAKER2. 

### Report
Search KEGG, EggNOG, and any MMseqs2 database for functional annotation of putative proteins.

Check the quality of your annotation using BUSCO.

## Installation

See <a href="https://github.com/cjneely10/EukMetaSanity/blob/main/INSTALLATION.md" target="_blank">INSTALLATION.md</a> 
for detailed installation instructions.

## Usage

After running `download-data.py`, config files will be available in the `data` 
directory. These can be edited to fit your needs.

```
usage: EukMetaSanity.py [-h] -f FASTA_DIRECTORY -c CONFIG_FILE [-x EXTENSIONS]
                        [-o OUTPUT] [-d]
                        command

Run EukMetaSanity pipeline

positional arguments:
  command               Select from run/report/refine

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
```

### Example usage

#### Run
Copy and edit the `run-config.yaml` config file to fit your analysis needs.

Ensure your input FASTA sequences do not have the pipe (`|`) character present.

For a directory of MAGs:
```
MAGs/
  |-- mag1.fna
  |-- mag2.fna
```

Generate initial ab initio and protein-based annotation models using the command:

```
EukMetaSanity -f MAGs -c run-config.yaml run
```

Add `-x ext[/ext]` if your extension does not match the default list, or if other files are present in the directory.

This will create a directory structure resembling:
```
out/
  |-- wdir/
  |-- run-eukmetasanity.log
  |-- results/
    |-- run/
        |-- run.pkl
        |-- mag1/
          |-- mask.final.gff3  # Repeat regions in GFF3 format
          |-- mask.final.tbl  # Summary table of identified repeats
          |-- mag1.mask.fna  # Masked FASTA sequence
          |-- tax-report.txt  # Taxonomy report from MMseqs2
          |-- mag1.augustus.gff3  # Augustus predictions
          |-- mag1.gmes.gff3  # GeneMark predictions
          |-- mag1.all.gff3  # Combined prediction tracks
          |-- mag1.all.tiern.nr.gff3  # Tiered output predictions in GFF3 format
          |-- mag1.all.tiern.cds.fna  # Tiered output CDS predictions in FASTA format
          |-- mag1.all.tiern.faa  # Tiered output protein predictions in FASTA format
        |-- mag2/
          .. 
```

### Note on running:
**EukMetaSanity** will not re-run already completed steps within a given pipeline. If you would like to re-do a particular
portion of the pipeline, simply delete its directories in the project structure. For example, to redo the `taxonomy` step
of the `run` pipeline for all MAGs, run the following command to delete all existing data:

```
rm -r out/wdir/*/taxonomy*
```

#### Refine (optional)
Copy and edit the `refine-config.yaml` config file to fit your analysis needs. Pay close attention to the input format
for RNA-seq and transcriptomes that is required by the config file:

```
# Paths to RNA-seq should be contained in a file with the format (excluding spaces around tab):
file-basename \t /path/to/r1.fq,/path/to/r2.fq;/path/to/r3.fq,/path/to/r4.fq
# Transcriptomes should be contained in a file with the format (excluding spaces around tab):
file-basename \t /path/to/tr1.fna,/path/to/tr2.fna
``` 

The listed paired-end or single-end reads will be mapped to the file that begins with `file-basename`, as will the list 
of transcriptomes.

Integrate RNAseq and transcriptomic evidence into annotation models using the command:

```
EukMetaSanity -c refine-config.yaml refine
```

This will update the directory structure:
```
out/
  |-- wdir/
  |-- refine-eukmetasanity.log
  |-- run-eukmetasanity.log
  |-- results/
      |-- refine/
          |-- refine.pkl
          |-- mag1/
              |-- mag1.nr.gff3  # Final predictions
              |-- mag1.cds.fna  # CDS sequences
              |-- mag1.faa  # Protein sequences
              |-- augustus.hints.gtf  # Augustus predictions
              |-- genemark.gtf  # GeneMark predictions
          |-- mag2/
              ..
      |-- run/
          |-- run.pkl
          |-- mag1/
              ..
          |-- mag2/
              .. 
```

#### Report (optional)
Copy and edit the `report-config.yaml` config file to fit your analysis needs. Set the `INPUT/base` section to be either
`run` or `refine`, depending on which set of predictions you want to annotate. You may also adjust the tier you wish to 
annotate.

Annotate gene models using the command:

```
EukMetaSanity -c report-config.yaml report
```

This will update the directory structure:
```
out/
  |-- wdir/
  |-- report-eukmetasanity.log
  |-- refine-eukmetasanity.log
  |-- run-eukmetasanity.log
  |-- results/
      |-- report/
          |-- report.pkl
          |-- mag1/
              ... (results based on annotation programs run)
          |-- mag2/
              ..
      |-- refine/
          |-- refine.pkl
          |-- mag1/
              ..
          |-- mag2/
              ..
      |-- run/
          |-- run.pkl
          |-- mag1/
              ..
          |-- mag2/
              .. 
```

## Citations

Neely, C. J., & Tully, B. *EukMetaSanity*. Source code is available at [https://github.com/cjneely10/EukMetaSanity](). 
It is implemented in Python 3 under the GNU General Public License v3.0.

Also cite all dependencies that you used, as **EukMetaSanity** would not be possible were it not for the developers of 
these programs.

(citation list upcoming)
