[![Build Status](https://travis-ci.com/cjneely10/EukMetaSanity.svg?token=M4ut94Kepv6qucNU1mEy&branch=main)](https://travis-ci.com/cjneely10/EukMetaSanity)

# EukMetaSanity

## About
Eukaryotic genome annotation is a laborious and time-intensive process. **EukMetaSanity** provides a structural and 
functional annotation of MAGs in a highly-parallel fashion, allowing for quick and in-depth analyses. The software is
customizable - users may choose from several provided options based on their analysis needs, and power users with Python
experience can easily extend the **EukMetaSanity** code base to add to or create new pipelines!

This software suite is broken up into several sub-programs

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

### Run
Copy and edit the `run-config.yaml` config file to fit your analysis needs. 

```shell
cp $EukMS_run/run-config.yaml ./
```

Activate your `EukMS_run` conda environment.

```shell
conda activate EukMS_run
```

Ensure your input FASTA sequences do not have the pipe (`|`) character present.

For a directory of MAGs:
```
MAGs/
  |-- mag1.fna
  |-- mag2.fna
```

Generate initial ab initio and protein-based annotation models using the command:

```
yapim run -i MAGs -c run-config.yaml -p $EukMS_run
```

This will create a directory structure resembling:
```
out/
  |-- wdir/
  |-- run-eukmetasanity.log
  |-- results/
    |-- run/
        |-- run.pkl
        |-- mag1/
          |-- mag1.n.Tier.gff3  # Tier n final results (GFF3)
          |-- mag1.n.Tier.faa  # Tier n final results (FASTA)
          |-- mag1.AbinitioAugustus.gff3  # Augustus results (GFF3)
          |-- mag1.AbinitioAugustus.faa  # Augustus results (FASTA)
          |-- mag1.AbinitioGeneMark.gff3  # GeneMark results (GFF3)
          |-- mag1.AbinitioGeneMark.faa  # GeneMark results (FASTA)
          |-- mag1.MetaEukEV.gff3  # MetaEuk results (GFF3)
          |-- mag1.MetaEukEV.faa  # MetaEuk results (FASTA)
          |-- mask.final.Repeats.gff3  # Repeats locations (GFF3)
          |-- mask.final.Repeats.tbl  # Summary of repeats annotation
          |-- mag1.mask.Repeats.fna  # Masked input genome (FASTA)
          |-- tax-report.Taxonomy.txt  # Taxonomy assignment summary
        |-- mag2/
          .. 
```

#### Note on running:
**EukMetaSanity** will not re-run already completed steps within a given pipeline. If you would like to re-do a particular
portion of the pipeline, simply delete its directories in the project structure. For example, to redo the `Taxonomy` step
of the `run` pipeline for all MAGs, run the following command to delete all existing data:

```
yapim clean -p $EukMS_run Taxonomy
```

### Refine (optional)
Copy and edit the `refine-config.yaml` config file to fit your analysis needs. 

```shell
cp $EukMS_refine/refine-config.yaml ./
```

Activate your `EukMS_refine` conda environment.

```shell
conda activate EukMS_refine
```

Pay close attention to the input format
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
yapim run -c refine-config.yaml -p $EukMS_refine
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

### Report (optional)
Copy and edit the `report-config.yaml` config file to fit your analysis needs. 

```shell
cp $EukMS_report/report-config.yaml ./
```

Activate your `EukMS_report` conda environment.

```shell
conda activate EukMS_report
```

Set the `INPUT/base` section to be either
`run` or `refine`, depending on which set of predictions you want to annotate. Activate your `EukMS_report` conda environment.

Annotate gene models using the command:

```
yapim run -c report-config.yaml -p $EukMS_report
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
