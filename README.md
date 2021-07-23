# EukMetaSanity

## About
Eukaryotic genome annotation is a laborious and time-intensive process. **EukMetaSanity** provides a structural and 
functional annotation of MAGs in a highly-parallel fashion, allowing for quick and in-depth analyses.

This software suite is broken up into several sub-programs

### Run
Perform high-quality structural annotation using protein evidence recruited from the Orthologous Database of Proteins and the Marine 
Microbial Eukaryotic Transcriptomic Sequencing Project databases.

Performs:

- Taxonomy prediction using proteins identified by MetaEuk
- DNA repeat modeling and masking using RepeatModeler and RepeatMasker, respectively
- *Ab initio* prediction with Augustus and/or GeneMark

### Refine
Map RNA-seq (using HISAT2) and assembled transcriptome (using GMAP) evidence from closely related organisms (same 
organism or species) to the genome to add additional evidence using BRAKER2. 

### Report
Search KEGG, EggNOG, and any MMseqs2 database for functional annotation of putative proteins.

Check the quality of your annotation using BUSCO.

## Installation

See <a href="https://github.com/cjneely10/EukMetaSanity/blob/main/INSTALLATION.md" target="_blank">INSTALLATION.md</a> 
for detailed installation instructions.

---

We have provided a test set of data for use in validating installation, or as a means of better understanding the EukMetaSanity implementation.

These files are present in the directory `tests/data`:

```
tests/
  |-- data/
    |-- NAO-all-DCM-20-180-00_bin-1.fna
    |-- NAO-all-DCM-20-180-00_bin-2.fna
    |-- NAO-all-DCM-20-180-00_bin-19.fna
```

If you are using your own input set, ensure that your FASTA files' extensions are one of the following:

```
.fna
.fa
.fasta
.faa
.fas
```

## Run

### Config setup
Copy and edit the `run-config.yaml` config file to fit your analysis needs.

```shell
cp $EukMS_run/run-config.yaml ./
```

In the `GLOBAL` section, you will want to set the number of `MaxThreads` you will use to run the analysis, as well as the `MaxMemory` to be used.
In the `SLURM` section, set `USE_CLUSTER` to `true` if running on slurm, and provide run configuration details (such as qos, job_name, partition, account, etc.).

In each subsequent section, you may adjust the `threads`, `memory`, and `FLAGS` that are passed to the program. Be sure to set the
time allocation for each step if running pipeline on `SLURM`.

If you choose to omit running `AbinitioAugustus` or `AbinitioGeneMark`, set its respective `skip` flag to be `true`.

```yaml
---  # document start

###########################################
## Pipeline input section
INPUT:
  root: all

## Global settings
GLOBAL:
  # Maximum threads/cpus to use in analysis
  MaxThreads: 20
  # Maximum memory to use (in GB)
  MaxMemory: 120

###########################################

SLURM:
  ## Set to True if using SLURM
  USE_CLUSTER: false
  ## Pass any flags you wish below
  ## DO NOT PASS the following:
  ## --nodes, --ntasks, --mem, --cpus-per-task
  --qos: unlim
  --job-name: EukMS
  user-id: uid

MetaEukEV:
  # Number of threads task will use
  threads: 16
  # Amount of memory task will use (in GB)
  memory: 100
  time: "12:00:00"
  dependencies:
    MetaEuk:
      program: metaeuk
      data:
        /path/to/odb-mmetsp_db
      # Pass any flags to metaeuk required
      FLAGS:
        --min-length 30
        --metaeuk-eval 0.0001
        -s 5
        --cov-mode 0
        -c 0.3
        -e 100
        --max-overlap 0
        --remove-tmp-files

Taxonomy:
  # Number of threads task will use
  threads: 16
  # Amount of memory task will use (in GB)
  memory: 100
  time: "8:00:00"
  dependencies:
    MMSeqsCreateDB:
      program: mmseqs
      threads: 1

    MMSeqsTaxonomy:
      program: mmseqs
      cutoff: 8.0  # Minimum % of mapped reads to tax level
      data:
        /path/to/odb-mmetsp_db
      # Pass any flags to mmseqs required
      FLAGS:
        --remove-tmp-files
        -s 6.5
        --min-seq-id 0.40
        -c 0.3
        --cov-mode 0

Repeats:
  # Number of threads task will use
  threads: 16
  # Amount of memory task will use (in GB)
  memory: 16
  time: "24:00:00"
  dependencies:
    RModBuildDatabase:
      time: "10:00"
      threads: 1
      program: BuildDatabase

    RModRepeatModeler:
      program: RepeatModeler
      skip: false

    RMaskRepeatMasker:
      program: RepeatMasker
      level: family
      data:
        "" # Comma-separated list of repeat models to incorporate
      FLAGS:
        -nolow

    RMaskProcessRepeats:
      time: "30:00"
      threads: 1
      program: ProcessRepeats
      FLAGS:
        -nolow

    RMaskRMOut:
      time: "10:00"
      threads: 1
      program: rmOutToGFF3.pl

AbinitioGeneMark:
  # Number of threads task will use
  threads: 16
  # Amount of memory task will use (in GB)
  memory: 16
  time: "16:00:00"
  skip: false
  dependencies:
    MMSeqsFilterTaxSeqDB:
      threads: 1
      time: "10:00"
      program: mmseqs
      level: order
      data:
        /path/to/odb-mmetsp_db

    GeneMarkProtHint:
      program: prothint.py

    GeneMarkPETAP:
      program: gmes_petap.pl
      FLAGS:
        --min_contig 100
        --max_contig 1000000000
        --max_gap 5000
        --max_mask 5000
        --min_contig_in_predict 100
        --min_gene_in_predict 10
        --gc_donor 0.001
        --max_intron 10000
        --max_intergenic 50000
        --soft_mask auto

AbinitioAugustus:
  # Number of threads task will use
  threads: 16
  # Amount of memory task will use (in GB)
  memory: 8
  time: "24:00:00"
  skip: false
  dependencies:
    MMSeqsCreateDB:
      time: "10:00"
      program: mmseqs
      threads: 1

    MMSeqsSearch:
      memory: 20
      data:
        /path/to/odb-mmetsp_db
      program: mmseqs
      subname: linsearch  # Can use `search` if linear index is not present for database
      FLAGS:
        --cov-mode 0
        -c 0.6
        -e 0.01
        --remove-tmp-files

    MMSeqsConvertAlis:
      time: "2:00:00"
      data:
        /path/to/odb-mmetsp_db
      program: mmseqs
      FLAGS:
        --format-output query,target,pident,taxid,taxname,taxlineage

    Augustus:
      program: augustus
      cutoff: 25.0
      rounds: 1

Tier:
  # Number of threads task will use
  threads: 1
  # Amount of memory task will use (in GB)
  memory: 2
  time: "10:00"
  program: locus_solver
  tier: 1

...  # document end
```

### Call pipeline

Activate your `EukMS_run` conda environment.

```shell
conda activate EukMS_run
```

Ensure your input FASTA sequences do not have the pipe (`|`) character present.

Run the pipeline using the command:

```
yapim run -i /path/to/EukMetaSanity/tests/data -c run-config.yaml -p $EukMS_run
```

If run with default config parameters, the analysis should complete within about 4 hours. This can be sped up by running the analysis on `SLURM` and/or increasing the thread count.

### Run output

This will create a directory structure resembling:
```
out/
  |-- wdir/
  |-- input/
  |-- run-eukmetasanity.log
  |-- results/
    |-- run/
        |-- run.pkl
        |-- NAO-all-DCM-20-180-00_bin-1
          |-- NAO-all-DCM-20-180-00_bin-1.n.Tier.gff3  # Tier n final results (GFF3)
          |-- NAO-all-DCM-20-180-00_bin-1.n.Tier.faa  # Tier n final results (FASTA)
          |-- NAO-all-DCM-20-180-00_bin-1.AbinitioAugustus.gff3  # Augustus results (GFF3)
          |-- NAO-all-DCM-20-180-00_bin-1.AbinitioAugustus.faa  # Augustus results (FASTA)
          |-- NAO-all-DCM-20-180-00_bin-1.AbinitioGeneMark.gff3  # GeneMark results (GFF3)
          |-- NAO-all-DCM-20-180-00_bin-1.AbinitioGeneMark.faa  # GeneMark results (FASTA)
          |-- NAO-all-DCM-20-180-00_bin-1.MetaEukEV.gff3  # MetaEuk results (GFF3)
          |-- NAO-all-DCM-20-180-00_bin-1.MetaEukEV.faa  # MetaEuk results (FASTA)
          |-- NAO-all-DCM-20-180-00_bin-1.Repeats.gff3  # Repeats locations (GFF3)
          |-- NAO-all-DCM-20-180-00_bin-1.Repeats.tbl  # Summary of repeats annotation
          |-- NAO-all-DCM-20-180-00_bin-1.Repeats.fna  # Masked input genome (FASTA)
          |-- NAO-all-DCM-20-180-00_bin-1.Taxonomy.txt  # Taxonomy assignment summary
        |-- NAO-all-DCM-20-180-00_bin-19/
          ...
        |-- NAO-all-DCM-20-180-00_bin-2/
          ...
```

#### Note on running:
**EukMetaSanity** will not re-run already completed steps within a given pipeline. If you would like to re-do a particular
portion of the pipeline, simply delete its directories in the project structure. For example, to redo the `Taxonomy` step
of the `run` pipeline for all MAGs, run the following command to delete all existing data:

```
yapim clean -p $EukMS_run Taxonomy -o out
```

The preceding command will not only delete the results generated by the `Taxonomy` step, but will also remove all other steps
in the pipeline that use these results. If you only wish to delete a step and nothing else, run:

```shell
rm -r out/wdir/*/Taxonomy*
```

---

## Refine (optional)

For the `Refine` pipeline, we do not provide transcriptomic data for testing. If you wish to test this installation, you must provide your own set of test data.

### Config setup

Copy and edit the `refine-config.yaml` config file to fit your analysis needs.

```shell
cp $EukMS_refine/refine-config.yaml ./
```

As before, set resource usage in the `GLOBAL` settings as well as for each subsequent section. Also set the `SLURM` settings, if needed.

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

```yaml
---  # document start

###########################################
## Pipeline input section
INPUT:
  run: all

## Global settings
GLOBAL:
  # Maximum threads/cpus to use in analysis
  MaxThreads: 24
  # Maximum memory to use (in GB)
  MaxMemory: 100

###########################################

SLURM:
  ## Set to True if using SLURM
  USE_CLUSTER: false
  ## Pass any flags you wish below
  ## DO NOT PASS the following:
  ## --nodes, --ntasks, --mem, --cpus-per-task
  --qos: unlim
  --job-name: EukMS
  user-id: uid

CollectInput:
  # Number of threads task will use
  threads: 1
  # Amount of memory task will use (in GB)
  memory: 8
  time: "4:00:00"
  # Should be in format (excluding spaces around tab):
  # file-basename \t /path/to/tr1.fna[,/path/to/tr2.fna]
  transcriptomes: /path/to/transcriptome-mapping-file
  # Should be in format (excluding spaces around tab):
  # file-basename \t /path/to/r1.fq[,/path/to/r2.fq][;/path/to/r3.fq[,/path/to/r4.fq]]
  rnaseq: /path/to/rnaseq-mapping-file

GatherProteins:
  # Number of threads task will use
  threads: 8
  # Amount of memory task will use (in GB)
  memory: 8
  time: "4:00:00"
  dependencies:
    MMSeqsFilterTaxSeqDB:
      program: mmseqs
      level: order
      data:
        /path/to/odb-mmetsp_db

Transcriptomes:
  # Number of threads task will use
  threads: 8
  # Amount of memory task will use (in GB)
  memory: 8
  time: "4:00:00"
  dependencies:
    GMAPBuild:
      threads: 1
      program: gmapindex

    GMAP:
      program: gmap
      FLAGS:
        -B 5
        --input-buffer-size 1000000
        --output-buffer-size 1000000
        -f samse

RNASeq:
  # Number of threads task will use
  threads: 8
  # Amount of memory task will use (in GB)
  memory: 8
  time: "4:00:00"
  dependencies:
    Hisat2Build:
      threads: 1
      program: hisat2-build

    Hisat2:
      program: hisat2
      FLAGS:
        ""

MergeSams:
  # Number of threads task will use
  threads: 1
  # Amount of memory task will use (in GB)
  memory: 8
  time: "4:00:00"

ProcessMapping:
  # Number of threads task will use
  threads: 8
  # Amount of memory task will use (in GB)
  memory: 60
  time: "4:00:00"
  dependencies:
    SambambaView:
      program: sambamba

    SambambaSort:
      program: sambamba

RunBraker:
  # Number of threads task will use
  threads: 4
  # Amount of memory task will use (in GB)
  memory: 100
  time: "4:00:00"
  dependencies:
    Braker:
      program: braker.pl
      FLAGS:
      # Provide flags as desired

...  # document end
```

### Call pipeline

Activate your `EukMS_refine` conda environment.

```shell
conda activate EukMS_refine
```

Run pipeline with the command:

```
yapim run -i /path/to/EukMetaSanity/tests/data -c refine-config.yaml -p $EukMS_refine
```

### Refine output

This will update the directory structure:

```
out/
  |-- input/
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
              ...
      |-- run/
          ...
```

## Report (optional)

### Config setup

Copy and edit the `report-config.yaml` config file to fit your analysis needs.

```shell
cp $EukMS_report/report-config.yaml ./
```

As before, set resource usage in the `GLOBAL` settings as well as for each subsequent section. Also set the `SLURM` settings, if needed.

Set `skip` to `true` if you wish to skip any of the steps listed in the pipeline. If running `EggNog`, ensure to provide the path to the eggnog databases on your system.

If you wish to use the results of the `refine` pipeline instead of default results from the `run` pipeline, update the `INPUT` section as follows:

```yaml
INPUT:
  root: all
  refine: all
```

Otherwise, set the protein output you wish to annotate:

```yaml
---  # document start

###########################################
## Pipeline input section
INPUT:
  root: all
  run:
    prot: merged-prot  # or genemark-prot or aug-prot or evidence-prot

## Global settings
GLOBAL:
  # Maximum threads/cpus to use in analysis
  MaxThreads: 20
  # Maximum memory to use (in GB)
  MaxMemory: 100

###########################################

SLURM:
  ## Set to True if using SLURM
  USE_CLUSTER: false
  ## Pass any flags you wish below
  ## DO NOT PASS the following:
  ## --nodes, --ntasks, --mem, --cpus-per-task
  --qos: unlim
  --job-name: EukMS
  user-id: uid

Quality:
  # Number of threads task will use
  threads: 16
  # Amount of memory task will use (in GB)
  memory: 8
  time: "4:00:00"
  skip: false
  dependencies:
    Busco:
      program: busco
      mode: prot
      lineage: eukaryota
      FLAGS:
        -f

MMSeqs:
  # Number of threads task will use
  threads: 16
  # Amount of memory task will use (in GB)
  memory: 90
  time: "4:00:00"
  skip: false
  dependencies:
    MMSeqsCreateDB:
      program: mmseqs

    MMSeqsSearch:
      data:
        /path/to/odb-mmetsp_db
      program: mmseqs
      subname: linsearch
      FLAGS:
        -c 0.3
        --cov-mode 1
        --remove-tmp-files

    MMSeqsConvertAlis:
      data:
        /path/to/odb-mmetsp_db
      program: mmseqs

KOFamScan:
  # Number of threads task will use
  threads: 1
  # Amount of memory task will use (in GB)
  memory: 8
  time: "4:00:00"
  skip: true
  dependencies:
    KofamscanExecAnnotation:
      program: /path/to/kofamscan/exec_annotation
      kolist: /path/to/kofam/ko_list
      profiles: /path/to/profiles/eukaryote.hal

EggNog:
  # Number of threads task will use
  threads: 1
  # Amount of memory task will use (in GB)
  memory: 8
  time: "4:00:00"
  skip: true
  dependencies:
    EMapper:
      program: emapper.py
      FLAGS:
        # Provide flags as needed
        -m diamond
        --override
        --data_dir /location/of/eggnog-data

...  # document end
```

### Call pipeline

Activate your `EukMS_report` conda environment.

```shell
conda activate EukMS_report
```

Run pipeline using the command:

```
yapim run -c report-config.yaml -p $EukMS_report
```

Note that we do not need to provide the input directory for this analysis, as the pipeline will only annotate genomes that have completed the `Run` or `Refine` pipeline.

With default settings, the analysis should complete in less than 10 minutes.

### Report output

This will update the directory structure:
```
out/
  |-- wdir/
  |-- input/
  |-- report-eukmetasanity.log
  |-- refine-eukmetasanity.log
  |-- run-eukmetasanity.log
  |-- results/
      |-- report/
          |-- report.pkl
          |-- NAO-all-DCM-20-180-00_bin-1/
              ... (results based on annotation programs run)
          |-- NAO-all-DCM-20-180-00_bin-19/
              ...
          |-- NAO-all-DCM-20-180-00_bin-2/
              ...
      |-- refine/
          ...
      |-- run/
          ...
```

## Citations

"The high-throughput gene prediction of more than 1,700 eukaryote genomes using the software package EukMetaSanity" by Neely, Hu, Alexander, and Tully, 2021.

Also cite all dependencies that you used, as **EukMetaSanity** would not be possible were it not for the developers of 
these programs.

(citation list upcoming)
