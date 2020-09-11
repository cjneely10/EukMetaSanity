# Installation

`awk`, `sed`, `grep`, `cp`, `rm`, `gunzip`, `cat`, `conda` should be on your PATH.

Ensure you have `conda`&ge;4.8.3 installed, that you have conda activated, and that you are in your `(base)` conda 
environment. Then, run the following commands:

```
git clone https://github.com/cjneely10/EukMetaSanity.git
cd EukMetaSanity && ./INSTALL.sh
conda activate EukMS
```

Update your `PATH` and `PYTHONPATH` variables in your .bashrc file. Create a link to a directory on your PATH to make 
**EukMetaSanity** more easily callable:

```
echo export PATH="$(pwd)"/bin/:'$PATH' >> ~/.bashrc
echo export PYTHONPATH="$(pwd)"/:'$PYTHONPATH' >> ~/.bashrc
ln -s $(pwd)/EukMetaSanity.py ~/bin/EukMetaSanity
```

You may need to restart your shell for these changes to take effect.

**EukMetaSanity**'s conda installation is packaged with all (most) of the required dependencies.
Users who wish to use [GeneMark](http://topaz.gatech.edu/GeneMark/license_download.cgi), 
[eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper), or [kofamscan](https://www.genome.jp/tools/kofamkoala/) 
must install them separately.

*(Ensure that AUGUSTUS_CONFIG_PATH environmental variable is set and writable prior to running if you are using this ab 
initio predictor)*

*(Ensure that .gm_key is present in your home directory if you are using GeneMark as your ab initio predictor)*

### Installing RepeatMasker libraries and scripts
**RepeatMasker** can incorporate additional DFam updates. [Install these](http://www.repeatmasker.org/RMDownload.html)
in your conda installation directory. Make sure your `EukMS` conda environment is still active 
prior to updating.

```
cd /path/to/miniconda/envs/EukMS/share/RepeatMasker/Libraries/
wget https://www.dfam.org/releases/Dfam_3.2/families/Dfam.h5.gz
gunzip Dfam.h5.gz && cd ..
perl ./configure
```

The configure script should ask you to confirm the location of your installation, as well as to select your search 
engine. Select 2 for `RMBlast`, and provide the path as `/path/to/miniconda3/envs/EukMS/bin/` when requested 
(substituting for the proper path on your system).

Finally, due to a small bug in the conda `RepeatMasker` conda environment, run the following command:

```
cp util/rmOutToGFF3.pl ./
```

### Installing required databases

**The `download-data.py` script** is provided to download all other required base data. Run the script to download the 
required databases:

```
cd /path/to/EukMetaSanity
./download-data.py -t <threads> -m <max-mem> data
```

This will download the OrthoDB and RFAM databases for use in **EukMetaSanity**. Additionally, config files will 
automatically generate for use when running **EukMetaSanity**. You can speed up `mmseqs` searches by generating 
pre-computed indices with the `-x` flag. This can take a lot of space to store.


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

*(Ensure that AUGUSTUS_CONFIG_PATH environmental variable is set and writable prior to running if you are using this ab 
initio predictor)*

*(Ensure that .gm_key is present in your home directory if you are using GeneMark as your ab initio predictor)*

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
#### Refine

Install [sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-sort.html) AND
[BEDtools](https://github.com/arq5x/bedtools2/releases/tag/v2.29.2)

##### Refine step 1: Incorporate RNAseq data (optional)
Install [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml#building-from-source)

##### Refine step 2: Incorporate trancriptomic data (optional)
Install [GMAP](http://research-pub.gene.com/gmap/)

##### Refine step 3: Run BRAKER2
Install [BRAKER](https://github.com/Gaius-Augustus/BRAKER#installation)



