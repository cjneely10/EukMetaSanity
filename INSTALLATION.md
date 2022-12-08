# Installation

`sed`, `grep`, `conda`, `wget`, and `git` should be installed on your system.

Clone this repository:

```
git clone https://github.com/cjneely10/EukMetaSanity.git
cd EukMetaSanity
```

Ensure that you have `conda`&ge;4.9.2 installed, that you have conda activated, and that you are in your `(base)` conda environment.


A typical installation can be run using:

```shell
./INSTALL.sh -t <num-threads>
```

To specify a separate directory for storing database files, provide a path using the `-d` flag:

```shell
./INSTALL.sh -t <num-threads> -d /path/to/download-location
```

Be sure you have >128GB of storage space and 4-8 hours to complete the installation and database downloads.

Your `~/.bashrc` file will be modified to append updated environment variables. You may change this using the `-b` flag.
More information is available with the `-h` flag.

After running the `INSTALL.sh` script, you must restart your shell for these changes to take effect.

## Installing dependencies

**EukMetaSanity**'s `conda` installation is packaged with all (most) of the required dependencies.
Users who wish to use [GeneMark](http://topaz.gatech.edu/GeneMark/license_download.cgi), 
[eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper), or [kofamscan](https://www.genome.jp/tools/kofamkoala/) 
must install them separately. We highly suggest using of these software suites, but they are not directly required.

EggNOG users should download the software and its required databases using `pip` with their `EukMS_report` environment loaded.

### Configuring GeneMark 4.65_lic

If you choose to include GeneMark in your analysis pipeline, follow the installation instructions [on their webpage](http://topaz.gatech.edu/GeneMark/license_download.cgi) to download their software and accept their license agreements.

Ensure that your `.gm_key` file is present in your home directory if you are using GeneMark as your ab initio predictor. 
You also may need to run their accessory script from within your `EukMS` environment:

```
perl change_path_in_perl_scripts.pl "/usr/bin/env perl"
```

Ensure that your `gmes.cfg` file has parameters that are sufficient for your dataset (min contig, etc.).

The `gmes_linux_64` directory and its enclosed `ProtHint` directory should both be on your system path:

```shell
~$ echo $PATH | tr ":" "\n"
...
/path/to/gmes_linux_64
/path/to/gmes_linux_64/ProtHint/bin
/path/to/gmes_linux_64/ProtHint/dependencies
/path/to/gmes_linux_64/ProtHint
...
```

## **Your installation is complete!**
 
If you wish to download additional databases to use in the `report` step, use the 
`mmseqs databases` command to download and install them prior to running **EukMetaSanity**, and add their location to your 
`report-config.yaml` file in the `mmseqs.dependencies` section under `MMSeqs.MMSeqsSearch.data` and 
`MMSeqs.Convertalis.data`.


## Uninstalling EukMetaSanity

```
for f in run report refine; do
    conda remove --name EukMS_$f --all -y
done
conda remove mamba -y
```

You will also need to remove the 4 lines added to your `.bashrc` file. You may also wish to delete this repository and the database directory.
