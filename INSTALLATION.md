# Installation

`sed`, `grep`, `conda`, and `git` should be available on your system `PATH`.

Ensure that the drive in which your conda installation is contained has ~80GB of extra space.

Clone this repository to a drive with >100GB of storage space:

```
git clone https://github.com/cjneely10/EukMetaSanity.git
cd EukMetaSanity
```

Ensure that you have `conda`&ge;4.9 installed, that you have conda activated, and that you are in your `(base)` conda environment.


A typical installation can be run using:

```shell
./INSTALL.sh -t <num-threads>
```

Expect ~4 hours to complete the installation and database downloads.

After running the `INSTALL.sh` script, you must restart your shell.

### Updating from existing installation

Prior to updating to the most recent version of this software, users should edit their `~/.bashrc` file and remove old
EukMetaSanity-related export and `PATH` update statements. Additionally, remove old database and `bin` directories.
Next, restart your shell.

Then, either update this repository:

```shell
git restore .
git pull
```

Or clone it (as described above), and then run the installation script:
```shell
./INSTALL.sh -t <num-threads> --upgrade
```

## Installing dependencies

**EukMetaSanity**'s `conda` installation is packaged with all (most) of the required dependencies.
Users who wish to use [GeneMark](http://topaz.gatech.edu/GeneMark/license_download.cgi), 
[eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper), or [kofamscan](https://www.genome.jp/tools/kofamkoala/) 
must install them separately. We highly suggest using of these software suites (especially GeneMark), 
but they are not directly required.

EggNOG users should download the software and its required databases using `pip` with their `EukMS_report` environment loaded.

### Configuring GeneMark &ge;4.65_lic

EukMetaSanity is packaged with all dependencies that are needed to run GeneMark
If you choose to include GeneMark in your analysis pipeline, follow the installation instructions [on their webpage](http://topaz.gatech.edu/GeneMark/license_download.cgi) to download their software and accept their license agreements.

Ensure that your `.gm_key` file is present in your home directory. 
You also may need to run their accessory script from within your `EukMS_run` environment:

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

```shell
./INSTALL.sh --uninstall
```

Optionally, remove the installer library

```shell
conda remove mamba -y
```

You will also need to remove the 7 lines added to your `.bashrc` file. You may also wish to delete this repository.
