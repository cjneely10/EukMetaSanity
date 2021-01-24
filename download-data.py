#!/usr/bin/env python3
"""
Module downloads requisite data for official pipelines in EukMetaSanity
"""
import os
from typing import List
from pathlib import Path
from plumbum import cli
from EukMetaSanity.data.download_utils import download_data, manage_downloaded_data


def _parse_args(ap: ArgParse) -> ArgParse:
    """ Confirm valid types for command-line arguments

    :param ap: Passed arguments
    :raises: AssertionError if attempt to index without building database
    :raises: ValueError for improper type of threads or negative number
    :return: Reference to modified ArgParse object
    """
    if ap.args.index:
        assert ap.args.build, "Build must be set in order to index!"
    try:
        ap.args.threads = int(ap.args.threads)
    except ValueError as e:
        print(e)
        exit(1)
    if ap.args.threads < 1:
        raise ValueError
    return ap


# TODO: Entire function needs modularity to pass review
def run(ap: ArgParse, out_dir: str):
    """ Download all data with wget from the data.py script - Part of API

    :param ap: Parsed arguments
    :param out_dir: Directory in which to write all downloaded data
    :raises: plumbum.commands.processes.ProcessExecutionError
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)


def _generate_config_files(_file_names: List[str], _replace_strings: List[str], _outdir: str):
    """ Copy base config files to install directory and add in proper data path locations

    :param _file_names: List of files to modify
    :param _replace_strings: String replacements
    :param _outdir: Output directory to generate output
    :raises: plumbum.commands.processes.ProcessExecutionError
    """
    _config_directory = os.path.join(os.path.dirname(__file__), "EukMetaSanity/config")
    for _config_file in os.listdir(_config_directory):
        _new_file = os.path.join(_outdir, os.path.basename(_config_file))
        cp[os.path.join(_config_directory, _config_file), _new_file]()
        for _file_name, _replace_string in zip(_file_names, _replace_strings):
            _file_name = str(Path(_file_name).resolve()).replace("/", "\/")
            _print_and_run(
                sed[
                    "-i", "s/\/path\/to\/%s/%s/g" % (_replace_string, _file_name),
                    _new_file
                ]
            )


def _odb_tax_parse(mmseqs_db_path: str, outfile: str):
    """ Generate NCBI taxonomy mappings from generated .lookup file

    :param mmseqs_db_path: Path to generated mmseqs database
    :param outfile: Output path for parsed orthodb mapping file
    """
    output_p = open(outfile, "w")
    mmseqs_input_fp = open(mmseqs_db_path + ".lookup", "r")
    for line in mmseqs_input_fp:
        line = line.split()
        output_p.write(line[1] + "\t" + line[1].split("_")[0] + "\n")
    output_p.close()


if __name__ == "__main__":
    _ap = ArgParse(
        (
            (("path",),
             {"help": "Download path"}),
            (("-b", "--build"),
             {"help": "Generate required MMseqs2 databases and linear indices, default True", "default": True,
              "action": "store_false"}),
            (("-x", "--index"),
             {"help": "Generate search index (recommended, but takes a lot of space), default False",
              "default": False, "action": "store_true"}),
            (("-t", "--threads"),
             {"help": "Number of threads to use in database generation, default 1", "default": "1"}),
            (("-m", "--max_mem"),
             {"help": "Set max memory per split. E.g. 800B, 5K, 10M, 1G; default 8G", "default": "8G"}),
        ),
        description="Download required data and build MMseqs2 databases"
    )
    run(_parse_args(_ap), _ap.args.path)
