#!/usr/bin/env python3
"""
Module downloads requisite data for official pipelines in EukMetaSanity
"""
import os
import shutil
from typing import List
from pathlib import Path
from plumbum import local
from plumbum.machines.local import LocalCommand
from EukMetaSanity.utils.arg_parse import ArgParse
from EukMetaSanity.data.data import data_urls, instructions

# Dependencies
cp = local["cp"]
sed = local["sed"]
tar = local["tar"]
wget = local["wget"]
gunzip = local["gunzip"]
mmseqs = local["mmseqs"]


def _print_and_run(cmd: LocalCommand):
    """ Display bash command and run

    :param cmd: Plumbum command to run
    :raises: plumbum.commands.processes.ProcessExecutionError for failed commands
    """
    print(cmd)
    cmd()


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
    # Download each URL to folder
    out = []
    ids = []
    for _id, url in data_urls().items():
        # Download
        new_dir = os.path.join(out_dir, _id)
        if os.path.exists(new_dir):
            shutil.rmtree(new_dir)
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
        _file = os.path.join(new_dir, os.path.basename(url.url))
        _out = os.path.join(new_dir, _id)
        out.append(_out)
        ids.append(_id)
        if not os.path.exists(_out):
            _print_and_run(wget[url.url, "-O", _file])
            # Tar/gunzip
            if url.tar:
                _print_and_run(tar[url.flags, _file])
                _file = os.path.splitext(_file.replace(".tar", ""))[0]
            elif url.gz:
                _print_and_run(gunzip[_file])
                _file = os.path.splitext(_file)[0]
            # Generate MMseqs2 database
            if url.type == "FASTA":
                # Generate database
                _print_and_run(
                    mmseqs[
                        "createdb",
                        _file,
                        _out
                    ]
                )
                # Remove downloaded file
                # Add taxonomy info from ncbi
                _create_tax_db(_out, out_dir)
            elif url.type == "profile":
                _msa_db = _out[:-3] + "_msa_db"
                _print_and_run(
                    mmseqs["convertmsa", _file, _msa_db]
                )
                _print_and_run(
                    mmseqs["msa2profile", _msa_db, _out, "--match-mode", "1"]
                )
            os.remove(_file)
            if ap.args.build and url.type != "profile":
                # Generate linear index
                _print_and_run(
                    mmseqs[
                        "createlinindex",
                        _out,
                        "tmp",
                        "--threads", str(ap.args.threads),
                        "--split-memory-limit", ap.args.max_mem,
                        "--remove-tmp-files",
                    ]
                )
            # Create MMseqs2 index file
            if ap.args.index and url.type != "profile":
                _print_and_run(
                    mmseqs[
                        "createindex",
                        _out,
                        "tmp",
                        "--threads", str(ap.args.threads),
                        "--split-memory-limit", ap.args.max_mem,
                        "--remove-tmp-files",
                    ]
                )
    # Handle merge instructions
    for new_db, _instructions in instructions():
        if _instructions is not None:
            if "merge" in _instructions:
                _, merge_db1, merge_db2 = _instructions.split("|")
                new_dir = os.path.join(out_dir, merge_db1)
                if os.path.exists(new_dir):
                    shutil.rmtree(new_dir)
                if not os.path.exists(new_dir):
                    os.makedirs(new_dir)
                _out = os.path.join(new_dir, merge_db1)
                _handle_merge(_out, os.path.join(new_dir, merge_db2), os.path.join(new_dir, new_db))
                out = [_o if _o != new_dir else _out for _o in out]

    _generate_config_files(out, ids, ap.args.path)


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


def _create_tax_db(_out: str, out_dir: str):
    """ Create taxonomy database from NCBI current taxdump

    :param _out: Expected database (e.g. orthodb, etc.)
    :param out_dir: Directory in which to generate taxonomy database
    :raises: plumbum.commands.processes.ProcessExecutionError
    """
    # Create mmseqs input file
    if "ortho" in _out:
        _odb_tax_parse(_out, os.path.join(out_dir, "mmseqs.input"))
    else:
        return
    # Download tax info
    _print_and_run(
        wget["ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", "-O", os.path.join(out_dir, "taxdump.tar.gz")]
    )
    _print_and_run(tar["xzvf", os.path.join(out_dir, "taxdump.tar.gz"), "-C", out_dir])
    # Generate tax db
    _print_and_run(
        mmseqs[
            "createtaxdb",
            _out,
            "tmp",
            "--ncbi-tax-dump", out_dir,
            "--tax-mapping-file", os.path.join(out_dir, "mmseqs.input")
        ]
    )


def _handle_merge(db_name: str, merge_db_name: str, output_name: str):
    """ Merge db_name and merge_db_name into output_name

    :param db_name: First db
    :param merge_db_name: Second db
    :param output_name: New mmseqs db name to create from first two databases
    """
    _print_and_run(
        mmseqs[
            "mergedbs",
            db_name,
            output_name,
            merge_db_name
        ]
    )


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
