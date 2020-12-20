#!/usr/bin/env python

"""
EukMetaSanity.py
===============================================
EukMetaSanity - Generate structural/functional annotations for Eukaryotes

"""

import os
import json
import logging
from pathlib import Path
from signal import signal, SIGPIPE, SIG_DFL
from typing import Generator, List, Tuple, Dict
from Bio import SeqIO
from EukMetaSanity.utils.arg_parse import ArgParse
from EukMetaSanity.tasks.base.path_manager import PathManager
from EukMetaSanity.tasks.base.task_manager import TaskManager
from EukMetaSanity.tasks.base.config_manager import ConfigManager
from EukMetaSanity.tasks.manager.pipeline_manager import PipelineManager


# Logging initialize
def _initialize_logging(ap: ArgParse):
    """

    :param ap:
    :return:
    """
    # Initialize logging
    log_file = os.path.join(ap.args.output, "%s-eukmetasanity.log" % ap.args.command)
    if os.path.exists(log_file):
        os.remove(log_file)
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO, filename=log_file, filemode='w')
    print(
        "*" * 80, "",
        "All log statements are redirected to %s" % log_file, "",
        "*" * 80, "",
        "Displaying step summaries here:\n",
        sep="\n"
    )


# Gather all files to parse that match user-passed extensions
def _files_iter(ap: ArgParse, storage_dir: str) -> Generator[str, ArgParse, None]:
    for file in os.listdir(ap.args.fasta_directory):
        for ext in ap.args.extensions:
            if file.endswith(ext):
                yield _simplify_fasta(ap, file, storage_dir)


def _simplify_fasta(ap: ArgParse, file, storage_dir: str) -> str:
    # Simplify FASTA of complex-named sequences
    fasta_file = str(Path(os.path.join(ap.args.fasta_directory, file)).resolve())
    out_file = os.path.join(storage_dir, os.path.basename(os.path.splitext(fasta_file)[0]) + ".fna")
    if os.path.exists(out_file):
        return out_file
    SeqIO.write(SeqIO.parse(fasta_file, "fasta"), out_file, "fasta")
    return out_file


# Get program-needed list of files for this step in pipeline
def _get_list_of_files(summary_file: str) -> Tuple[List[str], List[Dict[str, Dict[str, object]]]]:
    data = json.load(open(summary_file, "r"))
    out_ids = sorted(list(data.keys()))
    out_dict_list = []
    for _id in out_ids:
        to_add = {"root": {}}
        for key, val in data[_id].items():
            if isinstance(val, dict):
                to_add["root"][key] = val
            else:
                to_add[key] = val
        out_dict_list.append(to_add)
    return out_ids, out_dict_list


# # Driver logic
def _main(ap: ArgParse, cfg: ConfigManager, is_continued: bool, tpm: PipelineManager):
    # Generate primary path manager
    pm = PathManager(ap.args.output)
    # Begin logging
    _initialize_logging(ap)
    # Gather list of files to analyze
    # TODO: Handle loading from wdir all required input for specific pipeline
    if is_continued:
        # Gather from existing data
        for f in (logging.info, print):
            f("Getting files from last run...")
        input_prefixes, input_files = _get_list_of_files(ap.args.fasta_directory)
    else:
        # Simplify FASTA files into working directory
        pm.add_dirs("MAGS")
        input_files = list(_file for _file in _files_iter(ap, pm.get_dir("MAGS")))
        # List of prefixes for tracking each file's progress
        input_prefixes = [os.path.basename(os.path.splitext(_file)[0]) for _file in input_files]
        input_files = [{"root": {"fna": input_file}} for input_file in input_files]
        # Create base dir for each file to analyze
        all([pm.add_dirs(_file) for _file in input_prefixes])

    # # Begin task list
    tm = TaskManager(tpm, cfg, pm, input_files, input_prefixes, ap.args.debug, ap.args.command)
    tm.run(ap.args.output)


# Parse user arguments
def _parse_args(ap: ArgParse, tm: PipelineManager) -> Tuple[ConfigManager, bool]:
    # Confirm path existence
    assert os.path.exists(ap.args.config_file)
    is_continued = False
    # A single file was provided - use this as input to run
    if ap.args.fasta_directory is None:
        is_continued = True
    # Ensure command is valid
    assert ap.args.command in tm.programs.keys()
    if ap.args.debug is True:
        ap.args.debug = 0
    else:
        ap.args.debug = 1
    # Determine file extensions to keep
    ap.args.extensions = ap.args.extensions.split("/")
    return ConfigManager(ap.args.config_file), is_continued


if __name__ == "__main__":
    # Redirect dask nanny errors
    signal(SIGPIPE, SIG_DFL)
    DEFAULT_EXTS = ".fna/.fasta/.fa"
    _tm = PipelineManager()
    _ap = ArgParse(
        (
            (("command",),
             {"help": "Select from %s" % "/".join(list(_tm.programs.keys()))}),
            (("-f", "--fasta_directory",),
             {"help": "Directory of FASTA files to annotate"}),
            (("-c", "--config_file"),
             {"help": "Config file", "required": True}),
            (("-x", "--extensions"),
             {"help": "Gather files matching list of extensions separated by '/', default %s" % DEFAULT_EXTS,
              "default": DEFAULT_EXTS}),
            (("-o", "--output"),
             {"help": "Output directory, default out", "default": "out"}),
            (("-d", "--debug"),
             {"help": "Developer mode: display all commands on single thread, default False", "default": False,
              "action": "store_true"}),
        ),
        description="Run EukMetaSanity pipeline"
    )
    _cfg, _is_cont = _parse_args(_ap, _tm)
    _main(_ap, _cfg, _is_cont, _tm)
    for func in (logging.info, print):
        func("\nEukMetaSanity %s pipeline complete!" % _ap.args.command)
