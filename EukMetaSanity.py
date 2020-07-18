#!/usr/bin/env python3
import os
import sys
import logging
from Bio import SeqIO
from pathlib import Path
from plumbum import local
from string import punctuation
from typing import Generator, List, Tuple
from signal import signal, SIGPIPE, SIG_DFL
from EukMetaSanity.utils.arg_parse import ArgParse
from EukMetaSanity.utils.path_manager import PathManager
from EukMetaSanity.utils.helpers import prefix as _prefix
from EukMetaSanity.utils.config_manager import ConfigManager
from EukMetaSanity.tasks.manager.task_manager import TaskManager

"""
EukMetaSanity - Generate structural/functional annotations for simple Eukaryotes

"""


# # Helper functions
# Logging initialize
def _initialize_logging(ap: ArgParse) -> None:
    # Initialize logging
    log_file = os.path.join(ap.args.output, "eukmetasanity.log")
    if os.path.exists(log_file):
        os.remove(log_file)
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO, filename=log_file, filemode='w')
    print(
        "*" * 80, "",
        "All log statements are redirected to %s" % log_file, "",
        "*" * 80, "",
        "Displaying step summaries here:\n\n",
        "Simplifying FASTA sequences",
        sep="\n"
    )


# Gather all files to parse that match user-passed extensions
def _files_iter(ap: ArgParse, storage_dir: str) -> Generator[str, ArgParse, None]:
    for file in os.listdir(ap.args.fasta_directory):
        for ext in ap.args.extensions:
            if file.endswith(ext):
                yield _simplify_fasta(ap, file, storage_dir)
    return None


def _simplify_fasta(ap: ArgParse, file, storage_dir: str) -> str:
    # Simplify FASTA of complex-named sequences
    fasta_file = str(Path(os.path.join(ap.args.fasta_directory, file)).resolve())
    out_file = os.path.join(storage_dir, os.path.basename(os.path.splitext(fasta_file)[0]) + ".fna")
    record_p = SeqIO.parse(fasta_file, "fasta")
    i: int = 0
    records = []
    sys.stderr.write(file + "\n")
    for record in record_p:
        _i = str(i)
        _record_id = record.id
        for val in punctuation:
            _record_id = _record_id.replace(val, "")
        # Shorten id to 16 characters
        _len = len(_record_id)
        if _len > 16:
            _len = 16
        # Truncate and add unique number
        _record_id = str(_record_id[:_len - len(_i)]) + _i
        # Write ids to stderr for user
        sys.stderr.write(_record_id + "\t" + str(record.id) + "\n")
        # Store id
        record.id = _record_id
        records.append(record)
        i += 1
    SeqIO.write(records, out_file, "fasta")
    return out_file


# Create softlinks of final files to output directory
# Write summary file of reults
def _link_final_output(_output_files_list: List[List[object]], files_prefixes: List[str], _final_output_dir: str):
    # Store output paths as file for easy loading
    _paths_output_file = open(os.path.join(_final_output_dir, "paths_summary.tsv"), "w")
    for _files, _file_prefix in zip(_output_files_list, files_prefixes):
        _sub_out = os.path.join(_final_output_dir, _file_prefix)
        if not os.path.exists(_sub_out):
            os.makedirs(_sub_out)
        # Copy results to results dir for easier access
        for _file in _files:
            # Write info to file
            if isinstance(_file, dict):
                sorted_keys = sorted(list(_file.keys()))
                _paths_output_file.write(
                    "".join((
                        "\t".join((
                            _file_prefix,  # Name of record
                            *(os.path.join(_sub_out, str(_file[_f])) for _f in sorted_keys)  # Files for record
                        )), "\n"
                    ))
                )
            # Generate link of path
            elif isinstance(_file, str):
                local["ln", "-srf"][_file, _sub_out]()
        # Write actual paths for use in API, and as a means for users to copy files
        _paths_output_file.write("".join(("\t".join((os.path.join(_sub_out, str(_file)) for _file in _files)), "\n")))
    _paths_output_file.close()


# Parse user arguments
def _parse_args(ap: ArgParse, tm: TaskManager) -> Tuple[ConfigManager, bool]:
    # Confirm path existence
    assert os.path.exists(ap.args.config_file)
    is_continued = False
    # FASTA directory provided by user
    if ap.args.fasta_directory is not None:
        # Confirm that path exists
        assert os.path.exists(ap.args.fasta_directory)
    # Directory not provided, based on output directory
    else:
        assert os.path.exists(ap.args.output)
        # Generate list of files to parse
        is_continued = True
    # Ensure command is valid
    assert ap.args.command in tm.programs
    if ap.args.debug is True:
        ap.args.debug = 0
    else:
        ap.args.debug = 1
    # Determine file extensions to keep
    ap.args.extensions = ap.args.extensions.split("/")
    return ConfigManager(ap.args.config_file), is_continued


# # Driver logic
def _main(ap: ArgParse, cfg: ConfigManager, is_continued: bool, tm: TaskManager):
    # Generate primary path manager
    pm = PathManager(ap.args.output)
    # Begin logging
    _initialize_logging(ap)
    logging.info("Creating working directory")
    logging.info("Simplifying FASTA sequences")
    # Gather list of files to analyze
    if is_continued:
        # Gather from existing data
        input_files = _get_list_of_files(
            os.path.join(ap.args.output, "results/paths_summary.tsv"),
            tm.input_type[ap.args.command],
        )
        input_prefixes = [_prefix(_file) for _file in input_files]
    else:
        # Simplify FASTA files into working directory
        pm.add_dirs("MAGS")
        input_files = list(_file for _file in _files_iter(ap, pm.get_dir("MAGS")))
        # List of prefixes for tracking each file's progress
        input_prefixes = [_prefix(_file) for _file in input_files]
        # Create base dir for each file to analyze
        all([pm.add_dirs(_file) for _file in input_prefixes])

    # # Begin task list
    # Generate first task from list
    task_list = tm.programs[ap.args.command]
    task = task_list[0](cfg, input_files, pm, input_prefixes, ap.args.debug)
    for i in range(1, len(task_list)):
        # Run task
        task.run()
        task = task_list[i](*task.output())
    # Must call output on last task to generate final summary statistics
    task.run()
    output = task.output()
    # Create summary softlinks
    _link_final_output(output[1], output[3], os.path.join(ap.args.output, "results"))


# Get program-needed list of files for this step in pipeline
def _get_list_of_files(summary_file: str, file_type: str) -> List[str]:
    file_fp = open(summary_file, "r")
    _col_idx = next(file_fp).rstrip("\r\n").split("\t").index(file_type)
    return [
        line.rstrip("\r\n").split()[_col_idx]
        for line in file_fp
    ]


if __name__ == "__main__":
    # Redirect dask nanny errors
    signal(SIGPIPE, SIG_DFL)
    DEFAULT_EXTS = ".fna/.fasta/.fa"
    _tm = TaskManager()
    _ap = ArgParse(
        (
            (("command",),
             {"help": "Select from %s" % "/".join(_tm.programs)}),
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
        func("\nEukMetaSanity pipeline complete!")
