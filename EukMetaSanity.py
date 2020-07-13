#!/usr/bin/env python3
import os
import logging
from Bio import SeqIO
from pathlib import Path
from typing import Generator
from string import punctuation
from signal import signal, SIGPIPE, SIG_DFL
from EukMetaSanity.utils.arg_parse import ArgParse
from EukMetaSanity.utils.path_manager import PathManager
from EukMetaSanity.utils.helpers import prefix as _prefix
from EukMetaSanity.utils.config_manager import ConfigManager
from EukMetaSanity.tasks.manager.task_manager import TaskManager

"""
EukMetaSanity - Generate structural/functional annotations for simple Eukaryotes

"""


# # Available programs
# Return task-list for run command
def _run_iter(tm: TaskManager, program: str) -> Generator[type, TaskManager, None]:
    task_list = tm.programs[program]
    for task in task_list:
        yield task


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
    for record in record_p:
        _i = str(i)
        for val in punctuation:
            record.id = record.id.replace(val, "")
        # Shorten id to 16 characters
        if len(str(record.id)) > 16:
            # Store old id in description
            record.description = str(record.id) + " " + record.description
            # Update new id
            record.id = str(record.id[:16 - len(_i)]) + _i
        records.append(record)
        i += 1
    SeqIO.write(records, out_file, "fasta")
    return out_file


# Parse user arguments
def _parse_args(ap: ArgParse, tm: TaskManager) -> ConfigManager:
    # Confirm path existence
    assert os.path.exists(ap.args.config_file)
    assert os.path.exists(ap.args.fasta_directory)
    # Ensure command is valid
    assert ap.args.command in tm.programs or ap.args.command == "new_cfg"
    if ap.args.debug is True:
        ap.args.debug = 0
    else:
        ap.args.debug = 1
    # Determine file extensions to keep
    ap.args.extensions = ap.args.extensions.split("/")
    return ConfigManager(ap.args.config_file)


# # Driver logic
def _main(ap: ArgParse, cfg: ConfigManager, tm: TaskManager):
    # Generate primary path manager
    pm = PathManager(ap.args.output)
    # Begin logging
    _initialize_logging(ap)
    logging.info("Creating working directory")
    logging.info("Simplifying FASTA sequences")
    # Gather list of files to analyze and simplify FASTA files into working directory
    pm.add_dirs("MAGS")
    input_files = list(_file for _file in _files_iter(ap, pm.get_dir("MAGS")))
    # List of prefixes for tracking each file's progress
    input_prefixes = [_prefix(_file) for _file in input_files]
    # Create base dir for each file to analyze
    all([pm.add_dirs(_file) for _file in input_prefixes])

    # # Begin task list
    # Generate first task from list
    run_iter = _run_iter(tm, ap.args.command)
    task = next(run_iter)(cfg, input_files, pm, input_prefixes, ap.args.debug)
    # Run task
    task.run()
    # Primary program loop
    while True:
        try:
            # Run next task
            task = next(run_iter)(*task.output())
            task.run()
        except StopIteration:
            break
    # Must call output on last task to generate final summary statistics
    task.output()


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
             {"help": "Directory of FASTA files to annotate", "required": True}),
            (("-c", "--config_file"),
             {"help": "Config file", "required": True}),
            (("-x", "--extensions"),
             {"help": "Gather files matching '/'-separated list of extensions, default %s" % DEFAULT_EXTS,
              "default": DEFAULT_EXTS}),
            (("-o", "--output"),
             {"help": "Output directory, default out", "default": "out"}),
            (("-d", "--debug"),
             {"help": "Developer mode: display all commands on single thread, default False", "default": False,
              "action": "store_true"}),
        ),
        description="Run EukMetaSanity pipeline"
    )
    _main(_ap, _parse_args(_ap, _tm), _tm)
    for func in (logging.info, print):
        func("\nEukMetaSanity pipeline complete!")
