#!/usr/bin/env python3
import os
import logging
from Bio import SeqIO
from pathlib import Path
from typing import Generator
from string import punctuation
from signal import signal, SIGPIPE, SIG_DFL
from EukMetaSanity.src.utils.arg_parse import ArgParse
from EukMetaSanity.src.utils.path_manager import PathManager
from EukMetaSanity.src.tasks.task_manager import TaskManager
from EukMetaSanity.src.utils.config_manager import ConfigManager

"""
EukMetaSanity - Generate structural/functional annotations for simple Eukaryotes

"""


# # Available programs
# Return task-list for run command
def _run_iter(tm: TaskManager, program: str) -> Generator[type, TaskManager, None]:
    task_list = tm.tasks[program]
    for task in task_list:
        yield task


# # Helper functions
# Get prefix of path - e.g. for /path/to/file_1.ext, return file_1
def _prefix(_path: str) -> str:
    return ".".join(_path.split("/")[-1].split(".")[:-1])


# Logging initialize
def _initialize_logging(ap: ArgParse) -> str:
    # Initialize logging
    log_file = os.path.join(ap.args.output, "eukmetasanity.log")
    if os.path.exists(log_file):
        os.remove(log_file)
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO, filename=log_file, filemode='w')
    return log_file


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
    for record in record_p:
        SeqIO.write(record_generator(record, str(i)), out_file, "fasta")
        i += 1
    return out_file


def record_generator(record: type, _i: str):
    # Remove problem characters
    for val in punctuation:
        record.id = record.id.replace(val, "")
    # Shorten id to 16 characters
    if len(str(record.id)) > 16:
        record.id = record.id[:16 - len(_i)] + _i
    yield record


# Parse user arguments
def _parse_args(ap: ArgParse, tm: TaskManager) -> ConfigManager:
    # Confirm path existence
    assert os.path.exists(ap.args.config_file)
    assert os.path.exists(ap.args.fasta_directory)
    # Ensure command is valid
    assert ap.args.command in tm.tasks
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
    logging_path = _initialize_logging(ap)
    print(
        "*" * 80, "",
        "All log statements are redirected to %s" % logging_path, "",
        "*" * 80, "",
        "Displaying step summaries here:\n\n",
        "Simplifying FASTA sequences",
        sep="\n"
    )
    logging.info("Creating working directory")
    logging.info("Simplifying FASTA sequences")
    # Gather list of files to analyze and simplify FASTA files
    # Move all input files to output directory
    pm.add_dirs("MAGS")
    input_files = list(_file for _file in _files_iter(ap, pm.get_dir("MAGS")))
    input_prefixes = [_prefix(_file) for _file in input_files]
    # Create base dir for each file to analyze
    all([pm.add_dirs(_file) for _file in input_prefixes])

    # # Begin task list
    # Generate first task from list
    run_iter = _run_iter(tm, ap.args.command)
    task = next(run_iter)(input_files, cfg, pm, input_prefixes, ap.args.debug)
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
    _ap = ArgParse(
        (
            (("command",),
             {"help": "Select from run/refine"}),
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

    _tm = TaskManager()
    _main(_ap, _parse_args(_ap, _tm), _tm)
    logging.info("EukMetaSanity pipeline complete!")
