#!/usr/bin/env python3
import os
import logging
from numba import jit, types
from dask.distributed import Client, wait
from signal import signal, SIGPIPE, SIG_DFL
from EukMetaSanity.src.utils.arg_parse import ArgParse
from EukMetaSanity.src.tasks.taxonomy import TaxonomyIter
from EukMetaSanity.src.utils.path_manager import PathManager
from EukMetaSanity.src.utils.config_manager import ConfigManager


# # Available programs
# Return task-list for run command
def _run_iter():
    task_list = (TaxonomyIter,)
    for task in task_list:
        yield task


# # Helper functions
# Get prefix of path - e.g. for /path/to/file_1.ext, return file_1
@jit(types.unicode_type(types.unicode_type), nopython=True, cache=True)
def _prefix(_path: str):
    return ".".join(_path.split("/")[-1].split(".")[:-1])


# Gather all files to parse that match user-passed extensions
def _files_iter(ap: ArgParse):
    for file in os.listdir(ap.args.fasta_directory):
        for ext in ap.args.extensions:
            if file.endswith(ext):
                yield os.path.join(ap.args.fasta_directory, file)


# Parse user arguments
def _parse_args(ap: ArgParse):
    # Confirm path existence
    assert os.path.exists(ap.args.config_file)
    assert os.path.exists(ap.args.fasta_directory)
    # Ensure command is valid
    assert ap.args.command in ("run", "refine")
    # Determine file extensions to keep
    ap.args.extensions = ap.args.extensions.split("/")
    return ConfigManager(ap.args.config_file)


# # Driver logic
def _main(ap: ArgParse, cfg: ConfigManager):
    # Generate primary path manager
    pm = PathManager(ap.args.output)
    # Simplify FASTA files
    # Create base dir for each
    input_files = list(_files_iter(ap))
    input_prefixes = [_prefix(_file) for _file in input_files]
    all([pm.add_dirs(_file) for _file in input_prefixes])
    # Populate and call each task sublist
    for T_Task in _run_iter():
        task = T_Task(
            [{"-in": _file} for _file in input_files],
            cfg,
            pm,
            input_prefixes
        )
        # Call task list tasks
        task.run()


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
        ),
        description="Run EukMetaSanity pipeline"
    )
    _main(_ap, _parse_args(_ap))
