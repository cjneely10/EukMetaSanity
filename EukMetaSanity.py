#!/usr/bin/env python3
import os
from dask.distributed import Client, wait
from signal import signal, SIGPIPE, SIG_DFL
from EukMetaSanity.src.tasks.taxonomy import Taxonomy
from EukMetaSanity.src.utils.arg_parse import ArgParse
from EukMetaSanity.src.utils.path_manager import PathManager
from EukMetaSanity.src.utils.config_manager import ConfigManager


# # Available programs

# Return task-list for run command
def run():
    task_list = (Taxonomy,)
    for task in task_list:
        yield task


# # Helper functions

# Gather all files to parse that match user-passed extensions
def files_iter(ap):
    for file in os.listdir(ap.args.fasta_directory):
        for ext in ap.args.extensions:
            if file.endswith(ext):
                yield os.path.join(ap.args.fasta_directory, file)


# Parse user arguments
def _parse_args(ap):
    # Confirm path existence
    assert os.path.exists(ap.args.config_file)
    assert os.path.exists(ap.args.fasta_directory)
    # Ensure command is valid
    assert ap.args.command in ("run", "refine")
    # Determine file extensions to keep
    ap.args.extensions = ap.args.extensions.split("/")
    return ConfigManager(ap.args.config_file)


# # Driver logic
def _main(ap, cfg):
    # Generate primary path manager
    pm = PathManager(ap.args.output)
    for _file in files_iter(ap):
        prefix = os.path.basename(os.path.splitext(_file)[0])
        pm.add_dir(prefix)
        task_list = []
        for task in run():
            task_list.append(
                task(
                    {"in": _file},
                    cfg,
                    pm,
                    prefix,
                )
            )


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
