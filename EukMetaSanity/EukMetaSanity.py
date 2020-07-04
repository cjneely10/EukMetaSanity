#!/usr/bin/env python3
import os
from signal import signal, SIGPIPE, SIG_DFL
from EukMetaSanity.src.tasks.taxonomy import Taxonomy
from EukMetaSanity.src.utils.arg_parse import ArgParse
from EukMetaSanity.src.utils.config_manager import ConfigManager


def build_task_list(_type, cfg, ap):
    globals()["confirm_" + _type](cfg)


def confirm_run(cfg):
    pass


def _parse_args(ap):
    # Confirm path existence
    assert os.path.exists(ap.args.config_file)
    assert os.path.exists(ap.args.fasta_directory)
    # Ensure command is valid
    assert ap.args.command in ("run", "refine")
    # Determine file extensions to keep
    ap.args.extensions = ap.args.extensions.split("/")
    return ConfigManager(ap.args.config_file)


def _main(ap, cfg):
    # Generate primary fasta manager

    pass


if __name__ == "__main__":
    # Redirect dask nanny errors
    signal(SIGPIPE, SIG_DFL)
    DEFAULT_EXTS = ".fna/.fasta/.fa"
    _ap = ArgParse(
        (
            (("command",),
             {"help": "Select from run/refine"}),
            (("fasta_directory",),
             {"help": "Directory of FASTA files to annotate"}),
            (("-c", "--config_file"),
             {"help": "Config file"}),
            (("-x", "--extensions"),
             {"help": "Gather files matching '/'-separated list of extensions, default %s" % DEFAULT_EXTS,
              "default": DEFAULT_EXTS}),
            (("-o", "--output"),
             {"help": "Output directory, default out", "default": "out"}),
        ),
        description="Run EukMetaSanity pipeline"
    )
    _main(_ap, _parse_args(_ap))
