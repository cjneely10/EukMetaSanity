#!/usr/bin/env python

"""
EukMetaSanity.py
===============================================
EukMetaSanity - Generate structural/functional annotations for Eukaryotes

"""

import os
import logging
from signal import signal, SIGPIPE, SIG_DFL
from EukMetaSanity.utils.arg_parse import ArgParse
from EukMetaSanity.tasks.base.path_manager import PathManager
from EukMetaSanity.tasks.base.task_manager import TaskManager
from EukMetaSanity.tasks.base.input_manager import InputManager
from EukMetaSanity.tasks.base.config_manager import ConfigManager
from EukMetaSanity.tasks.manager.pipeline_manager import PipelineManager


def _initialize_logging(ap: ArgParse):
    """ Logging initialize

    :param ap: ArgParse object
    """
    # Initialize logging
    log_file = os.path.join(ap.args.output, "%s-eukmetasanity.log" % ap.args.command)
    if os.path.exists(log_file):
        os.remove(log_file)
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO, filename=log_file, filemode='w')
    print(
        "*" * 80, "",
        "Primary log statements are redirected to %s" % log_file,
        "Task-level log statements are redirected to subdirectory log files", "",
        "*" * 80, "",
        "Displaying step summaries here:\n",
        sep="\n"
    )


def _main(ap: ArgParse, cfg: ConfigManager, tpm: PipelineManager, pm: PathManager, im: InputManager):
    """ Driver logic

    :param ap: ArgParse object
    :param cfg: Loaded config file manager
    :param tpm: Loaded pipeline manager
    """
    # Begin logging
    _initialize_logging(ap)
    # # Begin task list
    tm = TaskManager(tpm, cfg, pm, im.input_files, im.input_prefixes, ap.args.debug, ap.args.command)
    tm.run(ap.args.output)


# Parse user arguments
def _parse_args(ap: ArgParse, tm: PipelineManager) -> ConfigManager:
    """ Parse file paths. Ensure pipeline requested by user is a valid pipeline

    :param ap: ArgParse object
    :param tm: PipelineManager object
    :return: Successfully loaded ConfigManager object
    """
    # Confirm path existence
    assert os.path.exists(ap.args.config_file)
    # Ensure command is valid
    assert ap.args.command in tm.programs.keys()
    if ap.args.debug is True:
        ap.args.debug = 0
    else:
        ap.args.debug = 1
    # Determine file extensions to keep
    ap.args.extensions = ap.args.extensions.split("/")
    return ConfigManager(ap.args.config_file)


if __name__ == "__main__":
    # Redirect dask nanny errors
    signal(SIGPIPE, SIG_DFL)
    DEFAULT_EXTS = ".fna/.fasta/.fa"
    _tm = PipelineManager()
    # TODO: UI methods - clean step, soft-run to show what will launch
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
    _cfg: ConfigManager = _parse_args(_ap, _tm)
    # Generate primary path manager
    _pm = PathManager(_ap.args.output)
    # Gather list of files to analyze
    _im = InputManager(_ap.args.output, _ap.args.fasta_directory, _pm, _cfg, _ap.args.extensions)
    # Run main program logic
    _main(_ap, _cfg, _tm, _pm, _im)
    # Display final output line
    for func in (logging.info, print):
        func("\nEukMetaSanity %s pipeline complete!" % _ap.args.command)
