#!/usr/bin/env python3
from EukMetaSanity.src.utils.arg_parse import ArgParse


def _parse_args(ap):
    pass


def _main(ap):
    pass


if __name__ == "__main__":
    _ap = ArgParse(
        (
            (("command",),
             {"help": "Select from run/refine"}),
            (("fasta_directory",),
             {"help": "Directory of FASTA files to annotate"}),
            (("-c", "--config_file"),
             {"help": "Config file"}),
            ((), {}),
            ((), {}),
            ((), {}),
        ),
        description="Run EukMetaSanity pipeline"
    )
    _parse_args(_ap)
