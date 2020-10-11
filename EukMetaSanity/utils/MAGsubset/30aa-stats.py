#!/usr/bin/env python3
import os
import squarify
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from EukMetaSanity.utils.arg_parse import ArgParse


def validate(_ap: ArgParse):
    for _file in (_ap.args.file_a, _ap.args.file_b, _ap.args.rbh_file):
        if not os.path.exists(_file):
            print(_file + " does not exist!")
            exit(1)


if __name__ == "__main__":
    ap = ArgParse(
        (
            (("rbh_file",),
             {"help": "Path to mmseqs rbh output file"}),
            (("-a", "--file_a"),
             {"help": "Path to file matching ids in first column of rbh_file"}),
            (("-b", "--file_b"),
             {"help": "Path to file matching ids in second column of rbh_file"}),
        ),
        description="Summarize amino acids between old and new comparisons"
    )
    validate(ap)
    output_file = os.path.join(os.path.splitext(Path(ap.args.rbh_file).resolve())[0], ".cmp.tsv")
