#!/usr/bin/env python3
import os
import glob
import numpy as np
from typing import List
import matplotlib.pyplot as plt
from EukMetaSanity.utils.arg_parse import ArgParse


def parse_bitscores(infile: str) -> List[float]:
    return [float(line.rstrip("\r\n").split()[-1]) for line in open(infile, "r")]


if __name__ == "__main__":
    ap = ArgParse(
        (
            (("-a", "--a_files"),
             {"help": "Comparison file list a", "nargs": "+", "required": True}),
            (("-b", "--b_files"),
             {"help": "Comparison file list b", "nargs": "+"}),
        ),
        description="Generate matplotlib boxplots comparing two lists of files"
    )

    a = []
    b = []
    for file in ap.args.a_files:
        a.extend(parse_bitscores(file))
    for file in ap.args.b_files:
        b.extend(parse_bitscores(file))

    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    b = ax.boxplot(np.array([a, b]))
    print([item.get_ydata()[1] for item in b["whiskers"]])
    plt.show()
