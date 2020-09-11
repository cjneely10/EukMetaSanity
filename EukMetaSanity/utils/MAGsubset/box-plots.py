#!/usr/bin/env python3
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
    _b = []
    for file in ap.args.a_files:
        a.extend(parse_bitscores(file))
    for file in ap.args.b_files:
        _b.extend(parse_bitscores(file))
    out = [a]
    if len(_b) > 0:
        out.append(_b)

    fig = plt.figure()
    b = plt.boxplot(np.array([a, _b]))
    print([item.get_ydata()[1] for item in b["whiskers"]])
    print([item.get_ydata() for item in b["boxes"]])
    print([item.get_ydata()[1] for item in b["medians"]])
    plt.xticks([0, 1, 2], ["", "GeneMark", "+MetaEuk"])
    plt.ylabel("Bitscore")
    plt.show()
