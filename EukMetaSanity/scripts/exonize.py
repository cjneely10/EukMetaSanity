#!/usr/bin/env python3
import os
from Bio import SeqIO
from typing import List
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from EukMetaSanity.utils.arg_parse import ArgParse

init_mapping = {
    **{key: 0 for key in ("A", "C", "G", "T")},
    "N": 1,
}


def _parse_args(ap: ArgParse):
    for _path in (ap.args.fasta_file, ap.args.gff3_file):
        assert os.path.exists(_path)


def exonize(fasta_file: str, gff3_file: str, output_file: str):
    w = open(output_file, "w")
    _gff3_fp = open(gff3_file, "r")
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    _data = defaultdict(list)
    parse_data(_data, _gff3_fp)
    _out_data = defaultdict(list)


def parse_data(data, gff3_fp):
    for _line in gff3_fp:
        if _line.startswith("#"):
            continue
        line = _line.rstrip("\r\n").split("\t")
        data[line[0]].append((line[2], int(line[3]), int(line[4]), line[6]))


def flatten(data, out_data):
    for _id, feature_list in data.items():
        for feature_tuple in feature_list:
            pass


def generate_initial_region(record: SeqRecord, _region: List[int]):
    for i, val in enumerate(record.seq):
        _region[i] = init_mapping[val]


if __name__ == "__main__":
    _ap = ArgParse(
        (
            (("-f", "--fasta_file"),
             {"help": "Input FASTA file", "required": True}),
            (("-g", "--gff3_file"),
             {"help": "Input GFF3 file", "required": True}),
            (("-o", "--output_file"),
             {"help": "Output path, default stdout", "default": "/dev/stdout"}),
        ),
        description="Convert EukMetaSanity .merged.gff3 into exonized .nr.gff3 file"
    )
    _parse_args(_ap)
    exonize(_ap.args.fasta_file, _ap.args.gff3_file, _ap.args.output_file)
