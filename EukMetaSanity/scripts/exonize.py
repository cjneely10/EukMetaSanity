#!/usr/bin/env python3
import os
from Bio import SeqIO
from EukMetaSanity.utils.arg_parse import ArgParse


def _parse_args(ap: ArgParse):
    for _path in (ap.args.fasta_file, ap.args.gff3_file):
        assert os.path.exists(_path)


def exonize(fasta_file: str, gff3_file: str, output_file: str):
    w = open(output_file, "w")
    gff3_fp = open(gff3_file, "r")
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    w.write("##gff-version 3\n")
    for line in gff3_fp:
        if line.startswith("#"):
            continue
        line = line.rstrip("\r\n").split("\t")



if __name__ == "__main__":
    _ap = ArgParse(
        (
            (("-f", "--fasta_file"),
             {"help": "Input FASTA file", "required": True}),
            (("-g", "--gff3_file"),
             {"help": "Input GFF3 file", "required": True}),
            (("-o", "--output"),
             {"help": "Output path, default stdout", "default": "/dev/stdout"}),
        ),
        description="Convert EukMetaSanity .merged.gff3 into exonized .nr.gff3 file"
    )
    _parse_args(_ap)
    exonize(_ap.args.fasta_file, _ap.args.gff3_file, _ap.args.output_file)
