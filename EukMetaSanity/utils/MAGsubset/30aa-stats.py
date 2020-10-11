#!/usr/bin/env python3
import os
from Bio import SeqIO
from pathlib import Path
from decimal import Decimal
from dataclasses import dataclass
from collections import defaultdict
from typing import Optional, Dict, List
from EukMetaSanity.utils.arg_parse import ArgParse


# # Record represents amino acid sequence and its respective metadata
@dataclass
class Record:
    size: int
    annotation: str = ""
    euk_ms_match: str = ""
    annotation_bitscore: float = -1.0
    annotated: bool = False
    in_eukms_repeat: bool = False
    annotation_evalue: Decimal = Decimal("-1.0")

    def __eq__(self, other):
        return self.annotation == other.annotation

    def set_annotation(self, annotation_str: str):
        self.annotation = annotation_str
        self.annotated = True


class RecordSet:
    def __init__(self, rbh_path: str):
        self._data = defaultdict(Record)
        self._rbh = RecordSet.load_rbh_dict(rbh_path)

    # Primary logic
    def generate(self, fasta_file: str):
        # Insert base record information
        for record in SeqIO.parse(fasta_file, "fasta"):
            self._insert(record.id, Record(len(record.seq)))
        # Check if RBH present to EukMS
        return -1

    @property
    def data(self):
        return self._data

    def _insert(self, key: str, val: Record):
        self._data[key] = val

    def _find(self, key: str) -> Optional[Record]:
        return self._data.get(key, None)

    @staticmethod
    def load_rbh_dict(rbh_path: str):
        out = defaultdict(dict)
        for line in open(rbh_path, "r"):
            line = line.rstrip("\r\n").split("\t")
            out[line[0]].update({line[1]: {
                "pident": line[2],
                "bitscore": line[-1],
                "evalue": line[-2],
            }})
        return out

    @staticmethod
    def load_annotation_summary_file(annot_path: str):
        out = defaultdict(dict)
        for line in open(annot_path, "r"):
            pass


# Validate provided paths
def validate(_ap: ArgParse):
    for _file in (
            _ap.args.rbh_file,
            _ap.args.fasta_file,
            _ap.args.annotation_summary_file,
            _ap.args.repeats_gff3_file,
            _ap.args.file_a,
            _ap.args.file_b,
    ):
        if not os.path.exists(_file):
            print(_file + " does not exist!")
            exit(1)


if __name__ == "__main__":
    ap = ArgParse(
        (
            (("rbh_file",), {"help": "Path to mmseqs rbh output file"}),
            (("fasta_file",), {"help": "Path to amino acid FASTA file"}),
            (("annotation_summary_file",), {"help": "Path to .summary output file"}),
            (("repeats_gff3_file",), {"help": "Path to repeats .gff3 file"}),
            (("-a", "--file_a"), {"help": "Path to file matching ids in first column of rbh_file"}),
            (("-b", "--file_b"), {"help": "Path to file matching ids in second column of rbh_file"}),
        ),
        description="Summarize amino acids between old and new comparisons"
    )
    validate(ap)
    output_file = os.path.join(os.path.splitext(Path(ap.args.rbh_file).resolve())[0], ".cmp.tsv")
    mag_data = RecordSet(
        ap.args.rbh_file
    )

