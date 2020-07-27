#!/usr/bin/env python3
import os
from Bio import SeqIO
from typing import List, Dict
from collections import namedtuple
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from EukMetaSanity.utils.arg_parse import ArgParse

GFFCoord = namedtuple("GFFCoord", ("evidence", "loc_type", "start", "end", "strand", "parent_id"))


class Coordinate:
    __slots__ = (
        "is_repeat_region",
        "is_exon",
        "is_intron",
        "parent_id",
        "evidence",
        "strand",
        "nucleotide",
    )

    def __init__(self, nucleotide: str):
        if nucleotide == "N":
            self.is_repeat_region = True
        else:
            self.is_repeat_region = False
        self.is_exon = False
        self.is_intron = False
        self.parent_id = None
        self.evidence = []
        self.strand = None
        self.nucleotide = nucleotide


# Convert gff3 file to dictionary
def gff3_to_dict(gff3_file: str) -> Dict[str, List[GFFCoord]]:
    gff3_fp = open(gff3_file, "r")
    out_data = defaultdict(list)
    for _line in gff3_fp:
        if _line.startswith("#"):
            continue
        line = _line.rstrip("\r\n").split()
        out_data[line[0]].append(
            GFFCoord(
                evidence=line[1],
                loc_type=line[2],
                start=int(line[3]),
                end=int(line[4]),
                strand=line[6],
                parent_id=line[8].split(";")[0].split("=")[1]
            )
        )
    return dict(out_data)


# Create list of Coordinate objects that describe each location on a chromosome
def generate_initial_region(record: SeqRecord) -> List[Coordinate]:
    return [
        Coordinate(val)
        for val in record.seq
    ]


# # Program driver logic

def _parse_args(ap: ArgParse):
    for _path in (ap.args.fasta_file, *ap.args.gff3_files):
        assert os.path.exists(_path)


def exonize(fasta_file: str, gff3_files: List[str], output_file: str):
    w = open(output_file, "w")
    # Get FASTA file as dict
    record_p = SeqIO.parse(fasta_file, "fasta")
    # Generate list of GFF data dictionaries
    gff_dict_list = [
        gff3_to_dict(_file) for _file in gff3_files
    ]
    # Iterate over each record
    for record in record_p:
        # Create bare region
        region = generate_initial_region(record)
        for gff_dict in gff_dict_list:
            pass

    w.close()


if __name__ == "__main__":
    _ap = ArgParse(
        (
            (("-f", "--fasta_file"),
             {"help": "Input FASTA file", "required": True}),
            (("-g", "--gff3_files"),
             {"help": "Input GFF3 files", "required": True, "nargs": "+"}),
            (("-o", "--output_file"),
             {"help": "Output path, default stdout", "default": "/dev/stdout"}),
        ),
        description="Convert EukMetaSanity .merged.gff3 into exonized .nr.gff3 file"
    )
    _parse_args(_ap)
    exonize(_ap.args.fasta_file, _ap.args.gff3_files, _ap.args.output_file)
