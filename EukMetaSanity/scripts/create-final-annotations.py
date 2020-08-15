#!/usr/bin/env python3
import itertools
import os
import re
from operator import itemgetter

from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Optional, Tuple, Generator, List, Callable, Dict
from EukMetaSanity.utils.arg_parse import ArgParse


class Gff3Parser:
    def __init__(self, gff3_file: str, priority: str = "metaeuk"):
        self.fp = open(gff3_file, "r")
        self.priority = getattr(Gff3Parser, priority, lambda _: _)

    def next_gene(self) -> Generator[Dict[str, object], None, None]:
        _line: str
        line: List[str]
        try:
            line = next(self.fp).rstrip("\r\n").split("\t")
        except StopIteration:
            return
        while True:
            if line[0][0] == "#":
                line = next(self.fp).rstrip("\r\n").split("\t")
                continue
            elif line[0][0] == ">":
                break
            if line[2] == "locus":
                # Putative gene
                gene_data = {
                    "fasta-id": line[0],
                    "start": line[3],
                    "end": line[4],
                    "strand": line[6],
                }
                transcripts = []
                line = next(self.fp).rstrip("\r\n").split("\t")
                while line[2] != "locus":
                    # Read in transcript info
                    transcripts.append(
                        [line[1], line[3], line[4], []]  # First line is a transcript: ID,source,start,end
                    )
                    line = next(self.fp).rstrip("\r\n").split("\t")
                    # Add exon to current info
                    while line[2] not in ("transcript", "locus"):
                        transcripts[-1][-1].append(
                            (line[3], line[4])
                        )
                        line = next(self.fp).rstrip("\r\n").split("\t")
                # Filter for specific transcripts
                gene_data["transcripts"] = self.priority(transcripts)
                yield gene_data

    @staticmethod
    def metaeuk(line: List[List]) -> List[List]:
        return Gff3Parser.filter_specific(line, "metaeuk")

    @staticmethod
    def ab_initio(line: List[List]) -> List[List]:
        return Gff3Parser.filter_specific(line, "ab-initio")

    @staticmethod
    def filter_specific(line: List[List], name: str) -> List[List]:
        for _l in line:
            if _l[0] == name:
                return _l
        return []

    @staticmethod
    def merge(line: List[List]) -> List[List]:
        # Sort coordinates by start value
        ranges_in_coords = sorted(itertools.chain(*[l[-1] for l in line]), key=itemgetter(0))
        # Will group together matching sections into spans
        # Return list of these spans at end
        # Initialize current span and list to return
        spans_in_coords = [list(ranges_in_coords[0]), ]
        for coords in ranges_in_coords[1:]:
            # The start value is within the farthest range of current span
            # and the end value extends past the current span
            if coords[0] <= spans_in_coords[-1][1] < coords[1]:
                spans_in_coords[-1][1] = coords[1]
            # The start value is past the range of the current span
            # Append old span to list to return
            # Reset current span to this new range
            elif coords[0] > spans_in_coords[-1][1]:
                spans_in_coords.append(list(coords))
        return spans_in_coords

    def _gene_to_string(self) -> str:
        pass

    def __str__(self) -> str:
        pass

    def __repr__(self) -> str:
        return self.__str__()

    def __iter__(self):
        return self.next_gene()


def convert_final_gff3(gff3_file: str, fasta_file: str, filter_function: str):
    gff3 = Gff3Parser(gff3_file, filter_function)
    for gene in gff3:
        print(gene)


def find_orf(record: SeqRecord) -> Tuple[Optional[SeqRecord], int]:
    longest = (0,)
    if record.description[-1] == "+":
        nuc = str(record.seq)
    else:
        nuc = str(record.reverse_complement().seq)
    for m in re.finditer("ATG", nuc):
        pro = Seq(nuc)[m.start():].translate(to_stop=True)
        if len(pro) > longest[0]:
            longest = (len(pro), m.start(), str(pro))
    if longest[0] > 0:
        return SeqRecord(
            seq=Seq(str(longest[2]) + "*"),
            id=record.id,
            description=record.description
        ), longest[1]
    return None, 0


if __name__ == "__main__":
    ap = ArgParse(
        (
            (("command", ),
             {"help": "Filter command, select from metaeuk/abinitio/merge"}),
            (("-g", "--gff3_file"),
             {"help": ".tmp.nr.gff3 file", "required": True}),
            (("-f", "--fasta_file"),
             {"help": "FASTA file", "required": True}),
        ),
        description="Convert .tmp.nr.gff3 to final .gff3 output"
    )
    for _file in (ap.args.gff3_file, ap.args.fasta_file):
        assert os.path.exists(_file)
    convert_final_gff3(ap.args.gff3_file, ap.args.fasta_file, ap.args.command)

if __name__ == "__main__":
    pass
