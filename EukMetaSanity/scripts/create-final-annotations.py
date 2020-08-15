#!/usr/bin/env python3
import os
import re
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Optional, Tuple, Generator, List, Callable, Dict
from EukMetaSanity.utils.arg_parse import ArgParse


class Gff3Parser:
    def __init__(self, gff3_file: str, priority: str = "metaeuk"):
        self.fp = open(gff3_file, "r")
        self.priority = priority

    def next_gene(self) -> Generator[Dict[str, object], None, None]:
        _line: str
        line: List[str]
        try:
            line = next(self.fp).rstrip("\r\n").split("\t")
        except StopIteration:
            return
        while True:
            if line[2] == "locus":
                # Putative gene
                gene_data = {
                    "id": line[0],
                    "start": line[3],
                    "end": line[4],
                    "strand": line[6],
                    "overlapped": int(line[1].split("num_samples=")[1].split(";")[0]) > 1,
                }
                # Read in transcript info
                transcripts = []
                line = next(self.fp).rstrip("\r\n").split("\t")
                while line[2] != "locus":
                    line = next(self.fp).rstrip("\r\n").split("\t")
                
    def _gene_to_string(self) -> str:
        pass

    def __str__(self) -> str:
        pass

    def __repr__(self) -> str:
        return self.__str__()


def _parse_args(ap: ArgParse):
    for _file in (ap.args.gff3_file, ap.args.fasta_file):
        assert os.path.exists(_file)


def convert_final_gff3(gff3_file: str, fasta_file: str):
    pass


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
    _ap = ArgParse(
        (
            (("-g", "--gff3_file"),
             {"help": ".tmp.nr.gff3 file", "required": True}),
            (("-f", "--fasta_file"),
             {"help": "FASTA file", "required": True}),
        ),
        description="Convert .tmp.nr.gff3 to final .gff3 output"
    )
    _parse_args(_ap)
    convert_final_gff3(_ap.args.gff3_file, _ap.args.fasta_file)


if __name__ == "__main__":
    pass
