#!/usr/bin/env python3
import os
from Bio import SeqIO
from pathlib import Path
from typing import Optional, Tuple
from collections import defaultdict
from dataclasses import dataclass, field
from EukMetaSanity.utils.arg_parse import ArgParse


# # Record represents amino acid sequence and its respective metadata
@dataclass
class Record:
    size: int  #
    annotated: bool = False  #
    euk_ms_match: str = ""  #
    in_eukms_repeat: bool = False  #
    euk_ms_match_score: float = 0.0  #
    exons: list = field(default_factory=list)  #
    annotation: dict = field(default_factory=dict)  #
    euk_ms_annotation: dict = field(default_factory=dict)  #

    def __eq__(self, other):
        return self.annotation == other.annotation

    def set_annotation(self, annotation_str: str):
        self.annotation = annotation_str
        self.annotated = True


class RecordSet:
    def __init__(self, rbh_path: str, annot_path: str, repeats_gff3: str, annotation_gff3: str):
        self._data = defaultdict(Record)
        self._rbh = RecordSet.load_rbh_dict(rbh_path)
        self._annotations = RecordSet.load_annotation_summary_file(annot_path)
        self._repeats = RecordSet.load_regions(repeats_gff3)
        self._gene_regions = RecordSet.load_regions(annotation_gff3, "gene")

    # Primary logic
    def generate(self, fasta_file: str, eukms_annotation_file: str):
        # # Insert base record information
        # id and size
        for record in SeqIO.parse(fasta_file, "fasta"):
            self._insert(record.id, Record(len(record.seq)))
            # Location on contig
            for location_tuple in self._gene_regions.get(record.id, []):
                self._data[record.id].exons.append(location_tuple)
            # print(self._data[record.id].exons)
        # # Update matches to EukMS calls
        # EukMS max RBH match
        for _id, annotation_dict in self._rbh.items():
            for match, match_data in annotation_dict.items():
                match_val = float(match_data["pident"])
                if self._data[_id].euk_ms_match_score < match_val:
                    self._data[_id].euk_ms_match_score = match_val
                    self._data[_id].euk_ms_match = match
        # In a EukMS repeat
        for contig_id, repeat_regions in self._repeats.items():
            for _id in self._data.keys():
                if contig_id in _id:
                    for exon in self._data[_id].exons:
                        for repeat_region in repeat_regions:
                            if RecordSet.overlap(exon, repeat_region):
                                self._data[_id].in_eukms_repeat = True
                                break
        # Store annotation
        for _id in self._data.keys():
            annotation = self._annotations[_id]
            if annotation != dict({}):
                self._data[_id].annotation = annotation
                self._data[_id].annotated = True
        # Compare to EukMS annotation
        eukms_annotations = RecordSet.load_annotation_summary_file(eukms_annotation_file)
        for _id in self._data.keys():
            annotation = eukms_annotations[self._data[_id].euk_ms_match]
            if annotation != dict({}):
                self._data[_id].euk_ms_annotation = annotation

    # def write(self, output_path: str):
    #     W = open(output_file, "w")
    #     for _id
    #
    #     W.close()

    @property
    def data(self):
        return self._data.values()

    def _insert(self, key: str, val: Record):
        self._data[key] = val

    def _find(self, key: str) -> Optional[Record]:
        return self._data.get(key, None)

    @staticmethod
    def overlap(a: Tuple[int, int], b: Tuple[int, int]) -> bool:
        return max(a[0], b[0]) <= min(a[1], b[1])

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
        R = open(annot_path, "r")
        # Skip top counts
        next(R)
        # Parse first line for identifiers
        idx = next(R).rstrip("\r\n").split("\t")[1:]
        for line in R:
            line = line.rstrip("\r\n").split("\t")
            # Store as dict
            to_add = {idx[i]: line[i + 1] for i in range(len(idx))}
            if to_add != {idx[i]: '0' for i in range(len(idx))}:
                out[line[0]] = to_add
        return out

    @staticmethod
    def load_regions(repeats_gff3: str, filter_id: str = ""):
        out = defaultdict(set)
        # Store as dict of contig id to set of tuple ranges of repeated regions
        for line in open(repeats_gff3, "r"):
            if len(line) > 0 and line[0] != "#":
                if line[0] == '>':
                    break
                line = line.rstrip("\r\n").split("\t")
                if filter_id == "":
                    out[line[0]].add((int(line[3]), int(line[4])))
                elif line[2] == filter_id:
                    out[line[-1].replace("ID=", "").split(";")[0]].add((int(line[3]), int(line[4])))
        return out


# Validate provided paths
def validate(_ap: ArgParse):
    for _file in (
            _ap.args.rbh_file,
            _ap.args.fasta_file,
            _ap.args.repeats_gff3_file,
            _ap.args.file_a,
            _ap.args.file_b,
            _ap.args.gff3_file,
    ):
        if not os.path.exists(_file):
            print(_file + " does not exist!")
            exit(1)


if __name__ == "__main__":
    ap = ArgParse(
        (
            (("rbh_file",), {"help": "Path to mmseqs rbh output file"}),
            (("fasta_file",), {"help": "Path to amino acid FASTA file"}),
            (("repeats_gff3_file",), {"help": "Path to repeats .gff3 file"}),
            (("-a", "--file_a"), {"help": "Path to .summary file matching ids in first column of rbh_file (GM file)"}),
            (("-b", "--file_b"), {"help": "Path to .summary file matching ids in second column of rbh_file (EukMS file)"}),
            (("-g", "--gff3_file"), {"help": "Original (old) gff3 annotation file"}),
        ),
        description="Summarize amino acids between old and new comparisons"
    )
    validate(ap)
    output_file = os.path.join(os.path.splitext(Path(ap.args.rbh_file).resolve())[0], ".cmp.tsv")
    mag_data = RecordSet(ap.args.rbh_file, ap.args.file_a, ap.args.repeats_gff3_file, ap.args.gff3_file)
    mag_data.generate(ap.args.fasta_file, ap.args.file_b)
    for record in mag_data.data:
        if record.in_eukms_repeat:
            print(record)
    # mag_data.write(output_file)

