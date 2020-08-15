#!/usr/bin/env python3
import os
import re
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
from operator import itemgetter
from Bio.SeqRecord import SeqRecord
from typing import Optional, Tuple, Generator, List, Dict
from EukMetaSanity.utils.arg_parse import ArgParse


class Gff3Parser:
    def __init__(self, gff3_file: str, fasta_file: str, priority: str = "metaeuk"):
        assert priority in dir(Gff3Parser)
        self.fp = open(gff3_file, "r")
        self.priority = getattr(Gff3Parser, priority, lambda _: _)
        self.version = priority
        self.fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        self.count = 0

    def next_gene(self) -> Generator[Tuple[Dict[str, object], Optional[SeqRecord], Optional[SeqRecord]], None, None]:
        _line: str
        line: List[str]
        line = next(self.fp).rstrip("\r\n").split("\t")
        while True:
            if line[0][0] == "#":
                line = next(self.fp).rstrip("\r\n").split("\t")
                continue
            elif line[0][0] == ">":
                break
            if line[2] == "locus":
                self.count += 1
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
                if len(gene_data["transcripts"]) == 0:
                    continue
                seq = StringIO()
                orig_seq = self.fasta_dict[gene_data["fasta-id"]]
                if gene_data["strand"] == "-":
                    orig_seq = orig_seq.reverse_complement()
                orig_seq = str(orig_seq.seq)
                for transcript in gene_data["transcripts"]:
                    seq.write(orig_seq[int(transcript[0]) - 1: int(transcript[1])])
                record = SeqRecord(
                    seq=Seq(seq.getvalue()),
                    id="gene%s" % str(self.count),
                    description="strand=%s" % gene_data["strand"],
                    name="",
                )
                orf, offset = find_orf(record)
                yield (
                    self._gene_to_string(gene_data, str(offset)),
                    record,
                    orf,
                )

    @staticmethod
    def metaeuk(line: List[List]) -> List[List]:
        return Gff3Parser.filter_specific(line, "metaeuk")

    @staticmethod
    def abinitio(line: List[List]) -> List[List]:
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
        spans_in_coords = [list(ranges_in_coords[0]), ]
        for coords in ranges_in_coords[1:]:
            # The start value is within the farthest range of current span and the end value extends past current
            if coords[0] <= spans_in_coords[-1][1] < coords[1]:
                spans_in_coords[-1][1] = coords[1]
            # The start value is past the range of the current span, append old span to list to return
            elif coords[0] > spans_in_coords[-1][1]:
                spans_in_coords.append(list(coords))
        return spans_in_coords

    @staticmethod
    def create_cds(exon_list: List):
        pass

    def _gene_to_string(self, gene_data: Dict, offset: str) -> str:
        gene_id = "gene%s" % str(self.count)
        ss = StringIO()
        ss.write("".join((
            "\t".join((
                gene_data["fasta-id"],
                self.version,
                "gene",
                gene_data["transcripts"][0][0],
                gene_data["transcripts"][0][1],
                ".",
                gene_data["strand"],
                ".",
                "ID=%s" % gene_id
            )),
            "\n",
            "\t".join((
                gene_data["fasta-id"],
                self.version,
                "mRNA",
                gene_data["transcripts"][0][0],
                gene_data["transcripts"][0][1],
                ".",
                gene_data["strand"],
                ".",
                "ID=%s-mRNA" % gene_id
            )),
            "\n",
        )))
        for exon_tuple in gene_data["transcripts"]:
            ss.write("".join((
                "\t".join((
                    gene_data["fasta-id"],
                    self.version,
                    "exon",
                    str(exon_tuple[0]),
                    str(exon_tuple[1]),
                    ".",
                    gene_data["strand"],
                    ".",
                    "ID=%s-exon;Parent=%s" % (gene_id, gene_id)
                )),
                "\n",
                "\t".join((
                    gene_data["fasta-id"],
                    self.version,
                    "CDS",
                    str(exon_tuple[0]),
                    str(exon_tuple[1]),
                    ".",
                    gene_data["strand"],
                    offset,
                    "ID=%s-cds;Parent=%s" % (gene_id, gene_id)
                )),
                "\n",
            )))
        return ss.getvalue()

    def __iter__(self):
        return self.next_gene()


def convert_final_gff3(gff3_file: str, fasta_file: str, filter_function: str, out_file: str):
    gff3 = Gff3Parser(gff3_file, fasta_file, filter_function)
    out_cds = []
    out_prots = []
    gff3_fp = open(out_file + ".gff3", "w")
    for gene, cds, prot in gff3:
        gff3_fp.write(gene)
        if cds is not None:
            out_cds.append(cds)
        if prot is not None:
            out_prots.append(prot)
    SeqIO.write(out_cds, out_file + ".cds.fna", "fasta")
    SeqIO.write(out_prots, out_file + ".faa", "fasta")


def find_orf(record: SeqRecord) -> Tuple[Optional[SeqRecord], int]:
    longest = (0,)
    nuc = str(record.seq)
    for m in re.finditer("ATG", nuc):
        pro = Seq(nuc)[m.start():].translate(to_stop=True)
        if len(pro) > longest[0]:
            longest = (len(pro), m.start(), str(pro))
    if longest[0] > 0:
        return SeqRecord(
            seq=Seq(str(longest[2]) + "*"),
            id=record.id,
            description=record.description,
            name="",
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
        description="GFF3 output final annotations"
    )
    for _file in (ap.args.gff3_file, ap.args.fasta_file):
        assert os.path.exists(_file)
    convert_final_gff3(ap.args.gff3_file, ap.args.fasta_file, ap.args.command, os.path.splitext(ap.args.gff3_file)[0])

if __name__ == "__main__":
    pass
