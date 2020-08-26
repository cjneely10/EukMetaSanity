#!/usr/bin/env python3
import os
import itertools
import regex as re
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
from operator import itemgetter
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from EukMetaSanity.utils.arg_parse import ArgParse
from typing import Optional, Tuple, Generator, List, Dict


class Gff3Parser:
    def __init__(self, gff3_file: str, fasta_file: str, priority: str = "metaeuk"):
        assert priority in dir(Gff3Parser)
        self.fp = open(gff3_file, "r")
        self.priority = getattr(Gff3Parser, priority, lambda a: a)
        self.version = priority
        self.fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        self.count = 0

    def __iter__(self):
        return self.next_gene()

    def next_gene(self) -> Generator[Tuple[str, Optional[SeqRecord], Optional[SeqRecord]], None, None]:
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
                        [line[1], line[3], line[4], []]  # First line is a transcript: source,tstart,tend
                    )
                    line = next(self.fp).rstrip("\r\n").split("\t")
                    # Add exon to current info
                    while line[2] not in ("transcript", "locus", "gene"):
                        if line[2] == "CDS":
                            transcripts[-1][-1].append(
                                (int(line[3]), int(line[4]))  # source,exstart,exend
                            )
                        line = next(self.fp).rstrip("\r\n").split("\t")
                # Merge based on name
                data = defaultdict(list)
                for transcript in transcripts:
                    data[transcript[0]].extend(transcript[-1])
                # Filter for specific transcripts
                gene_data["transcripts"] = [data[_t[0]] for _t in transcripts]
                record, orf, offset = self.find_longest_orf(gene_data)
                yield (
                    self._gene_to_string(gene_data, offset),
                    record,
                    orf,
                )

    def create_cds(self, gene_data: Dict) -> SeqRecord:
        seq = StringIO()
        orig_seq = self.fasta_dict[gene_data["fasta-id"]]
        orig_seq = str(orig_seq.seq)
        for transcript in gene_data["transcripts"]:
            seq.write(orig_seq[int(transcript[0]) - 1: int(transcript[1])])
        seq = Seq(seq.getvalue())
        if gene_data["strand"] == "-":
            seq = seq.reverse_complement()
        return SeqRecord(
            seq=seq,
            id="gene%s" % str(self.count),
            description="strand=%s" % gene_data["strand"],
            name="",
        )

    def find_longest_orf(self, gene_data: Dict) -> Tuple[SeqRecord, SeqRecord, int]:
        max_idx: int = 0
        max_len: int = 0
        max_offset: int = 0
        max_cds: Optional[SeqRecord] = None
        max_prot: Optional[SeqRecord] = None
        search_data = [*gene_data["transcripts"], self.priority(gene_data["transcripts"])]
        for i, transcript in enumerate(search_data):
            seq = StringIO()
            orig_seq = self.fasta_dict[gene_data["fasta-id"]]
            orig_seq = str(orig_seq.seq)
            for transcr in transcript:
                seq.write(orig_seq[int(transcr[0]) - 1: int(transcr[1])])
            seq = Seq(seq.getvalue())
            if gene_data["strand"] == "-":
                seq = seq.reverse_complement()
            cds = SeqRecord(
                seq=seq,
                id="gene%s" % str(self.count),
                description="strand=%s" % gene_data["strand"],
                name="",
            )
            out = Gff3Parser.find_orf(cds)
            if out is not None:
                if max_len < len(out[0].seq):
                    max_idx = i
                    max_len = len(out[0].seq)
                    max_cds = cds
                    max_prot = out[0]
                    max_offset = out[1]
        gene_data["transcripts"] = search_data[max_idx]
        return max_cds, max_prot, max_offset

    def _gene_to_string(self, gene_data: Dict, offset: int) -> Optional[str]:
        gene_id = "gene%s" % str(self.count)
        ss = StringIO()
        if len(gene_data["transcripts"]) == 0:
            return
        ss.write("".join((
            "\t".join((
                gene_data["fasta-id"], self.version,
                "gene", gene_data["start"], gene_data["end"],
                ".", gene_data["strand"], ".", "ID=%s" % gene_id
            )),
            "\n",
            "\t".join((
                gene_data["fasta-id"], self.version,
                "mRNA", gene_data["start"], gene_data["end"],
                ".", gene_data["strand"], ".", "ID=%s-mRNA;Parent=%s" % (gene_id, gene_id)
            )),
            "\n",
        )))
        for j, exon_tuple in enumerate(gene_data["transcripts"], start=1):
            ss.write("".join((
                "\t".join((
                    gene_data["fasta-id"], self.version,
                    "exon", str(exon_tuple[0]), str(exon_tuple[1]),
                    ".", gene_data["strand"], ".", "ID=%s-exon%i;Parent=%s" % (gene_id, j, gene_id)
                )),
                "\n",
                "\t".join((
                    gene_data["fasta-id"], self.version,
                    "CDS", str(exon_tuple[0]), str(exon_tuple[1]),
                    ".", gene_data["strand"], str(offset), "ID=%s-cds%i;Parent=%s" % (gene_id, j, gene_id)
                )),
                "\n",
            )))
        return ss.getvalue()

    @staticmethod
    def metaeuk(line: List[List]) -> List[List]:
        return Gff3Parser.filter_specific(line, "metaeuk")

    @staticmethod
    def abinitio(line: List[List]) -> List[List]:
        return Gff3Parser.filter_specific(line, "ab-initio")

    @staticmethod
    def filter_specific(line: List[List], name: str) -> List[List]:
        for _l in line:
            if _l[0] == name or len(_l[-1]) < 2:
                return _l[-1]
        return []

    @staticmethod
    def merge(line: List[List]) -> List[Tuple]:
        # Sort coordinates by start value
        ranges_in_coords = sorted(
            [
                (_v[0], _v[1])
                for _v in itertools.chain(*line)
            ],
            key=itemgetter(0)
        )
        # Will group together matching sections into spans
        spans_in_coords = [list(ranges_in_coords[0]), ]
        for coords in ranges_in_coords[1:]:
            # The start value is within the farthest range of current span and the end value extends past current
            if coords[0] <= spans_in_coords[-1][1] < coords[1]:
                spans_in_coords[-1][1] = coords[1]
            # The start value is past the range of the current span, append old span to list to return
            elif coords[0] > spans_in_coords[-1][1]:
                spans_in_coords.append(list(coords))
        return [tuple(val) for val in spans_in_coords]

    @staticmethod
    def find_orf(record: SeqRecord) -> Optional[Tuple[SeqRecord, int]]:
        longest = (0,)
        nuc = str(record.seq)
        start_pos = re.compile("ATG")
        for m in start_pos.finditer(nuc, overlapped=True):
            pro = Seq(nuc[m.start():]).translate(to_stop=True)
            if len(pro) > longest[0]:
                longest = (len(pro), m.start(), pro)
        if longest[0] > 0:
            return (
                SeqRecord(
                    seq=longest[2],
                    id=record.id,
                    description=record.description,
                    name="",
                ),
                longest[1] % 3
            )
        return None


def convert_final_gff3(gff3_file: str, fasta_file: str, out_prefix: str):
    gff3 = Gff3Parser(gff3_file, fasta_file, "merge")
    out_cds = []
    out_prots = []
    with open(out_prefix + ".nr.gff3", "w") as gff3_fp:
        for gene, cds, prot in gff3:
            if gene is not None:
                gff3_fp.write(gene)
            if cds is not None:
                out_cds.append(cds)
            if prot is not None:
                out_prots.append(prot)
    SeqIO.write(out_cds, out_prefix + ".cds.fna", "fasta")
    SeqIO.write(out_prots, out_prefix + ".faa", "fasta")


if __name__ == "__main__":
    ap = ArgParse(
        (
            (("-g", "--gff3_file"),
             {"help": ".tmp.nr.gff3 file", "required": True}),
            (("-f", "--fasta_file"),
             {"help": "FASTA file", "required": True}),
        ),
        description="GFF3 output final annotations"
    )
    for _file in (ap.args.gff3_file, ap.args.fasta_file):
        assert os.path.exists(_file), _file
    convert_final_gff3(ap.args.gff3_file, ap.args.fasta_file, os.path.splitext(ap.args.gff3_file)[0])
