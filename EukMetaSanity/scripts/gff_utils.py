import os
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from EukMetaSanity.utils.arg_parse import ArgParse
from typing import List, Dict, Tuple, Generator, Optional


class Gene:
    def __init__(self, ab_initio_data: List):
        self._exons: List = ab_initio_data[-1]
        self.source, self.start, self.end = ab_initio_data[:-1]

    def add_evidence(self, evidence_data: List):
        for ab_exon in self._exons:
            is_found = False
            for exon in evidence_data[-1]:
                if Gene.in_exon(ab_exon, exon):
                    is_found = True
                    break
            if not is_found:
                del ab_exon
        to_add = []
        for exon in evidence_data[-1]:
            for ab_exon in self._exons:
                if not Gene.in_exon(exon, ab_exon):
                    to_add.append(exon)
        self._exons.extend(to_add)

    @property
    def exons(self) -> List:
        return self._exons

    # Returns if part of query coord overlaps target coord
    @staticmethod
    def in_exon(query_coord: Tuple[int, int], target_coord: Tuple[int, int]):
        return query_coord[0] <= target_coord[0] <= query_coord[1] or \
               query_coord[0] <= target_coord[1] <= query_coord[1]


class GffReader:
    def __init__(self, gff3_path: str):
        self.fp = open(gff3_path, "r")
        self._count = 0

    def __iter__(self):
        return self.next_gene()

    def next_gene(self) -> Generator[defaultdict, None, None]:
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
                self._count += 1
                # Putative gene
                gene_data = {
                    "geneid": self._count,
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
                                (int(line[3]), int(line[4]))  # exstart,exend,offset
                            )
                        line = next(self.fp).rstrip("\r\n").split("\t")
                # Merge based on name
                data = defaultdict(list)
                for transcript in transcripts:
                    data[transcript[0]].extend(transcript[-1])
                # Store in data and yield
                gene_data["transcripts"] = data
                yield gene_data


class GffMerge:
    def __init__(self, gff3_path: str, fasta_path: str):
        self.reader = GffReader(gff3_path)
        self.fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    def merge(self):
        gene_data: Dict
        for gene_data in self.reader:
            # Generate initial exon structure
            gene = Gene(gene_data["transcripts"]["ab-initio"])
            # Keep exons with evidence, add exons missed by ab-initio
            gene.add_evidence(gene_data["transcripts"]["metaeuk"])
            gene_data["transcripts"] = gene.exons
            # Return data to write and output FASTA records
            yield (
                gene_data,
                *self.create_cds(gene_data)
            )

    def create_cds(self, gene_data: dict) -> Tuple[SeqRecord, SeqRecord]:
        orig_seq = self.fasta_dict[gene_data["fasta-id"]]
        strand = self.fasta_dict["strand"]
        for exon in gene_data["transcripts"]:
            seq = StringIO()
            seq.write(orig_seq[exon[0] - 1: exon[0]])


class GffWriter:
    def __init__(self, in_gff3_path: str, out_gff3_path: str, fasta_file: str):
        self.fasta_path = fasta_file
        self.in_fp = open(in_gff3_path, "r")
        self.out_fp = open(out_gff3_path, "w")
        self.merger = GffMerge(in_gff3_path, fasta_file)

    def write(self):
        pass

    @staticmethod
    def _gene_dict_to_string(gene_data: Dict):
        gene_id = gene_data["geneid"]
        version = "EukMS"
        ss = StringIO()
        if len(gene_data["transcripts"]) == 0:
            return
        ss.write("".join((
            "\t".join((
                gene_data["fasta-id"], version,
                "gene", gene_data["start"], gene_data["end"],
                ".", gene_data["strand"], ".", "ID=%s" % gene_id
            )),
            "\n",
            "\t".join((
                gene_data["fasta-id"], version,
                "mRNA", gene_data["start"], gene_data["end"],
                ".", gene_data["strand"], ".", "ID=%s-mRNA;Parent=%s" % (gene_id, gene_id)
            )),
            "\n",
        )))
        for j, exon_tuple in enumerate(gene_data["transcripts"], start=1):
            ss.write("".join((
                "\t".join((
                    gene_data["fasta-id"], version,
                    "exon", str(exon_tuple[0]), str(exon_tuple[1]),
                    ".", gene_data["strand"], ".", "ID=%s-exon%i;Parent=%s" % (gene_id, j, gene_id)
                )),
                "\n",
                "\t".join((
                    gene_data["fasta-id"], version,
                    "CDS", str(exon_tuple[0]), str(exon_tuple[1]),
                    ".", gene_data["strand"], str(exon_tuple[2]), "ID=%s-cds%i;Parent=%s" % (gene_id, j, gene_id)
                )),
                "\n",
            )))
        return ss.getvalue()


if __name__ == "__main__":
    pass
