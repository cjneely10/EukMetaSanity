#!/usr/bin/env python3
import os
import regex as re
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from operator import itemgetter
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from EukMetaSanity.utils.arg_parse import ArgParse
from typing import List, Dict, Tuple, Generator, Optional


class Gene:
    def __init__(self, ab_initio_data: List):
        self.exons: List = ab_initio_data

    def add_evidence(self, evidence_data: List):
        if len(self.exons) == 0:
            return
        if len(evidence_data) == 0:
            return
        out_exons = self.exons
        # out_exons = [self.exons[0]]
        # if len(self.exons) > 1:
        #     out_exons.append(self.exons[-1])
        # for ab_exon in self.exons[1:-1]:
        #     is_found = False
        #     for exon in evidence_data:
        #         if Gene.in_exon(exon, ab_exon):
        #             is_found = True
        #             break
        #     if is_found:
        #         out_exons.append(ab_exon)
        for exon in evidence_data:
            is_found = False
            for ab_exon in out_exons:
                if Gene.in_exon(ab_exon, exon):
                    is_found = True
                    break
            if not is_found:
                out_exons.append(exon)
        self.exons = out_exons
        self.exons.sort(key=itemgetter(0))

    # Returns if part of query coord overlaps target coord
    @staticmethod
    def in_exon(query_coord: Tuple[int, int, int], target_coord: Tuple[int, int, int]) -> bool:
        return len(set(range(*target_coord[:-1])).intersection(set(range(*query_coord[:-1])))) > 0


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
            self._count += 1
            # Putative gene
            gene_data = {
                "geneid": "gene%i" % self._count,
                "fasta-id": line[0],
                "strand": line[6],
            }
            transcripts = []
            line = next(self.fp).rstrip("\r\n").split("\t")
            while line[2] != "locus":
                # Read in transcript info
                transcripts.append(
                    [line[1], []]  # First line is a transcript: source,tstart,tend
                )
                line = next(self.fp).rstrip("\r\n").split("\t")
                # Add exon to current info
                while line[2] not in ("transcript", "locus", "gene"):
                    if line[2] == "CDS":
                        transcripts[-1][-1].append(
                            (int(line[3]), int(line[4]), int(line[7]))  # exstart,exend,offset
                        )
                    line = next(self.fp).rstrip("\r\n").split("\t")
            # Merge based on name
            data = defaultdict(list)
            for transcript in transcripts:
                data[transcript[0]].extend(transcript[-1])
            for transcript in data.keys():
                data[transcript].sort(key=itemgetter(0))
            # Store in data and yield
            gene_data["transcripts"] = data
            yield gene_data


class GffMerge:
    def __init__(self, gff3_path: str, fasta_path: str):
        self.reader = GffReader(gff3_path)
        self.fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

    def merge(self) -> Generator[Tuple[Dict, SeqRecord, SeqRecord, List[int]], None, None]:
        gene_data: Dict
        for gene_data in self.reader:
            # Generate initial exon structure
            gene = Gene(gene_data["transcripts"]["ab-initio"])
            # Keep exons with evidence, add exons missed by ab-initio
            gene.add_evidence(gene_data["transcripts"]["metaeuk"])
            gene_data["transcripts"] = gene.exons
            # Return data to write and output FASTA records
            yield (gene_data, *self.create_cds(gene_data))

    def create_cds(self, gene_data: dict) -> Tuple[Optional[SeqRecord], Optional[SeqRecord], List[int]]:
        orig_seq = str(self.fasta_dict[gene_data["fasta-id"]].seq)
        strand = gene_data["strand"]
        out_cds: List[SeqRecord] = []
        offsets: List[int] = []
        seq = StringIO()
        for i, exon in enumerate(gene_data["transcripts"]):
            start, end = 0, 0
            if strand == "-":
                end = exon[2]
            else:
                start = exon[2]
            dist = exon[1] - end - exon[0] - start + 1
            if dist % 3 == 0:
                dist = 3
            if strand == "+":
                seq.write(orig_seq[exon[0] - 1 + start: exon[1] - end] + "N" * (3 - dist % 3))
            else:
                seq.write("N" * (3 - dist % 3) + orig_seq[exon[0] - 1 + start: exon[1] - end])
            record = SeqRecord(seq=Seq(seq.getvalue()))
            if strand == "-":
                record = record.reverse_complement()
            out_cds.append(record)
            offsets.append(exon[2])
        # out = GffMerge.find_orf(record, exon[2])
        # if out is not None:
        #     out_prots.append(out[0])
        #     out_cds.append(out[1])
        #     offsets.append(out[2])
        # else:
        #     offsets.append(-1)
        cds = Seq("".join(str(val.seq) for val in out_cds))
        return (
            SeqRecord(
                id=gene_data["geneid"],
                name="",
                description="strand=%s" % strand,
                seq=Seq(str(cds.translate()).replace("X", ""))
            ),
            SeqRecord(
                id=gene_data["geneid"],
                name="",
                description="strand=%s" % strand,
                seq=cds
            ),
            offsets
        )

    @staticmethod
    def find_orf(record: SeqRecord, offset: int):
        nuc = str(record.seq)
        return (
            SeqRecord(
                seq=Seq(nuc[offset:]).translate(),
                id=record.id,
                description=record.description,
                name=""
            )
        )

    # @staticmethod
    # def find_orf(record: SeqRecord) -> Optional[Tuple[SeqRecord, SeqRecord]]:
    #     longest = (0,)
    #     nuc = str(record.seq)
    #     for m in re.finditer("ATG", nuc):
    #         pro = Seq(nuc[m.start():]).translate(to_stop=True)
    #         if len(pro) > longest[0]:
    #             longest = (len(pro), pro, nuc[m.start():m.start() + len(pro) * 3 + 3])
    #     if longest[0] > 0:
    #         return (
    #             SeqRecord(
    #                 seq=longest[1],
    #                 id=record.id,
    #                 description=record.description,
    #                 name="",
    #             ),
    #             SeqRecord(
    #                 seq=Seq(longest[2]),
    #                 id=record.id,
    #                 description=record.description,
    #                 name="",
    #             )
    #         )
    #     return None

    @staticmethod
    def pad_seq(sequence):
        remainder = len(sequence) % 3
        return sequence if remainder == 0 else sequence + Seq('N' * (3 - remainder))


class GffWriter:
    def __init__(self, in_gff3_path: str, fasta_file: str):
        self.in_fp = open(in_gff3_path, "r")
        self.base = os.path.splitext(in_gff3_path)[0]
        self.out_fp = open(self.base + ".nr.gff3", "w")
        self.merger = GffMerge(in_gff3_path, fasta_file)

    def write(self):
        self.out_fp.write("# EukMetaSanity merged annotations\n")
        out_prots: List[SeqRecord] = []
        out_cds: List[SeqRecord] = []
        current_id = ""
        for gene_dict, prot, cds, offsets in self.merger.merge():
            if gene_dict["fasta-id"] != current_id:
                current_id = gene_dict["fasta-id"]
                self.out_fp.write("# Region %s\n" % current_id)
            gene = GffWriter._gene_dict_to_string(gene_dict, offsets)
            if gene is not None:
                self.out_fp.write(gene)
            if prot and len(prot.seq) > 0:
                out_prots.append(prot)
            if cds and len(cds.seq) > 0:
                out_cds.append(cds)
        SeqIO.write(out_prots, self.base + ".faa", "fasta")
        SeqIO.write(out_cds, self.base + ".cds.fna", "fasta")

    @staticmethod
    def _gene_dict_to_string(gene_data: Dict, offsets: List[int]) -> Optional[str]:
        gene_id = gene_data["geneid"]
        mrna_id = gene_id + "-mRNA"
        version = "EukMS"
        ss = StringIO()
        if len(gene_data["transcripts"]) == 0:
            return
        ss.write("".join((
            "\t".join((
                gene_data["fasta-id"], version,
                "gene", str(gene_data["transcripts"][0][0]), str(gene_data["transcripts"][-1][1]),
                ".", gene_data["strand"], ".", "ID=%s" % gene_id
            )),
            "\n",
            "\t".join((
                gene_data["fasta-id"], version,
                "mRNA", str(gene_data["transcripts"][0][0]), str(gene_data["transcripts"][-1][1]),
                ".", gene_data["strand"], ".", "ID=%s-mRNA;Parent=%s" % (mrna_id, gene_id)
            )),
            "\n",
        )))
        for j, exon_tuple in enumerate(gene_data["transcripts"], start=1):
            ss.write("".join((
                "\t".join((
                    gene_data["fasta-id"], version,
                    "exon", str(exon_tuple[0]), str(exon_tuple[1]),
                    ".", gene_data["strand"], ".", "ID=%s-exon%i;Parent=%s" % (mrna_id, j, mrna_id)
                )),
                "\n",
                "\t".join((
                    gene_data["fasta-id"], version,
                    "CDS", str(exon_tuple[0]), str(exon_tuple[1]),
                    ".", gene_data["strand"], str(offsets[j - 1]), "ID=%s-cds%i;Parent=%s" % (mrna_id, j, mrna_id)
                )),
                "\n",
            )))
        return ss.getvalue()


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
    writer = GffWriter(ap.args.gff3_file, ap.args.fasta_file)
    writer.write()
