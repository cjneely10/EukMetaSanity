#!/usr/bin/env python3
import os
import re
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from typing import List
from Bio.SeqRecord import SeqRecord
from EukMetaSanity.utils.arg_parse import ArgParse


def _parse_args(ap: ArgParse):
    for _file in (ap.args.gff3_file, ap.args.fasta_file):
        assert os.path.exists(_file)


def convert_final_gff3(gff3_file: str, fasta_file: str):
    output_prefix = os.path.splitext(gff3_file)[0]
    gff3_p = open(gff3_file, "r")
    out_p = open(output_prefix + ".nr.gff3", "w")
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    # Skip over header lines
    for _ in range(2):
        next(gff3_p)
    cds_list = []
    i = 1
    for line in gff3_p:
        if line.startswith("#"):
            out_p.write(line)
            continue
        elif line.startswith(">"):
            break
        line = line.split()
        if line[2] in ("gene", "transcript"):
            seq = StringIO()
            gene_id = "gene%i" % i
            mrna_id = gene_id + "-mRNA"
            out_p.write("\t".join((
                *line[:2],
                "gene",
                *line[3:-1],
                "ID=%s;Name=%s\n" % (gene_id, gene_id),
            )))
            out_p.write("\t".join((
                *line[:2],
                "mRNA",
                *line[3:-1],
                "ID=%s;Name=%s;Parent=%s\n" % (mrna_id, mrna_id, gene_id),
            )))
            strand = line[6]
            line = next(gff3_p).split()
            j = 1
            while line[2] not in ("gene", "transcript"):
                out_p.write("\t".join((
                    *line[:2],
                    "exon",
                    *line[3:-2],
                    ".",
                    "Parent=%s;ID=%s\n" % (mrna_id, mrna_id + "-exon%i" % j),
                )))
                out_p.write("\t".join((
                    *line[:2],
                    "CDS",
                    *line[3:-2],
                    "0",
                    "Parent=%s;ID=%s\n" % (mrna_id, mrna_id + "-cds%i" % j),
                )))
                start, end = int(line[3]), int(line[4])
                if strand == "+":
                    seq.write(str(fasta_dict[line[0]].seq[start - 1:end]))
                else:
                    seq.write(str(fasta_dict[line[0]].reverse_complement().seq[start - 1:end]))
                try:
                    line = next(gff3_p).split()
                except StopIteration:
                    break
                j += 1
            i += 1
            cds_list.append(
                SeqRecord(
                    id=gene_id,
                    seq=Seq(seq.getvalue()),
                    description="strand=%s" % strand,
                    name="",
                )
            )
    out_p.close()
    SeqIO.write(cds_list, output_prefix + ".cds.fna", "fasta")
    SeqIO.write(find_orfs(cds_list), output_prefix + ".faa", "fasta")


def find_orfs(cds_list: List[SeqRecord]):
    for record in cds_list:
        longest = (0,)
        if record.description[-1] == "+":
            nuc = str(record.seq)
        else:
            nuc = str(record.reverse_complement().seq)
        for m in re.finditer("ATG", nuc):
            pro = Seq(nuc)[m.start():].translate(to_stop=True)
            if len(pro) > longest[0]:
                longest = (len(pro), m.start(), str(pro))
        if longest[0] >= 30:
            yield SeqRecord(
                seq=Seq(str(longest[2]) + "*"),
                id=record.id,
                description=record.description
            )


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

