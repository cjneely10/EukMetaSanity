#!/usr/bin/env python3
import os
import re
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from typing import Optional, Tuple
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
    prot_list = []
    cds_ids = set()
    i = 1
    line = next(gff3_p).split()
    end = 1
    at_end = False
    while True:
        if line[0].startswith("#"):
            out_p.write("\t".join(line) + "\n")
            line = next(gff3_p).split()
            continue
        elif line[0].startswith(">"):
            break
        if line[2] in ("gene", "transcript"):
            _new_end = int(line[4])
            _new_start = int(line[3])
            if _new_start <= end:
                line = next(gff3_p).split()
                while line[2] not in ("gene", "transcript"):
                    try:
                        line = next(gff3_p).split()
                    except StopIteration:
                        break
                continue
            end = _new_end
            # Gather CDS sequence
            seq = StringIO()
            # Generate gene/mRNA ids
            gene_id = "gene%i" % i
            mrna_id = gene_id + "-mRNA"
            # Write locus info
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
            # Start reading in exon/CDS info
            strand = line[6]
            line = next(gff3_p).split()
            j = 1
            # Gather lines and CDS sequence region
            out_lines = []
            while line[2] not in ("gene", "transcript"):
                start, end = int(line[3]), int(line[4])
                if strand == "+":
                    seq.write(str(fasta_dict[line[0]].seq[start - 1:end]))
                else:
                    seq.write(str(fasta_dict[line[0]].reverse_complement().seq[start - 1:end]))
                out_lines.append(line)
                try:
                    line = next(gff3_p).split()
                except StopIteration:
                    at_end = True
                    break
            # Ensure CDS id is unique
            while gene_id in cds_ids:
                if "_" in gene_id:
                    _num = gene_id.split("_")[1]
                    gene_id = gene_id[:-1 * (len(_num) + 1)] + "_" + str(int(_num) + 1)
                else:
                    gene_id = gene_id + "_" + "1"
            # Add new id
            cds_ids.add(gene_id)
            # Create CDS record and save
            record = SeqRecord(
                id=gene_id,
                seq=Seq(seq.getvalue()),
                description="strand=%s" % strand,
                name="",
            )
            cds_list.append(record)
            # Read to find longest ORF
            # Translate and return offset
            prot_record, offset = find_orf(record)
            if prot_record is not None:
                prot_list.append(prot_record)
            # Write lines with offset info
            for _line in out_lines:
                out_p.write("\t".join((
                    *_line[:2],
                    "exon",
                    *_line[3:-2],
                    ".",
                    "Parent=%s;ID=%s\n" % (mrna_id, mrna_id + "-exon%i" % j),
                )))
                out_p.write("\t".join((
                    *_line[:2],
                    "CDS",
                    *_line[3:-2],
                    str(offset % 3),
                    "Parent=%s;ID=%s\n" % (mrna_id, mrna_id + "-cds%i" % j),
                )))
                j += 1
            i += 1
            if at_end:
                break

    out_p.close()
    # Write CDS/protein info
    SeqIO.write(cds_list, output_prefix + ".cds.fna", "fasta")
    SeqIO.write(prot_list, output_prefix + ".faa", "fasta")


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
