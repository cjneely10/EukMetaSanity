#!/usr/bin/env python3
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
from io import TextIOWrapper
from Bio.SeqRecord import SeqRecord
from EukMetaSanity.utils.arg_parse import ArgParse


def get_gene(fp: TextIOWrapper):
    line = next(fp).rstrip("\r\n").split("\t")
    exons = []
    _id = line[3]
    contig_id = line[0]
    exons.append((int(line[1]), int(line[2])))
    for line in fp:
        line = line.rstrip("\r\n").split("\t")
        if line[3] != _id:
            yield contig_id, _id, exons
            exons = []
            _id = line[3]
            contig_id = line[0]
        exons.append((int(line[1]), int(line[2])))
    yield contig_id, _id, exons

    # while True:
    #     exons = []
    #     _id = line[3]
    #     contig_id = line[0]
    #     exons.append((int(line[1]), int(line[2])))
    #     while True:
    #         try:
    #             line = next(fp).rstrip("\r\n").split("\t")
    #         except StopIteration:
    #             break
    #         if line[3] != _id:
    #             break
    #         exons.append((int(line[1]), int(line[2])))
    #     yield contig_id, _id, exons


def bed_to_gff3(bed_file: str, fasta_file: str, out_file: str, source: str):
    fp = open(bed_file, "r")
    out_fp = open(out_file, "w")
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    for contig_id, _gene_id, coords in get_gene(fp):
        _min = str(coords[0][0] + 1)
        _max = str(coords[-1][1])
        gene_id = "gene" + _gene_id
        # Build CDS
        seq = StringIO()
        for coord in coords:
            seq.write(str(fasta_dict[contig_id].seq[coord[0]: coord[1]]))
        # Find direction
        _dir = find_orfs(SeqRecord(
            id="1",
            seq=Seq(seq.getvalue()),
        ))
        if _dir is None:
            continue
        # Write gene/mRNA info
        out_fp.write("\t".join((
            contig_id,
            source,
            "gene",
            _min,
            _max,
            ".",
            _dir,
            ".",
            "ID=%s;Name=%s\n" % (gene_id, gene_id),
        )))
        mrna_id = gene_id + "-mRNA"
        out_fp.write("\t".join((
            contig_id,
            source,
            "mRNA",
            _min,
            _max,
            ".",
            _dir,
            ".",
            "ID=%s;Name=%s;Parent=%s\n" % (mrna_id, mrna_id, gene_id),
        )))
        i = 1
        for coord in coords:
            # Write exon info
            out_fp.write("\t".join((
                contig_id,
                source,
                "exon",
                str(coord[0] + 1),
                str(coord[1]),
                ".",
                _dir,
                ".",
                "ID=%s-exon;Name=%s-exon;Parent=%s\n" % (mrna_id + str(i), mrna_id + str(i), mrna_id),
            )))
            # Write CDS info
            out_fp.write("\t".join((
                contig_id,
                source,
                "CDS",
                str(coord[0] + 1),
                str(coord[1]),
                ".",
                _dir,
                "0",
                "ID=%s-CDS;Name=%s-CDS;Parent=%s\n" % (mrna_id + str(i), mrna_id + str(i), mrna_id),
            )))
            i += 1
    out_fp.close()


def find_orfs(record: SeqRecord):
    longest = (0,)
    for _dir, nuc in (("+", str(record.seq)), ("-", str(record.reverse_complement().seq))):
        for m in re.finditer("ATG", nuc):
            pro = Seq(nuc)[m.start():].translate(to_stop=True)
            if len(pro) > longest[0]:
                longest = (len(pro), m.start(), str(pro))
        if longest[0] > 0:
            return _dir


if __name__ == "__main__":
    _ap = ArgParse(
        (
            (("bed_file",),
             {"help": "Path to BED file"}),
            (("fasta_file",),
             {"help": "Input FASTA file"}),
            (("-o", "--output"),
             {"help": "Output path, default stdout", "default": "/dev/stdout"}),
            (("-s", "--source"),
             {"help": "Source, default est", "default": "est"})
        ),
        description="Convert BED to GFF3 format"
    )
    assert os.path.exists(_ap.args.bed_file)
    bed_to_gff3(_ap.args.bed_file, _ap.args.fasta_file, _ap.args.output, _ap.args.source)
