#!/usr/bin/env python3
import os
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
from random import random
from arg_parse import ArgParse
from Bio.SeqRecord import SeqRecord


# Calculate N50 for file
def n50(records):
    total = 0
    total_len = sum(len(rec.seq) for rec in records) // 2
    records.sort(key=lambda rec: len(rec.seq), reverse=True)
    for _i, record in enumerate(records):
        total += len(record.seq)
        if total > total_len:
            return len(records[_i].seq)


# Generate random numbers, but prefer higher numbers
def monte_carlo():
    r1 = random()
    while True:
        if random() < r1:
            return r1


# Subset each record in FASTA file using monte_carlo method
def get_subset(out_records, record, size, len_record, skip_size):
    start_pos = 0
    while start_pos < len_record:
        start_pos += int(skip_size * random())
        contig_len = int(size * monte_carlo())
        if start_pos + contig_len > len_record:
            contig_len = len_record - start_pos - 1
        # Store in output list
        out_records.append(
            SeqRecord(
                id=record.id + ".%s_%s" % (str(start_pos + 1), str(start_pos + contig_len + 1)),
                seq=Seq(str(record.seq[start_pos:start_pos + contig_len])),
                description=""
            )
        )
        start_pos += contig_len + 1


def main(ap):
    record_fp = SeqIO.parse(ap.args.fasta_file, "fasta")
    out_fp = open(os.path.join(ap.args.output_directory, "contig-stats.tsv"), "w")
    out_fp.write("".join([
        "\t".join(["ID", "N50", "Num Contigs", "Length ratio"]),
        "\n"
    ]))
    for record in record_fp:
        len_record = len(record.seq)
        for size in ap.args.sizes:
            # Do not subset if sequence is too small
            if size > len_record:
                continue
            out_records = []
            get_subset(out_records, record, size, len_record, ap.args.skip_size)
            out_fp.write("".join([
                "\t".join([
                    record.id,
                    str(n50(out_records)),
                    str(len(out_records)),
                    str(sum([len(rec.seq) for rec in out_records]) / len_record),
                ]),
                "\n"
            ]))
            SeqIO.write(
                out_records,
                os.path.join(
                    ap.args.output_directory,
                    os.path.basename(str(record.id) + ".%i" % size + ".fasta")
                ),
                "fasta"
            )


if __name__ == "__main__":
    _ap = ArgParse(
        (
            (("fasta_file",),
             {"help": "Input FASTA sequence"}),
            (("-s", "--sizes"),
             {"help": "Comma-separated list of contig lengths to simulate", "default": "100000"}),
            (("-o", "--output_directory"),
             {"help": "Output directory, default current working directory", "default": os.getcwd()}),
            (("-k", "--skip_size"),
             {"help": "Skip max bp between contigs, default 1000", "default": "1000"}),
        ),
        description="Subset a FASTA file"
    )
    # Check argument existence
    _ap.args.fasta_file = str(Path(_ap.args.fasta_file).resolve())
    assert os.path.exists(_ap.args.fasta_file)
    try:
        # Confirm passed lists consist of integer values
        _ap.args.sizes = _ap.args.sizes.split(",")
        for i, _arg in enumerate(_ap.args.sizes):
            _ap.args.sizes[i] = int(float(_arg))
        # Confirm numeric types
        _ap.args.skip_size = int(float(_ap.args.skip_size))
        # Make output directory
        if not os.path.exists(_ap.args.output_directory):
            os.makedirs(_ap.args.output_directory)
    except ValueError as e:
        print(e)
        exit(1)
    main(_ap)
