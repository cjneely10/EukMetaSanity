#!/usr/bin/env python3
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from arg_parse import ArgParse
from Bio.SeqRecord import SeqRecord

DEFAULT_OUT = "_____"


def main(ap):
    out_records = []
    record_fp = SeqIO.parse(ap.args.fasta_file, "fasta")
    for record in record_fp:
        # Split to contigs
        seqs = re.split('[nN]+', str(record.seq))
        # Gather contigs
        for i, seq in enumerate(seqs):
            out_records.append(
                SeqRecord(
                    # Mask output if requested
                    seq=(Seq(seq) if ap.args.mask_output else Seq(seq.upper())),
                    id=record.id + "_%i" % (i + 1),
                    description="",
                )
            )
    # Write results
    SeqIO.write(
        out_records,
        (ap.args.output + ".contigs.fasta" if not ap.args.mask_output else ap.args.output + ".masked.contigs.fasta"),
        "fasta"
    )


if __name__ == "__main__":
    _ap = ArgParse(
        (
            (("fasta_file",),
             {"help": "Input FASTA file"}),
            (("-o", "--output"),
             {"help": "Output prefix, default is file without extension", "default": DEFAULT_OUT}),
            (("-m", "--mask_output"),
             {"help": "Keep sequence masking, if present, default False", "default": False, "action": "store_true"}),
        ),
        description="Split each FASTA record in scaffold file into contigs"
    )
    # Confirm path existence
    assert os.path.exists(_ap.args.fasta_file)
    # Set output as file basename if none passed
    if _ap.args.output == DEFAULT_OUT:
        _ap.args.output = os.path.splitext(_ap.args.fasta_file)[0]
    main(_ap)
