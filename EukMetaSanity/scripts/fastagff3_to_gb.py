#!/usr/bin/env python3
import os
import sys
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import UnknownSeq
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation
from EukMetaSanity.utils.arg_parse import ArgParse
"""
https://github.com/chapmanb/bcbb/blob/master/gff/Scripts/gff/gff_to_genbank.py

"""


def _parse_args(ap):
    assert os.path.exists(ap.args.fasta_file)
    assert os.path.exists(ap.args.gff3_file)
    assert ap.args.format in ("genbank", "CDS")


def _fix_ncbi_id(fasta_iter):
    """GenBank identifiers can only be 16 characters; try to shorten NCBI.
    """
    i = 1
    for rec in fasta_iter:
        if len(rec.name) > 16 and rec.name.find("|") > 0:
            new_id = [x for x in rec.name.split("|") if x][-1]
            sys.stderr.write("%s\t%s\n" % (rec.id, new_id))
            rec.id = new_id
            rec.name = new_id
        elif len(rec.name) > 16:
            new_id = rec.id[:15 - len(str(i))] + "_" + str(i)
            sys.stderr.write("%s\t%s\n" % (rec.id, new_id))
            i += 1
            rec.id = new_id
            rec.name = new_id
        yield rec


def _check_gff(gff_iterator, _type=None):
    """Check GFF files before feeding to SeqIO to be sure they have sequences.
    """
    for rec in gff_iterator:
        if isinstance(rec.seq, UnknownSeq):
            print("Warning: FASTA sequence not found for '%s' in GFF file" % (
                rec.id))
            rec.seq.alphabet = generic_dna
        yield _flatten_features(rec)


def _flatten_features(rec):
    """Make sub_features in an input rec flat for output.
    GenBank does not handle nested features, so we want to make
    everything top level.
    """
    out = []
    for f in rec.features:
        cur = [f]
        while len(cur) > 0:
            nextf = []
            for curf in cur:
                out.append(curf)
                if len(curf.sub_features) > 0:
                    nextf.extend(curf.sub_features)
            cur = nextf
    for feat in out:
        if feat.type == "transcript":
            feat.type = "mRNA"
    rec.features = [
        SeqFeature(
            FeatureLocation(0, len(rec.seq)), type="source",
        ),
        *out
    ]
    return rec


def write_genbank(fasta_file, gff3_file, output_file):
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(fasta_file)
    if not os.path.exists(gff3_file):
        raise FileNotFoundError(gff3_file)
    SeqIO.write(
        _check_gff(
            _fix_ncbi_id(
                GFF.parse(gff3_file, SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", generic_dna))),
            ),
        ),
        output_file,
        "genbank",
    )


def write_cds(genbank_file, output_file):
    out = []
    for seq_record in SeqIO.parse(genbank_file, "genbank"):
        for seq_feature in seq_record.features:
            if seq_feature.type == "CDS":
                out.append(seq_feature.translate(seq_record.seq, cds=False))
    SeqIO.write(out, output_file, "fasta")


def main(ap):
    if ap.args.format == "genbank":
        write_genbank(ap.args.fasta_file, ap.args.gff3_file, ap.args.output)
    elif ap.args.format == "CDS":
        write_cds(ap.args.fasta_file, ap.args.output)


if __name__ == "__main__":
    _ap = ArgParse(
        (
            (("fasta_file",),
             {"help": "FASTA file"}),
            (("gff3_file",),
             {"help": "GFF3 mapping file"}),
            (("-o", "--output"),
             {"help": "Output path, default stdout", "default": "/dev/stdout"}),
            (("-f", "--format"),
             {"help": "Select from genbank/CDS, default genbank", "default": "genbank"})
        ),
        description="Convert FASTA and GFF3 file to Genbank format"
    )
    _parse_args(_ap)
    main(_ap)
