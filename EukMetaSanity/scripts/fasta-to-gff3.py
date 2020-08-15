#!/usr/bin/env python3
import os
from BCBio import GFF
from Bio import SeqIO
from decimal import Decimal
from collections import namedtuple
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from EukMetaSanity.utils.arg_parse import ArgParse
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Result as read in from DNA/contig side as subject
Result = namedtuple("Result", ("loc_type", "sstart", "send", "strand", "score"))


def metaeuk(metaeuk_file_path, data, *args, **kwargs):
    records_fp = SeqIO.parse(metaeuk_file_path, "fasta")
    for record in records_fp:
        # Metaeuk header contains pipes as delimiters - be sure to remove all in input seqs!
        line = record.id.split("|")
        # Retain gene info in nested lists
        recs = []
        # Determine strand
        if line[2] == "+":
            strand = 1
        elif line[2] == "-":
            strand = -1
        else:
            strand = 0
        # Add base record as first record in nested list
        recs.append(
            Result(
                loc_type="gene",
                sstart=int(line[6]),
                send=int(line[7]),
                strand=strand,
                score=line[3],
            )
        )
        # Store exon/CDS info from rest of header at nested list loc
        for coords in line[8:]:
            start, end, length = coords.split(":")
            if strand < 0:
                start, end = end, start
            for _type in ("CDS",):
                recs.append(
                    Result(
                        loc_type=_type,
                        sstart=int(start.split("[")[0]),
                        send=int(end.split("[")[0]),
                        strand=strand,
                        score=".",
                    )
                )
        # Add nested list
        data[line[1]].append(recs)


def _parse_metaeuk(i, ap, feature_data, record, rec):
    for features in feature_data.get(record.id, [[]]):
        if len(features) == 0:
            continue
        # First value is gene
        rec.features.append(
            SeqFeature(
                FeatureLocation(features[0].sstart, features[0].send), type=features[0].loc_type,
                strand=features[0].strand,
                qualifiers={"score": features[0].score, "source": "metaeuk", "ID": "gene%i" % i}
            )
        )
        # Remaining in list are features that describe the gene
        rec.features[-1].sub_features = []
        for feature in features[1:]:
            rec.features[-1].sub_features.append(
                SeqFeature(
                    FeatureLocation(feature.sstart, feature.send), type=feature.loc_type,
                    strand=features[0].strand,
                    qualifiers={"score": feature.score, "source": "metaeuk"}
                )
            )
        i += 1
    return i


# # Write functions
def _iter_gff(fasta_file, feature_data, ap):
    assert isinstance(feature_data, defaultdict)
    fasta_fp = SeqIO.parse(fasta_file, "fasta")
    i = 1
    for record in fasta_fp:
        # Each GFF entry based on FASTA record
        rec = SeqRecord(id=record.id, seq=record.seq)
        feature_list = feature_data.get(record.id, [])
        if len(feature_list) > 0:
            i = _parse_metaeuk(i, ap, feature_data, record, rec)
        yield rec


# # Driver logic
def _parse_args(ap):
    # Confirm path existence
    assert os.path.exists(ap.args.data_file)
    assert os.path.exists(ap.args.fasta_file)


def _main(ap):
    # Store results data
    data = defaultdict(list)
    # Call parsing function - ensured to exist
    metaeuk(open(ap.args.data_file), data, ap)
    # Write results
    GFF.write(_iter_gff(ap.args.fasta_file, data, ap), open(ap.args.output, "w"))


if __name__ == "__main__":
    _ap = ArgParse(
        (
            (("fasta_file",),
             {"help": "FASTA file on which to map results"}),
            (("data_file",),
             {"help": "Either diamond blastx results, or output protein FASTA from MetaEuk workflow"}),
            (("-o", "--output"),
             {"help": "Output location, default stdout", "default": "/dev/stdout"}),
        ),
        description="Parse FASTA and BLAST results to GFF3 format"
    )
    _parse_args(_ap)
    _main(_ap)
