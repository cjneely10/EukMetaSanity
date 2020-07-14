#!/usr/bin/env python3
import os
from BCBio import GFF
from Bio import SeqIO
from decimal import Decimal
from _io import TextIOWrapper
from operator import itemgetter
from collections import namedtuple
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from EukMetaSanity.utils.arg_parse import ArgParse
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Result as read in from DNA/contig side as subject
Result = namedtuple("Result", ("loc_type", "sstart", "send", "strand", "score"))


# # Available parsing types
def diamond(diamond_fp, data, ap):
    # Ensure is passed file pointer
    assert isinstance(data, defaultdict)
    assert isinstance(diamond_fp, TextIOWrapper)
    # Parse file
    for line in diamond_fp:
        line = line.rstrip("\r\n").split("\t")
        assert len(line) == 12, "Diamond file is improperly formatted"
        # Convert hit metrics
        evalue = Decimal(line[10])
        pident = float(line[2])
        aligned_length = int(line[3])
        sstart, send = int(line[8]), int(line[9])
        qstart, qend = int(line[6]), int(line[7])
        # Coverage based on mode
        if ap.args.mode == "subject":
            coverage = aligned_length / (send - sstart + 1)
        else:
            coverage = aligned_length / (qend - qstart + 1)
        # Only valid results
        if evalue <= ap.args.evalue and pident >= ap.args.pident and coverage >= ap.args.coverage:
            # Per --outfmt 6 default
            data[line[0]].append(
                Result(
                    loc_type="CDS",
                    sstart=sstart,
                    send=send,
                    strand=0,
                    score=".",
                )
            )


def _parse_diamond(i, ap, feature_data, record, rec):
    # Determine source and add score
    qualifiers = {"source": ap.args.source, "score": ".", "ID": "gene%i" % i}
    # Store gene/exon data
    spans = _reduce_span([(feature.sstart, feature.send) for feature in feature_data.get(record.id, [])])
    for span, feature in zip(spans, feature_data.get(record.id, [])):
        rec.features.append(
            SeqFeature(
                FeatureLocation(span[0], span[1]), type="gene", strand=feature.strand, qualifiers=qualifiers
            )
        )
        # Default retain codon and CDS information
        rec.features[-1].sub_features = [
            SeqFeature(
                FeatureLocation(span[0], span[0] + 3), type="start_codon", strand=feature.strand,
                qualifiers={"source": ap.args.source}
            ),
            SeqFeature(
                FeatureLocation(span[0], span[1]), type="CDS", strand=feature.strand,
                qualifiers={"source": ap.args.source}
            ),
            SeqFeature(
                FeatureLocation(span[1] - 3, span[1]), type="stop_codon", strand=feature.strand,
                qualifiers={"source": ap.args.source}
            ),
        ]
    i += 1
    return i


def metaeuk(metaeuk_file_path, data, ap):
    records_fp = SeqIO.parse(metaeuk_file_path, "fasta")
    for record in records_fp:
        # Metaeuk header contains pipes as delimiters - be sure to remove all in input seqs!
        line = record.id.split("|")
        # Retain gene info in nested lists
        recs = []
        if Decimal(line[4]) < ap.args.evalue:
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
                for _type in ("exon", "CDS"):
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
        # First value is gene
        rec.features.append(
            SeqFeature(
                FeatureLocation(features[0].sstart, features[0].send), type=features[0].loc_type,
                strand=features[0].strand,
                qualifiers={"score": features[0].score, "source": ap.args.source, "ID": "gene%i" % i}
            )
        )
        # Remaining in list are features that describe the gene
        rec.features[-1].sub_features = []
        for feature in features[1:]:
            rec.features[-1].sub_features.append(
                SeqFeature(
                    FeatureLocation(feature.sstart, feature.send), type=feature.loc_type,
                    strand=features[0].strand,
                    qualifiers={"score": feature.score, "source": ap.args.source}
                )
            )
        i += 1
    return i


# # Helper functions
def _reduce_span(coords_list):
    # Sort coordinates by start value
    ranges_in_coords = sorted(coords_list, key=itemgetter(0))
    # Will group together matching sections into spans
    # Return list of these spans at end
    # Initialize current span and list to return
    spans_in_coords = [list(ranges_in_coords[0]), ]
    for coords in ranges_in_coords[1:]:
        # The start value is within the farthest range of current span
        # and the end value extends past the current span
        if coords[0] <= spans_in_coords[-1][1] < coords[1]:
            spans_in_coords[-1][1] = coords[1]
        # The start value is past the range of the current span
        # Append old span to list to return
        # Reset current span to this new range
        elif coords[0] > spans_in_coords[-1][1]:
            spans_in_coords.append(list(coords))
    return tuple(map(tuple, spans_in_coords))


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
            i = globals()["_parse_" + ap.args.source](i, ap, feature_data, record, rec)
        yield rec


# # Driver logic
def _parse_args(ap):
    # Confirm path existence
    assert os.path.exists(ap.args.data_file)
    assert os.path.exists(ap.args.fasta_file)
    # Confirm valid mode and parsing source
    assert ap.args.mode in ("query", "subject")
    assert ap.args.source in globals().keys()
    # Convert numerics
    try:
        ap.args.evalue = Decimal(ap.args.evalue)
        ap.args.pident = float(ap.args.pident)
        ap.args.coverage = float(ap.args.coverage)
        if ap.args.coverage > 1.0:
            raise ValueError("Coverage must be less than 1.0")
    except ValueError as e:
        print(e)
        exit(1)


def _main(ap):
    # Store results data
    data = defaultdict(list)
    # Call parsing function - ensured to exist
    globals()[ap.args.source](open(ap.args.data_file), data, ap)
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
            (("-e", "--evalue"),
             {"help": "Max evalue to keep, default 1E-5", "default": "1E-5"}),
            (("-p", "--pident"),
             {"help": "Minimum percent identity to retain, default 90.0", "default": "90.0"}),
            (("-c", "--coverage"),
             {"help": "Minimum coverage to retain, default 0.60", "default": "0.60"}),
            (("-m", "--mode"),
             {"help": "Map to subject/query, default subject", "default": "subject"}),
            (("-s", "--source"),
             {"help": "Source, select from diamond/metaeuk default diamond", "default": "metaeuk"}),
        ),
        description="Parse FASTA and BLAST results to GFF3 format"
    )
    _parse_args(_ap)
    _main(_ap)
