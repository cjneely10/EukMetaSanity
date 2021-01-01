#!/usr/bin/env python3

# pylint: disable=invalid-name
"""
Parse metaeuk results to gff3 format
"""

import os
from typing import Dict, List, TextIO, Generator
from collections import namedtuple
from collections import defaultdict
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from EukMetaSanity.utils.arg_parse import ArgParse

# Result as read in from DNA/contig side as subject
Result = namedtuple("Result", ("loc_type", "sstart", "send", "strand"))


def metaeuk(metaeuk_file_ptr: TextIO, data: Dict[str, List[List[Result]]]):
    """ Read metaeuk file into data dictionary

    :param metaeuk_file_ptr: File ptr to parse
    :param data: Reference to dictionary to load based on file
    """
    records_fp = SeqIO.parse(metaeuk_file_ptr, "fasta")
    for record in records_fp:
        # Metaeuk header contains pipes as delimiters - be sure to remove all in input seqs!
        line = record.id.split("|")
        # Retain gene info in nested lists
        recs = []
        # Determine strand
        if line[2] == "+":
            strand = 1
        else:
            strand = -1
        # Add base record as first record in nested list
        recs.append(
            Result(
                loc_type="gene",
                sstart=int(line[6]),
                send=int(line[7]),
                strand=strand,
            )
        )
        # Store exon/CDS info from rest of header at nested list loc
        all_coords = line[8:]
        if strand == -1:
            all_coords.reverse()
        for coords in all_coords:
            start, end, _ = coords.split(":")
            for _type in ("CDS",):
                start = int(start.split("[")[1][:-1])
                end = int(end.split("[")[1][:-1])
                if strand < 0:
                    start, end = end + 1, start + 1
                recs.append(
                    Result(
                        loc_type=_type,
                        sstart=start,
                        send=end,
                        strand=strand,
                    )
                )
        recs.sort(key=lambda res: res.sstart)
        # Add nested list
        data[line[1]].append(recs)


def _parse_metaeuk(i: int, feature_data: Dict[str, List[List[Result]]], record: SeqRecord, rec: SeqRecord):
    """ Parse metaeuk features from dict into gff3 format

    :param i: Gene idx
    :param feature_data: Reference to dict of all gene data
    :param record: Associated FASTA record
    :param rec: Associated GFF3 full record
    :return: Next gene index after this group has been processed
    """
    for features in feature_data.get(record.id, [[]]):
        if len(features) == 0:
            continue
        # First value is gene
        rec.features.append(
            SeqFeature(
                FeatureLocation(features[0].sstart, features[0].send), type=features[0].loc_type,
                strand=features[0].strand,
                qualifiers={"source": "metaeuk", "ID": "gene%i" % i}
            )
        )
        # Remaining in list are features that describe the gene
        rec.features[-1].sub_features = []
        for feature in features[1:]:
            rec.features[-1].sub_features.append(
                SeqFeature(
                    FeatureLocation(feature.sstart, feature.send), type=feature.loc_type,
                    strand=features[0].strand,
                    qualifiers={"source": "metaeuk"}
                )
            )
        i += 1
    return i


def _iter_gff(fasta_file: str, feature_data: Dict[str, List[List[Result]]]) -> Generator[SeqRecord, None, None]:
    """ Load FASTA records and build new records that include gff3 feature data

    :param fasta_file: Path to fasta data
    :param feature_data: Feature dict to build
    :return: Iterator of records
    """
    assert isinstance(feature_data, defaultdict)
    fasta_fp = SeqIO.parse(fasta_file, "fasta")
    i = 1
    for record in fasta_fp:
        # Each GFF entry based on FASTA record
        rec = SeqRecord(id=record.id, seq=record.seq)
        feature_list = feature_data.get(record.id, [])
        if len(feature_list) > 0:
            i = _parse_metaeuk(i, feature_data, record, rec)
        yield rec


def _parse_args(ap: ArgParse):
    """ Confirm path existence

    :param ap: Reference to ArgParse object
    :raises: AssertionError for invalid paths
    """
    assert os.path.exists(ap.args.data_file)
    assert os.path.exists(ap.args.fasta_file)


def _main(ap):
    """ Main logic

    :param ap: Reference to ArgParse object
    """
    # Store results data
    data = defaultdict(list)
    # Call parsing function - ensured to exist
    metaeuk(open(ap.args.data_file), data)
    # Write results
    GFF.write(_iter_gff(ap.args.fasta_file, data), open(ap.args.output, "w"))


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
