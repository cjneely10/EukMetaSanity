import os
from collections import namedtuple, defaultdict
from typing import List, Union, Type, TextIO, Dict, Generator

from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from yapim import Task, DependencyInput, prefix

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
                send=int(line[7]) + 1,
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
                    start, end = end, start
                recs.append(
                    Result(
                        loc_type=_type,
                        sstart=start,
                        send=end + 1,
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


class MetaEuk(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "gff3": os.path.join(self.wdir, self.record_id + ".gff3"),
            "prot": os.path.join(self.wdir, self.record_id + ".faa")
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return []

    def run(self):
        """
        Run metaeuk
        """
        if len(self.data) == 0:
            return
        database = self.data[0]
        is_profile = []
        if "p:" in database:
            is_profile.append("--slice-search")
            database = database[2:]
        db_prefix = prefix(database)
        _outfile = os.path.join(self.wdir, "%s_%s" % (self.record_id, db_prefix))
        # Run MetaEuk
        self.parallel(
            self.program[
                "easy-predict",
                str(self.input["fasta"]),
                database,
                _outfile,
                os.path.join(self.wdir, "tmp"),
                "--threads", self.threads,
                (*self.added_flags),
                (*is_profile),
            ]
        )
        # Store results data
        data = defaultdict(list)
        # Call parsing function - ensured to exist
        metaeuk(open(_outfile + ".fas"), data)
        # Write results
        GFF.write(_iter_gff(str(self.input["fasta"]), data), open(str(self.output["gff3"]), "w"))
        # Rename output file
        os.replace(_outfile + ".fas", str(self.output["prot"]))
