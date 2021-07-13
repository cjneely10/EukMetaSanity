import os
from pathlib import Path
from typing import List

from BCBio import GFF
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


class GeneCount:
    def __init__(self):
        self._count = 0

    @property
    def n(self):
        self._count += 1
        return self._count


def merge(files: List[Path], output_path: Path, _round: int):
    output_data: List[SeqRecord] = []
    gene_count = GeneCount()
    for file in files:
        _collect(Path(str(file) + f".{_round}.gff"), output_data, gene_count)
        os.remove(file)
        os.remove(str(file) + f".{_round}.gff")
    with open(output_path, "w") as out_ptr:
        GFF.write(output_data, out_ptr)


def _collect(file: Path, output_data: List[SeqRecord], counter: GeneCount):
    with open(file, "r") as file_ptr:
        record: SeqRecord
        for record in GFF.parse(file_ptr):
            new_record = SeqRecord(
                id=record.id,
                seq=record.seq
            )
            feature: SeqFeature
            for feature in record.features:
                qualifiers = {
                    "source": "augustus",
                    "ID": f"gene{counter.n}-{str(feature.id)}"
                }
                top_feature = SeqFeature(
                    feature.location,
                    type=feature.type,
                    strand=feature.strand,
                    qualifiers=qualifiers
                )
                top_feature.sub_features = []
                for sub_feature in feature.sub_features:
                    sub_feature.qualifiers = {
                        "source": "augustus",
                    }
                    top_feature.sub_features.append(sub_feature)
                new_record.features.append(top_feature)
            output_data.append(new_record)