#!/usr/bin/env python
"""
Installation:

pip install numpy==1.19.5 bcbio-gff==0.6.6 biopython==1.78 plumbum==1.6.9
"""
import os.path
from pathlib import Path
from typing import Generator, List

import numpy as np
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from plumbum import cli


class FeatureDensityHandler(cli.Application):
    fasta_file: Path
    gff3_file: Path
    bucket_size: int = 100000
    feature: str

    @staticmethod
    def get_file(file: str, err_message: str):
        file = Path(file).resolve()
        if not file.exists():
            print(err_message)
            exit(1)
        return file

    @cli.switch(["-f", "--fasta"], str, mandatory=True)
    def set_fasta_file(self, fasta_file: str):
        """
        Set FASTA file path
        """
        self.fasta_file = FeatureDensityHandler.get_file(fasta_file, "FASTA file not found")

    @cli.switch(["-b", "--bucket-size"], int)
    def set_bucket_size(self, bucket_size: int):
        """
        Set bucket size, default 100000
        """
        if not isinstance(bucket_size, int) or bucket_size < 1:
            print("Bucket size must be > 0")
            exit(1)
        self.bucket_size = bucket_size

    @cli.switch(["-g", "--gff3"], str, mandatory=True)
    def set_gff3_file(self, gff3_file):
        """
        Set GFF3 file path
        """
        self.gff3_file = FeatureDensityHandler.get_file(gff3_file, "GFF3 file not found")

    @cli.switch(["-q", "--query-feature"], str, mandatory=True)
    def set_query_feature(self, feature):
        """
        Set query feature for which to populate data (3rd column, e.g. gene, CDS, dispersed_repeat, etc.
        """
        self.feature = feature

    def get_features(self) -> Generator[SeqRecord, None, None]:
        """
        Return records parsed together with features
        """
        return GFF.parse(open(self.gff3_file),
                         base_dict=SeqIO.to_dict(SeqIO.parse(open(self.fasta_file), "fasta")),
                         limit_info={"gff_type": [self.feature]})

    def parse_to_buckets(self, out_path: Path):
        out_ptr = open(out_path, "w")
        sorted_records: List[SeqRecord] = sorted(self.get_features(), key=lambda rec: len(rec.seq), reverse=True)
        genome_length: int = sum([len(record) for record in sorted_records])
        pos_in_genome: int = 0
        buckets: np.ndarray = np.zeros(genome_length // self.bucket_size + 1, dtype="int64")
        for record in sorted_records:
            for feature in record.features:
                for i in feature.location:
                    buckets[(pos_in_genome + i) // self.bucket_size] += 1
            pos_in_genome += len(record.seq)
        out_ptr.write("\t".join(map(str, buckets)))
        out_ptr.write("\n")
        out_ptr.close()

    def main(self, *args):
        self.parse_to_buckets(Path(f"{os.path.splitext(self.gff3_file)[0]}.{self.feature}.density").resolve())


if __name__ == "__main__":
    FeatureDensityHandler.run()
