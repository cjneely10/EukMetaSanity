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
        sorted_records: List[SeqRecord] = sorted(self.get_features(), key=lambda rec: len(rec.seq), reverse=True)
        genome_length: int = sum([len(record) for record in sorted_records])
        pos_in_genome: int = 0
        buckets: np.ndarray = np.zeros(genome_length // self.bucket_size + 1, dtype="int64")
        for record in sorted_records:
            for feature in record.features:
                feature_start = feature.location.start + pos_in_genome
                feature_end = feature.location.end + pos_in_genome
                start = self.start_bucket_idx(feature_start)
                end = self.end_bucket_idx(feature_end)
                if start - 1 != end:
                    buckets[start - 1] += self.amt_at_start(start, feature_start)
                    buckets[end] += self.amt_at_end(end, feature_end)
                    for i in self.overlapping_bucket_range(start, end):
                        buckets[i] += self.bucket_size
                else:
                    buckets[start - 1] += feature_end - feature_start
            pos_in_genome += len(record.seq)
        self.write_results(out_path, buckets)

    def start_bucket_idx(self, start: int):
        return start // self.bucket_size + 1

    def end_bucket_idx(self, end: int):
        return end // self.bucket_size

    def amt_at_start(self, start_idx: int, start: int):
        return self.bucket_size * start_idx - start

    def amt_at_end(self, end_idx: int, end: int):
        return end - self.bucket_size * end_idx

    def overlapping_bucket_range(self, start: int, end: int):
        return range(self.start_bucket_idx(start), self.end_bucket_idx(end))

    def write_results(self, out_path: Path, buckets: np.ndarray):
        out_ptr = open(out_path, "w")
        out_ptr.write(str(round(buckets[0] / self.bucket_size, 2)))
        for bucket in buckets[1:]:
            out_ptr.write("\t")
            out_ptr.write(str(round(bucket / self.bucket_size, 2)))
        out_ptr.write("\n")
        out_ptr.close()

    def main(self):
        self.parse_to_buckets(Path(f"{os.path.splitext(self.gff3_file)[0]}.{self.feature}.density").resolve())


if __name__ == "__main__":
    FeatureDensityHandler.run()
