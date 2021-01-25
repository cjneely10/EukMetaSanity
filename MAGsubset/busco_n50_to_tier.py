#!/usr/bin/env python
import os
import glob
from pathlib import Path

import numpy as np
from Bio import SeqIO
from plumbum import cli
import matplotlib.pyplot as plt


class BuscoN50TierGrapher(cli.Application):
    _directory: Path
    _tier_stats_file: Path

    @cli.switch(["-d", "--directory"], str, mandatory=True)
    def set_directory(self, directory):
        if not os.path.exists(directory):
            raise FileNotFoundError(directory + " does not exist")
        self._directory = Path(directory).resolve()

    @cli.switch(["-p", "--tier-file-path"], str, mandatory=True)
    def set_tier_summary_path(self, tier_file_path):
        if not os.path.exists(tier_file_path):
            raise FileNotFoundError(tier_file_path)
        self._tier_stats_file = Path(tier_file_path).resolve()

    @staticmethod
    def get_completion_score(file_path: str) -> float:
        with open(file_path, "r") as r:
            for line in r:
                line = line.split()
                if len(line) > 0:
                    if line[0][0] == "C":
                        line = next(r)
                        return float(line.split()[0])

    @staticmethod
    def calculate_n50(fasta_file: str, min_contig_size: int = 500) -> tuple:
        """ Determine n50 of a genome in a FASTA file given a minimum contig size

        :param fasta_file: Path to FASTA file
        :param min_contig_size: Minimum contig size to consider in n50 count
        :raises: FileNotFoundError/TypeError for improper type forming, Attribute error for empty file
        :return: n50 and length of genome in file
        """
        if not os.path.exists(fasta_file):
            raise FileNotFoundError
        if min_contig_size <= 0:
            raise TypeError("min_contig_size must be positive")

        # Get sorted lengths
        record_lengths = [len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
        # Confirm non-empty file
        if len(record_lengths) == 0:
            raise AttributeError("File contains no records")

        # Sort descending
        record_lengths.sort(reverse=True)

        # Calculate n50
        total_length = sum(record_lengths)
        half_len = total_length // 2
        total = 0
        pos = 0
        # Stop once passed halfway point
        while total <= half_len and record_lengths[pos] >= min_contig_size:
            total += record_lengths[pos]
            pos += 1
        return record_lengths[pos - 1], total_length, record_lengths[0]

    def main(self):
        n50s = []
        buscos = []
        tier_ratios = []
        tier_fp = open(self._tier_stats_file, "r")
        next(tier_fp)
        for line in tier_fp:
            line = line.split(" ")
            record_id = line[0].replace(".all.gff3", "")
            print(record_id)
            tier_ratios.append(float(line[2]) / float(line[1]))
            busco_summary_file = glob.glob(os.path.join(
                self._directory,
                record_id,
                "busco-" + record_id + ".all.1",
                "short_summary*"
            ))
            if len(busco_summary_file) == 0:
                busco_summary_file = glob.glob(os.path.join(
                    self._directory,
                    record_id,
                    "busco-" + record_id,
                    "short_summary*"
                ))
            busco_summary_file = busco_summary_file[0]
            mag_fasta_file = os.path.join(
                self._directory,
                record_id,
                record_id + ".fna"
            )
            buscos.append(BuscoN50TierGrapher.get_completion_score(busco_summary_file))
            n50s.append(BuscoN50TierGrapher.calculate_n50(mag_fasta_file)[0])

        buscos = np.array(buscos)
        tier_ratios = np.array(tier_ratios)
        n50s = np.array(n50s)

        plt.plot(buscos, tier_ratios, label="BUSCO scores", color="bx")
        plt.plot(n50s, tier_ratios, label="Assembly N50", color="ro")

        plt.savefig("TiervBuscoN50.png")


if __name__ == "__main__":
    BuscoN50TierGrapher.run()
