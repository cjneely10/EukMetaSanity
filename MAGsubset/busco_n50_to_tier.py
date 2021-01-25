#!/usr/bin/env python
import os
from pathlib import Path
from plumbum import cli


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

    def main(self):
        # Gather record ids

        # Gather busco summary paths

        # Gather
        tier_fp = open(self._tier_stats_file, "r")
        next(tier_fp)
        for line in tier_fp:
            pass


if __name__ == "__main__":
    BuscoN50TierGrapher.run()
