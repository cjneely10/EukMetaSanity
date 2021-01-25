#!/usr/bin/env python
import os
from pathlib import Path
from BCBio import GFF
from plumbum import cli
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


class TierCounter(cli.Application):
    @staticmethod
    def get_tier_values(file_path: Path) -> (dict, dict, int):
        file = open(file_path, "r")
        rec: SeqRecord
        tier_values = {1: 0, 2: 0}
        assignment_values = {"genemark": 0, "metaeuk": 0}
        total = 0
        for rec in GFF.parse(file):
            feature: SeqFeature
            for feature in rec.features:
                if feature.type == "locus":
                    total += 1
                    tier = set()
                    assignments = set()
                    for qualifier in feature.qualifiers["transcripts"]:
                        # metaeuk
                        if qualifier.startswith("gene"):
                            tier.add(1)
                            assignments.add(1)
                        # GeneMark
                        if qualifier.endswith("_t"):
                            tier.add(2)
                            assignments.add(2)
                        if len(tier) == 2:
                            break
                    for i in range(len(tier)):
                        tier_values[i + 1] += 1
                    if 1 in assignments:
                        assignment_values["metaeuk"] += 1
                    if 2 in assignments:
                        assignment_values["genemark"] += 1
        file.close()
        return tier_values, assignment_values, total

    def main(self):
        _in = Path(input()).resolve()
        print("ID", "Tier1", "Tier2", "GeneMark", "MetaEuk", "total_prots")
        while _in is not None:
            tier_values, assignment_values, total = TierCounter.get_tier_values(_in)
            print(os.path.basename(_in), tier_values[1], tier_values[2], assignment_values["genemark"],
                  assignment_values["metaeuk"], total)
            _in = Path(input()).resolve()


if __name__ == "__main__":
    TierCounter.run()
