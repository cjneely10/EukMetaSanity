#!/usr/bin/env python
"""
Installation:

pip install plumbum==1.6.9
"""
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, wait
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import Optional

from plumbum import cli


@dataclass
class Region:
    class LocusIdentifier(Enum):
        FP = auto()  # False positive
        FN = auto()  # False negative
        TP = auto()  # True positive
        NONE = auto()  # Default: None

    start: int
    end: int
    overlap: int = 0
    locus_identifier: LocusIdentifier = LocusIdentifier.NONE

    def __len__(self):
        return self.end - self.start + 1

    def __eq__(self, other: "Region") -> bool:
        return max(self.start, other.start) <= min(self.end, other.end)

    def overlaps(self, other: "Region") -> bool:
        return self.__eq__(other)

    def __ne__(self, other: "Region") -> bool:
        return not self == other

    def __gt__(self, other: "Region") -> bool:
        return self.start > other.end

    def __lt__(self, other: "Region") -> bool:
        return self.end < other.start

    def __str__(self) -> str:
        return f"{self.start} {self.end} {self.locus_identifier.value} {self.overlap}"

    def overlap_length(self, other: "Region"):
        return min(self.end, other.end) - max(self.start, other.start) + 1


class File:
    def __init__(self, file_path: Path, feature_type: str):
        self._data = defaultdict(list)
        self.load(file_path, feature_type)

    @property
    def data(self):
        return self._data

    def load(self, file: Path, feature_type: str):
        file_ptr = open(file, "r")
        for line in file_ptr:
            if "#" in line:
                continue
            line = line.split(maxsplit=5)
            if len(line) < 5 or feature_type not in line[2]:
                continue
            self._data[line[0]].append(Region(int(line[3]), int(line[4])))
        # Sort results in case out of order
        for key in self._data.keys():
            self._data[key].sort(key=lambda region: region.start)


class RegionOverlapper:
    def __init__(self, query_file: File, reference_file: File, threshold: float):
        self._query: File = query_file
        self._reference: File = reference_file
        self._keys = {*self._query.data.keys(), *self._reference.data.keys()}
        self._run(threshold)
        self.threshold = threshold
        for key in self._keys:
            for val in self.query[key]:
                if val.locus_identifier == Region.LocusIdentifier.NONE:
                    raise AttributeError(f"Query locus not set at {val}")
            for val in self.reference[key]:
                if val.locus_identifier == Region.LocusIdentifier.NONE:
                    raise AttributeError(f"Reference locus not set at {val}")

    @property
    def keys(self) -> set:
        return self._keys

    @property
    def query(self) -> defaultdict:
        return self._query.data

    @property
    def reference(self) -> defaultdict:
        return self._reference.data

    def _run(self, threshold: float):
        for contig_id in self._keys:
            self._assess_region(contig_id, threshold)

    def _assess_region(self, contig_id: str, threshold: float):
        query_pos = 0
        ref_pos = 0
        while query_pos < len(self.query[contig_id]) and ref_pos < len(self.reference[contig_id]):
            query_region: Region = self.query[contig_id][query_pos]
            ref_region: Region = self.reference[contig_id][ref_pos]
            if query_region.overlaps(ref_region):
                overlap_fraction = ref_region.overlap_length(query_region) / len(ref_region)
                if overlap_fraction >= threshold and ref_region.locus_identifier != Region.LocusIdentifier.TP:
                    query_region.locus_identifier = Region.LocusIdentifier.TP
                    ref_region.locus_identifier = Region.LocusIdentifier.TP
                    query_pos += 1
                    ref_pos += 1
                else:
                    query_region.locus_identifier = Region.LocusIdentifier.FP
                    query_pos += 1
            else:
                if query_region > ref_region:
                    ref_region.locus_identifier = Region.LocusIdentifier.FN
                    ref_pos += 1
                else:
                    query_region.locus_identifier = Region.locus_identifier.FP
                    query_pos += 1
        self._consume_ending(query_pos, ref_pos, contig_id)

    def _consume_ending(self, query_pos: int, ref_pos: int, contig_id: str):
        while query_pos != len(self.query[contig_id]):
            self._query.data[contig_id][query_pos].locus_identifier = Region.LocusIdentifier.FP
            query_pos += 1
        while ref_pos != len(self.reference[contig_id]):
            self._reference.data[contig_id][ref_pos].locus_identifier = Region.LocusIdentifier.FN
            ref_pos += 1

    @property
    def fp(self) -> int:
        out = 0
        for record_id in self.keys:
            for region in self.query[record_id]:
                if region.locus_identifier == Region.LocusIdentifier.FP:
                    out += 1
        return out

    @property
    def fn(self) -> int:
        out = 0
        for record_id in self.keys:
            for region in self.reference[record_id]:
                if region.locus_identifier == Region.LocusIdentifier.FN:
                    out += 1
        return out

    @property
    def tp(self) -> int:
        out = 0
        for record_id in self.keys:
            for region in self.reference[record_id]:
                if region.locus_identifier == Region.LocusIdentifier.TP:
                    out += 1
        return out

    @property
    def sensitivity(self) -> Optional[float]:
        tp = self.tp
        denominator = tp + self.fn
        if denominator != 0:
            return tp / denominator

    @property
    def precision(self) -> Optional[float]:
        tp = self.tp
        denominator = tp + self.fp
        if denominator != 0:
            return tp / denominator

    @property
    def f1(self) -> Optional[float]:
        precision = self.precision
        recall = self.sensitivity
        if precision is None or recall is None:
            return None
        if precision + recall == 0.0:
            return 0.0
        return 2 * (precision * recall) / (precision + recall)


class MatchRunner(cli.Application):
    """
    Compare GFF3 feature position for sensitivity/precision metrics
    """
    reference_file: Path
    query_file: Path
    reference_feature: str
    query_feature: str

    @cli.switch(["-qf", "--query-feature"], str)
    def set_query_feature(self, feature_type: str):
        """
        Set query feature type to track. If not set, will be the same as --rf
        """
        self.query_feature = feature_type

    @cli.switch(["-q", "--query"], str, mandatory=True)
    def set_query_file(self, query_path: str):
        """
        Path to query file
        """
        query_path = Path(query_path).resolve()
        if not query_path.exists():
            print("Cannot locate query file")
            exit(1)
        self.query_file = query_path

    @cli.switch(["-rf", "--reference-feature"], str, mandatory=True)
    def set_reference_feature(self, feature_type: str):
        """
        Set reference feature type to track
        """
        self.reference_feature = feature_type

    @cli.switch(["-r", "--reference"], str, mandatory=True)
    def set_reference_file(self, reference_path: str):
        """
        Path to reference file
        """
        reference_path = Path(reference_path).resolve()
        if not reference_path.exists():
            print("Cannot locate reference file")
            exit(1)
        self.reference_file = reference_path

    def main(self, record_id):
        if self.query_feature is None:
            self.query_feature = self.reference_feature
        out = []
        _range = [i * 0.1 for i in range(1, 11)]
        with ThreadPoolExecutor() as executor:
            futures = []
            for threshold in _range:
                query_file = File(self.query_file, self.query_feature)
                ref_file = File(self.reference_file, self.reference_feature)
                futures.append(executor.submit(RegionOverlapper, query_file, ref_file, threshold))
            wait(futures)
            for future in futures:
                overlaps = future.result()
                out.extend([overlaps.threshold, overlaps.sensitivity, overlaps.precision, overlaps.f1])
        print(" ".join(map(str, [record_id, *out])))


if __name__ == "__main__":
    MatchRunner.run()
