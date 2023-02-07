import os
from pathlib import Path
from typing import List, Iterator

from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


class _GeneCount:
    def __init__(self):
        self._count = 0

    @property
    def n(self):
        self._count += 1
        return self._count


class ContigSplitter:
    def __init__(self, fasta_file: Path, threads: int, wdir: Path, record_id: str):
        assert fasta_file.exists()
        self._files: List[str] = ContigSplitter._split_contigs(fasta_file, threads, wdir, record_id)

    def __iter__(self) -> Iterator[str]:
        return iter(self._files)

    def merge(self, output_path: Path, _round: int):
        output_data: List[SeqRecord] = []
        gene_count = _GeneCount()
        for file in self._files:
            gff_file = f"{file}.{_round}.gff"
            ContigSplitter._collect(Path(gff_file), output_data, gene_count)
            os.remove(gff_file)
        with open(output_path, "w") as out_ptr:
            GFF.write(output_data, out_ptr)
        output_data.clear()
        del output_data

    def finalize(self):
        for file in self._files:
            if os.path.exists(file):
                os.remove(file)

    @staticmethod
    def _split_contigs(fasta_file: Path, threads: int, wdir: Path, record_id: str) -> List[str]:
        records = list(SeqIO.parse(str(fasta_file), "fasta"))
        split_records: List[List[SeqRecord]] = ContigSplitter._split_data(records, threads)
        out = []
        for pos, records_list in enumerate(split_records):
            out_file = str(wdir.joinpath(f"{record_id}.{pos}.fasta"))
            with open(out_file, "w") as out_ptr:
                SeqIO.write(records_list, out_ptr, "fasta")
            out.append(out_file)
        for spl_rec in split_records:
            spl_rec.clear()
        split_records.clear()
        del split_records
        records.clear()
        del records
        return out

    @staticmethod
    def _split_data(data: List[SeqRecord], n: int) -> List[List[SeqRecord]]:
        """Accepts list of contig lengths, splits to n lists of indices"""
        buckets: List[List[SeqRecord]] = [[] for _ in range(n + 1)]
        sums = {i: 0 for i in range(len(buckets))}
        total_size = sum([len(rec.seq) for rec in data])
        bucket_size = total_size // n

        # Track contigs that are too big to fit into a bucket,
        too_large = []
        # and contigs that couldn't fit because buckets are all mostly full
        unable_to_fit = []
        for record in data:
            record_length = len(record.seq)
            # Contig too large for bucket size - place into own bucket
            if record_length > bucket_size:
                too_large.append([record])
                continue
            # Try to find a bucket that can fit contig
            pos = 0
            while pos < len(buckets):
                current_sum = sums[pos]
                if current_sum + record_length > bucket_size:
                    pos += 1
                else:
                    buckets[pos].append(record)
                    sums[pos] += len(record)
                    break
            # Reached end and did not find appropriate bucket, collect into own bucket
            if pos == len(buckets):
                unable_to_fit.append(record)

        # Add bucket of records that could not be placed
        buckets.append(unable_to_fit)
        # Extend using records that were too large to fit
        buckets.extend(too_large)
        # Remove empty buckets
        buckets = [bucket for bucket in buckets if len(bucket) > 0]
        # Sanity check
        assert total_size == sum([sum([len(record.seq) for record in bucket]) for bucket in buckets])
        return buckets

    @staticmethod
    def _collect(file: Path, output_data: List[SeqRecord], counter: _GeneCount):
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
