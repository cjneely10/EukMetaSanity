import glob
import os
import shutil
from collections import Counter
from pathlib import Path
from typing import List, Union, Type, Iterable

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from yapim import Task, DependencyInput, touch

from .merge_parallelized_output import merge
from .taxon_ids import augustus_taxon_ids


class Augustus(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "ab-gff3": self.wdir.joinpath(self.record_id + ".gff3"),
            "prot": self.wdir.joinpath(self.record_id + ".faa")
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [
            DependencyInput("MMSeqsConvertAlis")
        ]

    def run(self):
        """
        Run augustus
        """
        if self.parse_search_output(str(self.input["MMSeqsConvertAlis"]["results_files"][0])) == "":
            touch(str(self.output["ab-gff3"]))
            touch(str(self.output["prot"]))
            return
        # Initial training based on best species from taxonomy search
        out_gff = self._augustus(
            self.parse_search_output(str(self.input["MMSeqsConvertAlis"]["results_files"][0])), 1,
            str(self.input["fasta"])
        )
        if len(open(out_gff, "r").readlines()) < 200 or int(self.config["rounds"]) == 0:
            # Move any augustus-generated config stuff
            self._finalize_output(out_gff)
            return
        self._train_augustus(1, str(self.input["fasta"]), out_gff)
        # Remaining rounds of re-training on generated predictions
        for i in range(int(self.config["rounds"])):
            _last = False
            if i == int(self.config["rounds"]) - 1:
                _last = True
            out_gff = self._augustus(self.record_id + str(i + 1), i + 2, str(self.input["fasta"]), _last)
            if len(open(out_gff, "r").readlines()) < 200:
                break
            if i != int(self.config["rounds"]) - 1:
                self._train_augustus(i + 2, str(self.input["fasta"]), out_gff)
        # Move any augustus-generated config stuff
        self._finalize_output(out_gff)

    def _finalize_output(self, out_gff: str):
        self._handle_config_output()
        # Rename final file
        self.single(
            self.local["gffread"][out_gff, "-o", str(self.output["ab-gff3"]), "-G"],
            "5:00"
        )
        touch(str(self.output["prot"]))
        self.single(
            self.local["gffread"][
                str(self.output["ab-gff3"]), "-y", self.output["prot"], "-g", self.input["fasta"], "-S"
            ],
            "5:00"
        )

    @staticmethod
    def split_data(data: List[SeqRecord], n: int) -> List[List[SeqRecord]]:
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
        # Remore empty buckets
        buckets = [bucket for bucket in buckets if len(bucket) > 0]
        # Sanity check
        assert total_size == sum([sum([len(record.seq) for record in bucket]) for bucket in buckets])
        return buckets

    def _contig_splitter(self, created_files_list: List[str]) -> Iterable[str]:
        records = list(SeqIO.parse(self.input["fasta"], "fasta"))
        split_records: List[List[SeqRecord]] = Augustus.split_data(records, int(self.threads))
        for pos, records_list in enumerate(split_records):
            out_file = str(self.wdir.joinpath(f"{self.record_id}.{pos}.fasta"))
            with open(out_file, "w") as out_ptr:
                SeqIO.write(records_list, out_ptr, "fasta")
            created_files_list.append(out_file)
            yield out_file

    def _augustus(self, species: str, _round: int, _file: str, _last: bool = False) -> str:
        """ Run augustus training round

        :param species: Species string
        :param _round: Training round number
        :param _file: FASTA file
        :param _last: Is last training round
        :return: Path to output gff3 file
        """
        contig_files = []
        contig_files_iter = self._contig_splitter(contig_files)
        self.parallel(
            self.create_script(
                [self.program[
                     "--codingseq=on",
                     "--stopCodonExcludedFromCDS=false",
                     "--species=%s" % species,
                     "--outfile=%s" % contig_file + f".{_round}.gff",
                     ("--gff3=on" if _last else "--gff3=off"),
                     contig_file,
                 ] for contig_file in contig_files_iter], "augustus-runner.sh",
                parallelize=True
            )
        )

        out_gff = Path(os.path.join(
            self.wdir, Augustus.out_path(str(self.input["fasta"]), ".%i.gff" % _round)
        ))
        if out_gff.exists():
            self.local["rm"][out_gff]()
        merge(contig_files, out_gff, _round)
        return str(out_gff)

    def _handle_config_output(self):
        """
        Move the augustus training folders to their wdir folders
        """
        config_dir = os.path.join(
            os.path.dirname(os.path.dirname(Path(str(self.program)).resolve())),
            "config", "species"
        )
        for file in glob.glob(os.path.join(config_dir, self.record_id + "*")):
            shutil.move(
                file,
                os.path.join(self.wdir, os.path.basename(file))
            )

    def _train_augustus(self, _round: int, _file: str, out_gff: str):
        """ Run training set on augustus results

        :param _round: Round number of training
        :param _file: File being trained
        :param out_gff: Output gff3 path
        :return: Output gff3 path
        """
        # Remove old training directory, if needed
        config_dir = os.path.join(
            os.path.dirname(os.path.dirname(Path(str(self.program)).resolve())),
            "config", "species", self.record_id + str(_round)
        )
        if os.path.exists(config_dir):
            shutil.rmtree(config_dir)
        # Parse to genbank
        out_gb = os.path.join(self.wdir, Augustus.out_path(_file, ".%i.gb" % _round))
        self.single(
            self.local["gff2gbSmallDNA.pl"][
                out_gff,
                _file,
                "1000",
                out_gb
            ],
            "2:00"
        )
        species_config_prefix = self.record_id + str(_round)
        # Write new species config file
        self.single(
            self.local["new_species.pl"][
                "--species=%s" % species_config_prefix,
                out_gb
            ],
            "2:00"
        )
        # Run training
        self.single(
            self.local["etraining"][
                "--species=%s" % species_config_prefix,
                out_gb
            ],
            "10:00"
        )
        return out_gff

    @staticmethod
    def out_path(_file_name: str, _ext: str) -> str:
        """ Path join quick helper method to add new extension to file basename

        :param _file_name: Name of file
        :param _ext: New extension to give to file
        :return: New path
        """
        return os.path.basename(os.path.splitext(_file_name)[0]) + _ext

    def parse_search_output(self, search_results_file: str) -> str:
        """ Return optimal taxonomy from mmseqs search

        :param search_results_file: MMseqs results file
        :return: Augustus species with the most hits
        """
        augustus_ids_dict = augustus_taxon_ids()
        found_taxa = Counter()
        with open(search_results_file, "r") as _file:
            for line in _file:
                line = line.rstrip("\r\n").split()
                # Count those that pass the user-defined cutoff value
                if line[3] in augustus_ids_dict.keys() and float(line[2]) >= (float(self.config["cutoff"]) / 100.):
                    found_taxa[line[3]] += 1
        if len(found_taxa.most_common()) == 0:
            return ""
        return augustus_ids_dict[found_taxa.most_common()[0][0]]
