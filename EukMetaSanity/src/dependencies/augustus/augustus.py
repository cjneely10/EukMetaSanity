import glob
import os
import shutil
from collections import Counter
from copy import deepcopy
from pathlib import Path
from threading import Lock
from typing import List, Union, Type, Optional

from yapim import Task, DependencyInput, touch, clean

from .contig_splitter import ContigSplitter
from .taxon_ids import augustus_taxon_ids


class _UniqueIdentifiersFactory:
    _lock: Lock = Lock()
    _next: int = 256  # Begin after integers that are statically stored in memory

    @staticmethod
    def get_next() -> int:
        with _UniqueIdentifiersFactory._lock:
            _UniqueIdentifiersFactory._next += 1
            out = deepcopy(_UniqueIdentifiersFactory._next)
        return out


class Augustus(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "ab-gff3": self.wdir.joinpath(self.record_id + ".gff3"),
            "prot": self.wdir.joinpath(self.record_id + ".faa")
        }
        self._unique_id: int = _UniqueIdentifiersFactory.get_next()
        self._contig_splitter: Optional[ContigSplitter] = None

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return []

    @clean("*.gb", "*.gff")
    def run(self):
        self._contig_splitter = ContigSplitter(self.input["fasta"], int(self.threads), self.wdir, self.record_id)
        rounds = int(self.config["rounds"])
        if os.path.exists(self.input["search_results"]) and os.stat(self.input["search_results"]).st_size > 0:
            tax_search_results = self.parse_search_output(self.input["search_results"])
            if tax_search_results == "":
                touch(str(self.output["ab-gff3"]))
                touch(str(self.output["prot"]))
                return
            # Initial training based on best species from taxonomy search
            out_gff = self._augustus(tax_search_results, 1, str(self.input["fasta"]))
            rounds -= 1
            if rounds < 0:
                rounds = 0
        else:
            # Initial training based on provided gff3 file
            out_gff = self.input["gff3"]
        if rounds == 0 or self._line_count(out_gff) < 200:
            # Move any augustus-generated config stuff
            self._finalize_output(out_gff)
            return
        self._train_augustus(1, str(self.input["fasta"]), out_gff)
        # Remaining rounds of re-training on generated predictions
        for i in range(rounds):
            _last = False
            if i == rounds - 1:
                _last = True
            out_gff = self._augustus(self._species_config_prefix(i + 1), i + 2, str(self.input["fasta"]), _last)
            if Augustus._line_count(out_gff) < 200:
                break
            if i != rounds - 1:
                self._train_augustus(i + 2, str(self.input["fasta"]), out_gff)
        # Move any augustus-generated config stuff
        self._finalize_output(out_gff)
        # Remove temporarily-created files
        self._contig_splitter.finalize()

    @staticmethod
    def _line_count(file: str) -> int:
        with open(file, "r") as file_ptr:
            return sum(1 for _ in file_ptr)

    def _finalize_output(self, out_gff: str):
        self._handle_config_output()
        # Rename final file
        self.single(self.local["gffread"][out_gff, "-o", str(self.output["ab-gff3"]), "-G"], "5:00")
        touch(str(self.output["prot"]))
        self.single(
            self.local["gffread"][
                str(self.output["ab-gff3"]), "-y", self.output["prot"], "-g", self.input["fasta"], "-S"
            ],
            "5:00"
        )

    def _augustus(self, species: str, _round: int, _file: str, _last: bool = False) -> str:
        """ Run augustus training round

        :param species: Species string
        :param _round: Training round number
        :param _file: FASTA file
        :param _last: Is last training round
        :return: Path to output gff3 file
        """
        assert self._contig_splitter is not None
        self.parallel(
            self.create_script(
                [self.program[
                     "--codingseq=on",
                     "--stopCodonExcludedFromCDS=false",
                     "--species=%s" % species,
                     "--outfile=%s" % contig_file + f".{_round}.gff",
                     ("--gff3=on" if _last else "--gff3=off"),
                     contig_file,
                 ] for contig_file in self._contig_splitter], "augustus-runner.sh",
                parallelize=True
            )
        )

        out_gff = Path(os.path.join(self.wdir, Augustus.out_path(str(self.input["fasta"]), ".%i.gff" % _round)))
        if out_gff.exists():
            self.local["rm"][out_gff]()
        self._contig_splitter.merge(out_gff, _round)
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
            shutil.move(file, os.path.join(self.wdir, os.path.basename(file)))

    def _species_config_prefix(self, _round: int) -> str:
        """
        Generate species prefix from training round

        :param _round: training round
        :return: species identifier corresponding to provided `_round`
        """
        return f"EukMS-{self.record_id}-{self._unique_id}-{_round}"

    def _train_augustus(self, _round: int, _file: str, out_gff: str):
        """ Run training set on augustus results

        :param _round: Round number of training
        :param _file: File being trained
        :param out_gff: Output gff3 path
        :return: Output gff3 path
        """
        species_config_prefix = self._species_config_prefix(_round)
        # Remove old training directory, if needed
        config_dir = os.path.join(
            os.path.dirname(os.path.dirname(Path(str(self.program)).resolve())),
            "config", "species", species_config_prefix
        )
        if os.path.exists(config_dir):
            shutil.rmtree(config_dir)
        # Parse to genbank
        out_gb = os.path.join(self.wdir, Augustus.out_path(_file, ".%i.gb" % _round))
        self.single(self.local["gff2gbSmallDNA.pl"][out_gff, _file, "1000", out_gb], "2:00")
        # Write new species config file
        self.single(self.local["new_species.pl"]["--species=%s" % species_config_prefix, out_gb], "2:00")
        # Run training
        self.single(self.local["etraining"]["--species=%s" % species_config_prefix, out_gb], "10:00")
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
