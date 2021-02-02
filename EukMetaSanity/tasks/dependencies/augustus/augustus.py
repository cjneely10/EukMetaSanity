"""
Module holds augustus build functionality
"""
import glob
import os
import shutil
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from EukMetaSanity import Task, TaskList, program_catch, DependencyInput, touch, set_complete
from EukMetaSanity.tasks.dependencies.augustus.taxon_ids import augustus_taxon_ids


class AugustusIter(TaskList):
    """ TaskList iterates over augustus tasks

    name: augustus

    requires:

    depends: mmseqs.convertalis

    expects: fasta[Path]

    output: ab-gff3[Path]

    config:
        augustus:
          program: augustus
          cutoff: 25.0
          rounds: 3

    """
    name = "augustus"
    requires = []
    depends = [DependencyInput("mmseqs.convertalis")]

    class Augustus(Task):
        """
        Task class handles augustus task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "ab-gff3": os.path.join(self.wdir, self.record_id + ".augustus.gff3")
            }

        @program_catch
        def run(self):
            """
            Run augustus
            """
            # Initial training based on best species from taxonomy search
            out_gff = self._augustus(
                self.parse_search_output(str(self.input["mmseqs.convertalis"]["results_files"][0])), 1,
                str(self.dependency_input["fasta"])
            )
            self._train_augustus(1, str(self.dependency_input["fasta"]), out_gff)
            # Remaining rounds of re-training on generated predictions
            for i in range(int(self.config["rounds"])):
                _last = False
                if i == int(self.config["rounds"]) - 1:
                    _last = True
                out_gff = self._augustus(self.record_id + str(i + 1), i + 2, str(self.dependency_input["fasta"]), _last)
                if i != int(self.config["rounds"]) - 1:
                    self._train_augustus(i + 2, str(self.dependency_input["fasta"]), out_gff)
            # Move any augustus-generated config stuff
            self._handle_config_output()
            # Rename final file
            os.replace(out_gff, str(self.output["ab-gff3"]))

        def _augustus(self, species: str, _round: int, _file: str, _last: bool = False) -> str:
            """ Run augustus training round

            :param species: Species string
            :param _round: Training round number
            :param _file: FASTA file
            :param _last: Is last training round
            :return: Path to output gff3 file
            """
            out_gff = os.path.join(
                self.wdir, AugustusIter.Augustus.out_path(str(self.dependency_input["fasta"]), ".%i.gb" % _round)
            )
            # Chunk file predictions
            record_p = SeqIO.parse(_file, "fasta")
            progs = []
            out_files = []
            out_gffs = []
            for record in record_p:
                out_file_path = os.path.join(self.wdir, str(record.id) + ".fna")
                _out_gff = out_gff + "-" + str(record.id)
                out_files.append(out_file_path)
                out_gffs.append(_out_gff)
                SeqIO.write([record], out_file_path, "fasta")
                # Run prediction
                progs.append(
                    self.program[
                        "--codingseq=on",
                        "--stopCodonExcludedFromCDS=false",
                        "--species=%s" % species,
                        "--outfile=%s" % _out_gff,
                        ("--gff3=on" if _last else "--gff3=off"),
                        out_file_path,
                    ]
                )
            self.batch(progs)
            touch(out_gff + ".tmp")
            for out_g in out_gffs:
                (self.local["cat"][out_g] >> out_gff + ".a.tmp")()
            # Combine files
            self.single(
                self.local["gffread"]["-o", out_gff + ".tmp", "-F", "-G", "--keep-comments", out_gff + ".a.tmp"]
            )
            # Make ids unique
            self._make_unique(out_gff)
            # Remove intermediary files
            all([os.remove(_file) for _file in out_files])
            all([os.remove(_file) for _file in out_gffs])
            return out_gff

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

        @program_catch
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
            out_gb = os.path.join(self.wdir, AugustusIter.Augustus.out_path(_file, ".%i.gb" % _round))
            self.single(
                self.local["gff2gbSmallDNA.pl"][
                    out_gff,
                    _file,
                    "1000",
                    out_gb
                ]
            )

            species_config_prefix = self.record_id + str(_round)
            # Write new species config file
            self.single(
                self.local["new_species.pl"][
                    "--species=%s" % species_config_prefix,
                    out_gb
                ]
            )
            # Run training
            self.single(
                self.local["etraining"][
                    "--species=%s" % species_config_prefix,
                    out_gb
                ]
            )
            return out_gff

        @staticmethod
        def _make_unique(out_gff: str):
            """ Make gene identifiers in GFF3 file unique

            :param out_gff: Path to output gff3 file (less the .tmp identifier)
            """
            gff_fp = open(out_gff + ".tmp", "r")
            out_fp = open(out_gff, "w")
            i = 1
            line = next(gff_fp)
            while True:
                if line.startswith("#"):
                    out_fp.write(line)
                else:
                    line = line.split("\t")
                    if line[2] == "transcript":
                        out_fp.write("\t".join((
                            line[0],
                            "augustus",
                            *line[2:-1],
                            "ID=gene%i\n" % i
                        )))
                        try:
                            line = next(gff_fp).split("\t")
                        except StopIteration:
                            break
                        while line[2] != "transcript":
                            out_fp.write("\t".join((
                                line[0],
                                "augustus",
                                *line[2:-1],
                                "Parent=gene%i\n" % i
                            )))
                            try:
                                line = next(gff_fp).split("\t")
                            except StopIteration:
                                break
                        i += 1
                try:
                    line = next(gff_fp)
                except StopIteration:
                    break
            out_fp.close()

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
            return augustus_ids_dict[found_taxa.most_common()[0][0]]

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(AugustusIter.Augustus, AugustusIter.name, *args, **kwargs)
