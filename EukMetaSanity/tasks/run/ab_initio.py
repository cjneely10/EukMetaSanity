import os
import shutil
from Bio import SeqIO
from pathlib import Path
from typing import Tuple, List
from collections import Counter
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.utils.helpers import augustus_taxon_ids

"""
Perform ab initio gene identification using either Augustus or GeneMark

"""


class FailedAugustusIdentification(ChildProcessError):
    pass


class AbInitioIter(TaskList):
    class AbInitio(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # Only looking for final trained ab initio prediction
            self.output = [
                os.path.join(self.wdir, self.record_id + ".gff3"),  # Output gff3 ab initio predictions, final round
                self.input[1],  # Forward masked mmseqs-db to initial evidence step
                self.input[2],  # Original file,
                self.input[3],  # Tax file
                self.input[0],  # Repeat-masked FASTA file
                self.input[5],  # Summarized repeats file
                self.input[6],  # Repeats gff3 file
            ]
            self.rounds = int(self.rounds)

        def run(self) -> None:
            super().run()

        def run_1(self):
            # Call protocol method
            getattr(self, self.protocol)()

        def augustus(self):
            # Initial training based on best species from taxonomy search
            out_gff = self._augustus(self._augustus_tax_ident(), 1, self.input[0])
            self._train_augustus(1, self.input[0], out_gff)
            # Remaining rounds of re-training on generated predictions
            for i in range(self.rounds):
                _last = False
                if i == self.rounds - 1:
                    _last = True
                out_gff = self._augustus(self.record_id + str(i + 1), i + 2, self.input[0], _last)
                if i != self.rounds - 1:
                    self._train_augustus(i + 2, self.input[0], out_gff)
            # Move any augustus-generated config stuff
            self._handle_config_output()
            # Rename final file
            os.replace(out_gff, os.path.join(self.wdir, self.record_id + ".gff3"))

        def _augustus(self, species: str, _round: int, _file: str, _last: bool = False):
            out_gff = os.path.join(self.wdir, AbInitioIter.AbInitio._out_path(self.input[1], ".%i.gff3" % _round))
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
                    self.program_augustus[
                        "--codingseq=on",
                        "--stopCodonExcludedFromCDS=false",
                        "--species=%s" % species,
                        "--outfile=%s" % _out_gff,
                        ("--gff3=on" if _last else "--gff3=off"),
                        out_file_path,
                    ]
                )
            self.batch(progs)
            # Combine files
            (
                self.local["cat"][out_gffs] |
                self.local["gffread"]["-o", out_gff + ".tmp", "-F", "-G", "--keep-comments"]
            )()
            # Make ids unique
            self._make_unique(out_gff)
            # Remove intermediary files
            all([os.remove(_file) for _file in out_files])
            all([os.remove(_file) for _file in out_gffs])
            return out_gff

        @program_catch
        def _augustus_tax_ident(self) -> str:
            tax_db = os.path.join(self.wdir, self.record_id + "-augustus_db")
            seq_db = self.input[4]
            if not os.path.exists(tax_db + ".m8"):
                # Run taxonomy search
                self.log_and_run(
                    self.program_mmseqs[
                        "linsearch",
                        seq_db,  # Input FASTA sequence db
                        self.data,  # Input augustus-db
                        tax_db,  # Output tax db
                        os.path.join(self.wdir, "tmp"),
                        (*self.added_flags),
                        "--threads", self.threads,
                    ]
                )
                # Output results
                self.log_and_run(
                    self.program_mmseqs[
                        "convertalis",
                        seq_db,  # Input FASTA sequence db
                        self.data,  # Input augustus-db
                        tax_db,  # Input tax db
                        tax_db + ".m8",  # Output results file
                        "--threads", self.threads,
                        "--format-output", "query,target,pident,taxid,taxname,taxlineage",
                    ]
                )
            # Return optimal taxonomy
            augustus_ids_dict = augustus_taxon_ids()
            found_taxa = Counter()
            with open(tax_db + ".m8", "r") as R:
                for line in R:
                    line = line.rstrip("\r\n").split()
                    # Count those that pass the user-defined cutoff value
                    if line[3] in augustus_ids_dict.keys() and float(line[2]) >= (float(self.cutoff) / 100.):
                        found_taxa[line[3]] += 1
            return augustus_ids_dict[found_taxa.most_common()[0][0]]

        @program_catch
        def _train_augustus(self, _round: int, _file: str, out_gff: str):
            # Remove old training directory, if needed
            config_dir = os.path.join(
                os.path.dirname(os.path.dirname(Path(str(self.program_augustus)).resolve())),
                "config", "species", self.record_id + str(_round)
            )
            if os.path.exists(config_dir):
                shutil.rmtree(config_dir)
            # Parse to genbank
            out_gb = os.path.join(self.wdir, AbInitioIter.AbInitio._out_path(_file, ".%i.gb" % _round))
            self.log_and_run(
                self.local["gff2gbSmallDNA.pl"][
                    out_gff,
                    _file,
                    "1000",
                    out_gb
                ]
            )

            species_config_prefix = self.record_id + str(_round)
            # Write new species config file
            self.log_and_run(
                self.program_new_species_pl[
                    "--species=%s" % species_config_prefix,
                    out_gb
                ]
            )
            # Run training
            self.log_and_run(
                self.program_etraining[
                    "--species=%s" % species_config_prefix,
                    out_gb
                ]
            )
            return out_gff

        @staticmethod
        def _out_path(_file_name: str, _ext: str) -> str:
            return os.path.basename(os.path.splitext(_file_name)[0]) + _ext

        @staticmethod
        def _make_unique(out_gff):
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
                            *line[0:-1],
                            "ID=gene%i\n" % i
                        )))
                        try:
                            line = next(gff_fp).split("\t")
                        except StopIteration:
                            break
                        while line[2] != "transcript":
                            out_fp.write("\t".join((
                                *line[0:-1],
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

        def _handle_config_output(self):
            # Move the augustus training folders to their wdir folders
            config_dir = os.path.join(
                os.path.dirname(os.path.dirname(Path(str(self.program_augustus)).resolve())),
                "config", "species"
            )
            for i in range(1, int(self.rounds)):
                shutil.move(
                    os.path.join(config_dir, self.record_id + str(i + 1)),
                    self.wdir
                )

        @program_catch
        def gmes(self):
            assert os.path.exists(self.gmes_cfg)
            # Copy base config file to working dir
            def_cfg = os.path.join(self.wdir, "gmes.default.cfg")
            new_cfg = os.path.join(self.wdir, "gmes.cfg")
            self.local["cp"][self.gmes_cfg, def_cfg]()
            # Update working directory path in cfg file
            AbInitioIter.AbInitio.update_cfg(
                def_cfg,
                [
                    ("def_cfg", new_cfg),
                    ("work_dir", self.wdir),
                ],
                new_cfg
            )
            # Run GeneMark with updated config file
            self.log_and_run(
                self.program_gmes[
                    "--sequence", self.input[0],
                    "--ES", "--cores", self.threads, (*self.added_flags),
                    "--usr_cfg", new_cfg
                ]
            )
            # Move program to match required output name
            if self.mode == 1:
                os.replace(
                    os.path.join(self.wdir, "genemark.gtf"),
                    self.output[0],
                )

        @staticmethod
        def update_cfg(in_path: str, replace_tuple: List[Tuple[str, str]], out_path: str):
            in_fp = open(in_path, "r")
            out_fp = open(out_path, "w")
            for line in in_fp:
                for _tuple in replace_tuple:
                    if _tuple[0] not in line:
                        out_fp.write(line)
                    else:
                        out_fp.write("  %s:  %s" % _tuple)
            out_fp.close()

    def __init__(self, *args, **kwargs):
        super().__init__(AbInitioIter.AbInitio, "abinitio", *args, **kwargs)


if __name__ == "__main__":
    pass
