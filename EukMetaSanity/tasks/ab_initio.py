import os
import shutil
from Bio import SeqIO
from pathlib import Path
from collections import Counter
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.utils.helpers import augustus_taxon_ids
from EukMetaSanity.scripts.fastagff3_to_gb import write_genbank


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
            ]

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
            for i in range(int(self.rounds)):
                out_gff = self._augustus(self.record_id + str(i + 1), i + 2, self.input[0])
                self._train_augustus(1, self.input[0], out_gff)
            # Move any augustus-generated config stuff
            self._handle_config_output()
            # Rename final file
            shutil.copy(out_gff, os.path.join(self.wdir, self.record_id + ".gff3"))

        def _augustus(self, species: str, _round: int, _file: str):
            out_gff = os.path.join(self.wdir, AbInitioIter.AbInitio._out_path(self.input[1], ".%i.gff3" % _round))
            # Chunk file predictions
            record_p = SeqIO.parse(_file, "fasta")
            for record in record_p:
                out_file_path = str(record.id) + ".fna"
                SeqIO.write([record], out_file_path, "fasta")
                # Run prediction
                self.log_and_run(
                    self.program_augustus[
                        "--codingseq=on",
                        "--stopCodonExcludedFromCDS=true",
                        "--species=%s" % species,
                        "--outfile=%s" % out_gff + "-" + str(record.id),
                        out_file_path,
                    ]
                )
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
                    if line[3] in augustus_ids_dict.keys() and float(line[2]) > (float(self.cutoff) / 100.):
                        found_taxa[line[3]] += 1
            return augustus_ids_dict[found_taxa.most_common()[0][0]]

        @program_catch
        def _train_augustus(self, _round: int, _file: str, out_gff: str):
            # Parse to genbank
            out_gb = os.path.join(self.wdir, AbInitioIter.AbInitio._out_path(self.input[1], ".%i.gb" % _round))
            write_genbank(
                _file,
                out_gff,
                out_gb
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

        def _handle_config_output(self):
            # Move the augustus training folders to their wdir folders
            config_dir = os.path.join(
                os.path.dirname(os.path.dirname(Path(str(self.program_augustus)).resolve())),
                "config"
            )
            for i in range(int(self.rounds)):
                shutil.move(
                    os.path.join(config_dir, self.record_id + str(i + 2)),
                    self.wdir
                )

        @program_catch
        def gmes(self):
            self.log_and_run(
                self.program_gmes[
                    "--sequence", self.input[0],
                    "--ES", "--cores", self.threads, (*self.added_flags)
                ]
            )
            # Move program to match required output name
            if self.mode == 1:
                os.replace(
                    os.path.join(self.pm.get_dir(self.record_id, self.name), "genemark.gtf"),
                    self.output[0],
                )

    def __init__(self, *args, **kwargs):
        super().__init__(AbInitioIter.AbInitio, "abinitio", *args, **kwargs)


if __name__ == "__main__":
    pass
