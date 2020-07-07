import os

from EukMetaSanity.src.tasks.tasks.taxonomy import TaxonomyIter
from EukMetaSanity.src.utils.data import Data
from EukMetaSanity.bin.fastagff3_to_gb import write_genbank
from EukMetaSanity.src.utils.config_manager import ConfigManager
from EukMetaSanity.src.tasks.task_class import TaskList, Task, program_catch


class AbInitioIter(TaskList):
    class AbInitio(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def run(self) -> None:
            # Only looking for final trained ab initio prediction
            self.output = {Data.Type.OUT: [
                os.path.join(self.wdir, self.record_id + ".gff3"),  # Output gff3 ab initio predictions, final round
            ]}
            super().run()

        def run_1(self):
            # Call protocol method
            getattr(self, self.config[ConfigManager.PROTOCOL])()

        def augustus(self):
            out = self._augustus(self._augustus_tax_ident(), 1)
            for i in range(int(self.config[ConfigManager.ROUNDS]) - 1):
                out = self._augustus(self.record_id + str(i + 2), i + 2)
            if self.mode == 1:
                os.replace(out, self.output[Data.Type.OUT][0])

        @program_catch
        # TODO: Provide implementation for searching for optimal augustus species
        def _augustus_tax_ident(self) -> str:
            tax_db = os.path.join(self.wdir, self.record_id + "-tax_db")
            seq_db = self.output[Data.Type.OUT][2]
            # Run taxonomy search
            self.log_and_run(
                self.program[
                    "taxonomy",
                    seq_db,  # Input FASTA sequence db
                    self.input[Data.Type.ACCESS][0],  # Input OrthoDB
                    tax_db,  # Output tax db
                    os.path.join(self.wdir, "tmp"),
                    (*self.added_flags),
                    "--threads", self.threads,
                ]
            )
            # Output results
            self.log_and_run(
                self.program[
                    "taxonomyreport",
                    self.input[Data.Type.ACCESS][0],  # Input OrthoDB
                    tax_db,  # Input tax db
                    tax_db + ".taxreport"  # Output results file
                ]
            )
            # Return optimal taxonomy
            return TaxonomyIter.Taxonomy.get_taxonomy(tax_db + ".taxreport")

        @program_catch
        def _augustus(self, species: str, _round: int):
            out_gff = AbInitioIter.AbInitio._out_path(self.input[Data.Type.IN][1], ".%i.gff3" % _round)
            # Run prediction
            self.log_and_run(
                self.program[
                    "--codingseq=on",
                    "--stopCodonExcludedFromCDS=true",
                    "--species=%s" % species,
                    self.input[Data.Type.IN][1],
                    "--outfile", out_gff
                ]
            )
            # Parse to genbank
            out_gb = AbInitioIter.AbInitio._out_path(self.input[Data.Type.IN][1], ".%i.gb" % _round)
            write_genbank(
                self.input[Data.Type.IN][1],
                out_gff,
                out_gb
            )
            species_config_prefix = self.record_id + str(_round)
            # Write new species config file
            self.log_and_run(
                self.program2[
                    "--species=%s" % species_config_prefix,
                    out_gb
                ]
            )
            # Run training
            self.log_and_run(
                self.program3[
                    "--species=%s" % species_config_prefix,
                    out_gb
                ]
            )
            return out_gff

        @program_catch
        def gmes(self):
            self.log_and_run(
                self.program[
                    "--sequence", self.input[Data.Type.IN][1],
                    "--ES", "--cores", self.threads, (*self.added_flags)
                ]
            )
            # Move program to match required output name
            if self.mode == 1:
                os.replace(
                    os.path.join(self.pm.get_dir(self.record_id, self.name), "genemark.gtf"),
                    self.output[Data.Type.OUT][0],
                )

        @staticmethod
        def _out_path(_file_name: str, _ext: str) -> str:
            _file_name = _file_name.split(".")
            return ".".join(_file_name[:-1]) + _ext

    def __init__(self, *args, **kwargs):
        super().__init__(AbInitioIter.AbInitio, "abinitio", *args, **kwargs)


if __name__ == "__main__":
    pass
