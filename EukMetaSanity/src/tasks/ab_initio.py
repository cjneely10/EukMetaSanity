import os
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.src.tasks.taxonomy import TaxonomyIter
from EukMetaSanity.src.scripts.fastagff3_to_gb import write_genbank


class AbInitioIter(TaskList):
    class AbInitio(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # Only looking for final trained ab initio prediction
            self.output = [
                os.path.join(self.wdir, self.record_id + ".gff3"),  # Output gff3 ab initio predictions, final round
                self.input[1],  # Forward masked mmseqs-db to initial evidence step
                self.input[2],  # Original file
            ]

        def run(self) -> None:
            super().run()

        def run_1(self):
            # Call protocol method
            getattr(self, self.protocol)()

        def augustus(self):
            # Initial training based on best species from taxonomy search
            out, _file = self._augustus(self._augustus_tax_ident(), 1, self.input[0])
            # Remaining rounds of re-training on generated predictions
            for i in range(int(self.rounds)):
                out, _file = self._augustus(self.record_id + str(i + 2), i + 2, _file)
            if self.mode == 1:
                os.replace(out, self.output[0])

        @program_catch
        def _augustus_tax_ident(self) -> str:
            tax_db = os.path.join(self.wdir, self.record_id + "-tax_db")
            seq_db = self.input[1]
            # Run taxonomy search
            self.log_and_run(
                self.program_mmseqs[
                    "taxonomy",
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
                    "taxonomyreport",
                    self.data,  # Input augustus-db
                    tax_db,  # Input tax db
                    tax_db + ".taxreport"  # Output results file
                ]
            )
            # Return optimal taxonomy
            return TaxonomyIter.Taxonomy.get_taxonomy(tax_db + ".taxreport", float(self.cutoff))[0]

        @program_catch
        def _augustus(self, species: str, _round: int, _file: str):
            out_gff = AbInitioIter.AbInitio._out_path(self.input[1], ".%i.gff3" % _round)
            # Run prediction
            self.log_and_run(
                self.program[
                    "--codingseq=on",
                    "--stopCodonExcludedFromCDS=true",
                    "--species=%s" % species,
                    _file,
                    "--outfile", out_gff
                ]
            )
            # Parse to genbank
            out_gb = AbInitioIter.AbInitio._out_path(self.input[1], ".%i.gb" % _round)
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
            return out_gff, out_gb

        @program_catch
        def gmes(self):
            self.log_and_run(
                self.program[
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

        @staticmethod
        def _out_path(_file_name: str, _ext: str) -> str:
            return os.path.splitext(_file_name)[0] + _ext

    def __init__(self, *args, **kwargs):
        super().__init__(AbInitioIter.AbInitio, "abinitio", *args, **kwargs)


if __name__ == "__main__":
    pass
