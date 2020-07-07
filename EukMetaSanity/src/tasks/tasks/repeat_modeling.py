import os
from EukMetaSanity.src.utils.data import Data
from EukMetaSanity.src.utils.config_manager import ConfigManager
from EukMetaSanity.src.tasks.task_class import Task, TaskList, program_catch

"""
Model the repeated regions of a FASTA sequence

"""


class RepeatsIter(TaskList):
    class Repeats(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def run(self) -> None:
            masked_db_path = os.path.join(self.wdir, self.record_id + "-mask_db")
            masked_fa_path = masked_db_path[:-3] + ".fna"
            self.output = {Data.Type.OUT: [
                self.input[Data.Type.IN][0],  # Taxonomic results for ab initio prediction
                masked_fa_path,  # Input FASTA file for ab initio
                masked_db_path,  # MMseqs database for use in metaeuk
            ]}
            super().run()

        def run_1(self):
            # Call protocol method
            getattr(self, self.config[ConfigManager.PROTOCOL])()

        # Simple repeat masking using mmseqs
        @program_catch
        def simple(self):
            masked_db_path = self.output[Data.Type.OUT][2]
            # Generate the masked sequence
            self.log_and_run(
                self.program[
                    "masksequence",
                    self.input[Data.Type.IN][2],
                    masked_db_path,
                    "--threads", self.threads,
                ]
            )
            # Output as FASTA file
            self.log_and_run(
                self.program[
                    "convert2fasta",
                    masked_db_path,
                    self.output[Data.Type.OUT][1],
                ]
            )

        # Complete masking using RepeatModeler/Masker
        def full(self):
            self._model()
            self._mask()

        @program_catch
        def _model(self):
            # Build database
            self.log_and_run(
                self.program[
                    "-name", self.output[Data.Type.OUT][2],
                    (*self.added_flags),
                    self.input[Data.Type.IN][1],
                ]
            )
            # Run RepeatModeler
            self.log_and_run(
                self.program2[
                    "-pa", self.threads,
                    (*self.added_flags),
                    "-database", self.output[Data.Type.OUT][2],
                ]
            )

        @program_catch
        def _mask(self):
            pass

    def __init__(self, *args, **kwargs):
        super().__init__(RepeatsIter.Repeats, "repeats", *args, **kwargs)


if __name__ == "__main__":
    pass
