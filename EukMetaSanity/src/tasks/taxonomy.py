import os
import logging
from typing import Dict, List, Tuple
from EukMetaSanity.src.utils.data import Data
from EukMetaSanity.src.utils.path_manager import PathManager
from plumbum.commands.processes import ProcessExecutionError
from EukMetaSanity.src.tasks.task_class import Task, TaskList
from EukMetaSanity.src.utils.config_manager import ConfigManager

"""
Determine the taxonomy of the Eukaryotic MAG

"""


class TaxonomyIter(TaskList):
    class Taxonomy(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def run(self) -> None:
            seq_db = os.path.join(self.wdir, self.record_id + "_db")
            results_file = os.path.join(self.wdir, self.record_id + "-tax-report.txt")
            # Expected output
            self.output = {Data.OUT: [
                results_file,  # Taxonomic results for ab initio prediction
                self.input[Data.IN][0],  # Input FASTA file for repeat masking
                seq_db,  # MMseqs database for use in metaeuk or repeat masking
            ]}
            super().run()

        def results(self) -> Dict[str, List[str]]:
            return super().results()

        def run_1(self):
            tax_db = os.path.join(self.wdir, self.record_id + "-tax_db")
            seq_db = self.output[Data.OUT][2]
            try:
                # Create sequence database
                super().log_and_run(
                    self.program[
                        "createdb",
                        self.input[Data.IN],  # Input FASTA file
                        seq_db,  # Output FASTA sequence db
                    ],
                    self.mode
                )
                # Run taxonomy search
                super().log_and_run(
                    self.program[
                        "taxonomy",
                        seq_db,  # Input FASTA sequence db
                        self.input[Data.ACCESS],  # Input OrthoDB
                        tax_db,  # Output tax db
                        os.path.join(self.wdir, "tmp"),
                        (*self.cfg.get_added_flags(Data().taxonomy()[0])),
                        "--threads", self.threads,
                    ],
                    self.mode
                )
                # Output results
                super().log_and_run(
                    self.program[
                        "taxonomyreport",
                        self.input[Data.ACCESS],  # Input OrthoDB
                        tax_db,  # Input tax db
                        self.output[Data.OUT][0]  # Output results file
                    ],
                    self.mode
                )
            except ProcessExecutionError as e:
                logging.info(e)

    def __init__(self, input_paths: List[List[str]], cfg: ConfigManager, pm: PathManager,
                 record_ids: List[str], mode: int):
        dt = Data()
        super().__init__(TaxonomyIter.Taxonomy, input_paths, record_ids, dt.taxonomy, cfg, pm, mode,
                         {Data.ACCESS: [dt.taxonomy()[1]]})

    def run(self):
        super().run()

    def output(self) -> Tuple[List[List[str]], ConfigManager, PathManager, List[str], int]:
        return super().output()


if __name__ == "__main__":
    pass
