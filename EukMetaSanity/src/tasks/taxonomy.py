#!/usr/bin/env python3
import os
from plumbum import local
from typing import Dict, List, Tuple
from EukMetaSanity.src.utils.data import Data
from EukMetaSanity.src.utils.helpers import log_and_run
from EukMetaSanity.src.utils.path_manager import PathManager
from plumbum.commands.processes import ProcessExecutionError
from EukMetaSanity.src.tasks.task_class import Task, TaskList
from EukMetaSanity.src.utils.config_manager import ConfigManager

"""
Determine the taxonomy of the Eukaryotic MAG

"""

mmseqs = local["mmseqs"]


class TaxonomyIter(TaskList):
    class Taxonomy(Task):
        def __init__(self, input_path_dict: Dict[str, str], cfg: ConfigManager, pm: PathManager, record_id: str):
            super().__init__(input_path_dict, cfg, pm, record_id, Data().taxonomy()[0], [Data.IN, Data.ACCESS])
            # Set required data
            self.required_data = [Data.OUT]

        def run(self):
            # Call superclass run method
            super().run()

        def results(self) -> Dict[str, str]:
            # Call superclass results method
            return super().results()

        def parse_output(self, output_files: List[str]) -> List[Dict[str, str]]:
            pass

        def run_1(self):
            name, ident = Data().taxonomy()
            # Create sequence database
            seq_db = os.path.join(self.wdir, self.record_id + "_db")
            tax_db = os.path.join(self.wdir, self.record_id + "-tax_db")
            results_file = os.path.join(self.wdir, self.record_id + "-tax-report.txt")
            try:
                log_and_run(
                    mmseqs[
                        "createdb",
                        self.input[Data.IN],  # Input FASTA file
                        seq_db,  # Output FASTA sequence db
                    ]
                )
                # Run taxonomy search
                log_and_run(
                    mmseqs[
                        "taxonomy",
                        seq_db,  # Input FASTA sequence db
                        self.input[Data.ACCESS],  # Input OrthoDB
                        tax_db,  # Output tax db
                        os.path.join(self.wdir, "tmp"),
                        (*self.cfg.get_added_flags(name)),
                        "--threads", self.threads,
                    ]
                )
                # Output results
                log_and_run(
                    mmseqs[
                        "taxonomyreport",
                        self.input[Data.ACCESS],  # Input OrthoDB
                        tax_db,  # Input tax db
                        results_file  # Output results file
                    ]
                )
            except ProcessExecutionError as e:
                print(e)
            # DB path
            self.output_paths_dict = {Data.OUT: results_file}

    def __init__(self, input_paths: List[str], cfg: ConfigManager, pm: PathManager,
                 record_ids: List[str]):
        name, ident = Data().taxonomy()
        workers = int(cfg.config.get(name, ConfigManager.WORKERS))
        super().__init__(
            # List of tasks
            [
                TaxonomyIter.Taxonomy(
                    {Data.IN: input_path, Data.ACCESS: cfg.config[ConfigManager.DATA][ident]},
                    cfg,
                    pm,
                    record_id
                )
                for input_path, record_id in zip(input_paths, record_ids)
            ],
            # Logging statement
            "Running mmseqs to identify taxonomy using %i workers and %i threads per worker" % (
                workers,
                int(cfg.config.get(name, ConfigManager.THREADS)),
            ),
            # Dask workers
            workers,
            cfg,
            pm,
        )

    def run(self):
        super().run()

    def results(self):
        return super().results()

    def output(self) -> Tuple[List[str], ConfigManager, PathManager, List[str]]:
        # Run task list
        return (
            [result[Data.OUT] for result in self.results()],
            self.cfg,
            self.pm,
            [task.record_id for task in self.tasks]
        )


if __name__ == "__main__":
    pass
