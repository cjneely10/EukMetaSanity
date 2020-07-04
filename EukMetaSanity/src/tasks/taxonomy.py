#!/usr/bin/env python3
from plumbum import local
from typing import Dict, List
from dask.distributed import Client, wait
from EukMetaSanity.src.utils.data import Data
from EukMetaSanity.src.utils.path_manager import PathManager
from EukMetaSanity.src.tasks.task_class import Task, TaskList
from EukMetaSanity.src.utils.config_manager import ConfigManager

mmseqs = local["mmseqs"]


class TaxonomyIter(TaskList):
    class Taxonomy(Task):
        def __init__(self, input_path_dict: Dict[str, str], cfg: ConfigManager, pm: PathManager, record_id: str):
            super().__init__(input_path_dict, cfg, pm, record_id, "taxonomy", [Data.IN, Data.ACCESS])

        def run(self):
            # Set required data
            self.required_data = [Data.OUT]
            # Logic for mmseqs

            # Parse for output
            self.output_paths_dict = {Data.OUT: "MAGS"}
            # Call superclass run method
            super().run()

        def results(self) -> Dict[str, str]:
            # Call superclass results method
            return super().results()

        def parse_output(self, output_files: List[str]) -> List[Dict[str, str]]:
            pass

        def run_mmseqs(self):
            # Create sequence database
            # Run taxonomy search
            # Output results
            pass

    def __init__(self, input_paths: List[str], cfg: ConfigManager, pm: PathManager,
                 record_ids: List[str]):
        super().__init__(
            [
                TaxonomyIter.Taxonomy(
                    {Data.IN: input_path, Data.ACCESS: cfg.config[Data.DATA][Data().taxonomy()]},
                    cfg,
                    pm,
                    record_id
                )
                for input_path, record_id in zip(input_paths, record_ids)
            ],
            "Running mmseqs to identify taxonomy"
        )

    def run(self):
        super().run()

    def results(self):
        return super().results()
