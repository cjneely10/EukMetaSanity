#!/usr/bin/env python3
from plumbum import local
from typing import Dict, List
from dask.distributed import Client, wait
from EukMetaSanity.src.utils.path_manager import PathManager
from EukMetaSanity.src.tasks.task_class import Task, TaskList
from EukMetaSanity.src.utils.config_manager import ConfigManager

mmseqs = local["mmseqs"]

_TAXONOMY_NAME = "taxonomy"


class TaxonomyIter(TaskList):
    INFO = "Running mmseqs to identify taxonomy"

    class Taxonomy(Task):
        def __init__(self, input_paths_dict: Dict[str, str], cfg: ConfigManager, pm: PathManager, record_id: str):
            super().__init__(input_paths_dict, cfg, pm, record_id, _TAXONOMY_NAME, ["-in", "-db"])

        def run(self):
            # Set required data
            self.required_data = ["-out"]
            # Logic for mmseqs

            # Parse for output
            self.output_paths_dict = {"-out": self.record_id}
            # Call superclass run method
            super().run()

        def results(self):
            # Call superclass results method
            return super().results()

        def parse_output(self, output_files: List[str]) -> List[Dict[str, str]]:
            pass

    def __init__(self, input_paths_dict: List[Dict[str, str]], cfg: ConfigManager, pm: PathManager,
                 record_ids: List[str]):
        super().__init__(
            [TaxonomyIter.Taxonomy(input_path_dict, cfg, pm, record_id)
             for input_path_dict, record_id in zip(input_paths_dict, record_ids)],
            TaxonomyIter.INFO
        )

    def run(self):
        super().run()

    def results(self):
        return super().results()
