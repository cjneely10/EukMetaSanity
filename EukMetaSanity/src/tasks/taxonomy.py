#!/usr/bin/env python3
from plumbum import local
from typing import Dict, List
from EukMetaSanity.src.utils.path_manager import PathManager
from EukMetaSanity.src.tasks.task_class import Task, TaskList
from EukMetaSanity.src.utils.config_manager import ConfigManager

mmseqs = local["mmseqs"]

_TAXONOMY_NAME = "taxonomy"


class TaxonomyIter(TaskList):
    class Taxonomy(Task):
        def __init__(self, input_paths_dict: Dict[str, str], cfg: ConfigManager, pm: PathManager, record_id: str):
            super().__init__(input_paths_dict, cfg, pm, record_id, _TAXONOMY_NAME, ["-in", "-db"])

        def run(self):
            # Generate genome directory
            # Call superclass run method
            self.required_data = ["-out"]
            self.output_paths_dict = {"-out": self.record_id}
            super().run()

        def results(self):
            # Call superclass results method
            return super().results()

    def __init__(self, input_paths_dict: List[Dict[str, str]], cfg: ConfigManager, pm: PathManager,
                 record_ids: List[str]):
        super().__init__(
            [TaxonomyIter.Taxonomy(input_path_dict, cfg, pm, record_id)
             for input_path_dict, record_id in zip(input_paths_dict, record_ids)]
        )

    def run(self):
        super().run()

    def results(self):
        return super().results()
