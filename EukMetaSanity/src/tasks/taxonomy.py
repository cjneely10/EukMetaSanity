#!/usr/bin/env python3
from plumbum import local
from EukMetaSanity.src.tasks.task_class import Task

mmseqs = local["mmseqs"]

_TAXONOMY_NAME = "taxonomy"


class Taxonomy(Task):
    def __init__(self, input_paths_dict, cfg, pm, record_id):
        super().__init__(input_paths_dict, cfg, pm, record_id, _TAXONOMY_NAME)
        print(self.wdir)

    def run(self):
        # Generate genome directory
        # Call superclass run method
        super().run()

    def results(self):
        # Call superclass results method
        return super().results()
