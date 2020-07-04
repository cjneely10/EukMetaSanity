#!/usr/bin/env python3
from plumbum import local
from numba import types, jit
from EukMetaSanity.src.tasks.task_class import Task
from EukMetaSanity.src.utils.config_manager import ConfigManager

mmseqs = local["mmseqs"]

_TAXONOMY_NAME = "TAXONOMY"


class Taxonomy(Task):
    @jit(types.void(types.Tuple(types.string), ConfigManager), nopython=True, cache=True)
    def __init__(self, input_fasta_list, pm, cfg):
        super().__init__(input_fasta_list, cfg, pm, _TAXONOMY_NAME)

    def run(self):
        # Generate genome directory
        # Call superclass run method
        super().run()

    def results(self):
        # Call superclass results method
        return super().results()
