from plumbum import local
from numba import types, jit
from EukMetaSanity.tasks.task_class import Task
from EukMetaSanity.src.utils.path_manager import PathManager

mmseqs = local["mmseqs"]


class Taxonomy(Task):
    @jit(types.void(PathManager))
    def __init__(self, _pm):
        super().__init__(
            ()
        )

    def results(self):
        pass

    def run(self):
        pass
