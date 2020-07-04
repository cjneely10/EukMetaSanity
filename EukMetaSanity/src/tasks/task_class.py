from abc import ABC, abstractmethod
from numba import types, jitclass, jit
from EukMetaSanity.src.utils.path_manager import PathManager


@jitclass((
    ("_input_paths", types.Tuple(types.string)),
    ("output_paths", types.Tuple(types.string)),
    ("_threads_pw", types.int8),
    ("_workers", types.int8),
    ("_pm", PathManager),

))
class Task(ABC):
    def __init__(self, input_fasta_list, cfg, pm, dict_name):
        self._input_paths = input_fasta_list
        self.output_paths = None
        self._threads_pw = int(cfg.config.get(dict_name, {}).get("THREADS", 1)),
        self._workers = int(cfg.config.get(dict_name, {}).get("WORKERS", 1)),
        self._pm = pm
        super().__init__()

    @property
    def threads(self):
        return self._threads_pw

    @property
    def workers(self):
        return self._workers

    @property
    def pm(self):
        return self._pm

    @abstractmethod
    def run(self):
        assert self.output_paths is not None

    @abstractmethod
    @jit(types.ListType(types.string)())
    def results(self):
        if self.output_paths is None:
            self.run()
        return self.output_paths


if __name__ == "__main__":
    pass
