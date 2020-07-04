from abc import ABC, abstractmethod
from numba import types, jitclass, jit


@jitclass((
    ("_input_paths", types.Tuple(types.string)),
    ("output_paths", types.Tuple(types.string)),
    ("_threads_pw", types.int8),
    ("_workers", types.int8),

))
class Task(ABC):
    def __init__(self, input_paths, threads_pw, workers):
        self._input_paths = input_paths
        self.output_paths = []
        self._threads_pw = threads_pw
        self._workers = workers
        super().__init__()

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    @jit(types.ListType(types.string)())
    def results(self):
        if len(self.output_paths) == 0:
            self.run()
        return self.output_paths


if __name__ == "__main__":
    pass
