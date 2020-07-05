import os
import logging
from typing import Dict, List, Tuple
from EukMetaSanity.src.utils.data import Data
from EukMetaSanity.src.utils.helpers import log_and_run
from EukMetaSanity.src.utils.path_manager import PathManager
from plumbum.commands.processes import ProcessExecutionError
from EukMetaSanity.src.tasks.task_class import Task, TaskList
from EukMetaSanity.src.utils.config_manager import ConfigManager

"""
Model the repeated regions of a FASTA sequence

"""


class RepeatsIter(TaskList):
    class Repeats(Task):
        def __init__(self, input_path_dict: Dict[str, List[str]], cfg: ConfigManager, pm: PathManager, record_id: str,
                     mode: int):
            # Needs an input file only
            super().__init__(input_path_dict, cfg, pm, record_id, Data().taxonomy()[0], [Data.IN], mode)

        def run(self) -> None:
            super().run()

        def results(self) -> Dict[str, List[str]]:
            return super().results()

        def run_1(self):
            try:
                pass
            except ProcessExecutionError as e:
                logging.info(e)

    def __init__(self, input_paths: List[str], cfg: ConfigManager, pm: PathManager,
                 record_ids: List[str], mode: int):
        name, ident = Data().repeat_modeling()
        workers = int(cfg.config.get(name, ConfigManager.WORKERS))
        # Parse for taxonomy identification
        super().__init__(
            [
            ],
            "Running mmseqs to identify taxonomy using %i workers and %i threads per worker" % (
                workers,
                int(cfg.config.get(name, ConfigManager.THREADS)),
            ),
            workers,
            cfg,
            pm,
            mode
        )

    def simple(self):
        pass

    def full(self):
        pass

    def run(self):
        super().run()

    def output(self) -> Tuple[List[List[str]], ConfigManager, PathManager, List[str], int]:
        return super().output()


if __name__ == "__main__":
    pass
