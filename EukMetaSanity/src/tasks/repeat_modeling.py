import os
import logging
from io import TextIOWrapper
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
            super().__init__(input_path_dict, cfg, pm, record_id, Data().repeat_modeling()[0], [Data.IN], mode)

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
        # Get name of config section and required data
        name, _ = Data().repeat_modeling()
        # Ensure protocol is passed and valid
        assert ConfigManager.PROTOCOL in cfg.config[name].keys()
        protocol = cfg.config.get(name, ConfigManager.PROTOCOL)
        assert protocol in dir(self)
        workers = int(cfg.config.get(name, ConfigManager.WORKERS))
        # Call protocol function to geenrate list for superclass initialization
        super().__init__(
            getattr(self, protocol)(input_paths, record_ids, cfg, pm, mode),
            "Running %s repeat identification using %i workers and %i threads per worker" % (
                protocol,
                workers,
                int(cfg.config.get(name, ConfigManager.THREADS)),
            ),
            workers,
            cfg,
            pm,
            mode
        )

    # Simple repeat annotation initialization using mmseqs
    def simple(self, inputs: List[str], r_ids: List[str], cfg: ConfigManager, pm: PathManager, mode: int) -> List[Task]:
        return [
            RepeatsIter.Repeats(
                {},
                cfg,
                pm,
                record_id,
                mode
            )
            for input_path, record_id in zip(inputs, r_ids)
        ]

    # Full repeat annotation using RepeatModeler/RepeatMasker
    def full(self, inputs: List[str], r_ids: List[str], cfg: ConfigManager, pm: PathManager, mode: int) -> List[Task]:
        pass

    def run(self):
        super().run()

    def output(self) -> Tuple[List[List[str]], ConfigManager, PathManager, List[str], int]:
        return super().output()


if __name__ == "__main__":
    pass
