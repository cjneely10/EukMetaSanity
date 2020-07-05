import os
import logging
from typing import Dict, List, Tuple
from EukMetaSanity.src.utils.data import Data
from EukMetaSanity.src.utils.path_manager import PathManager
from plumbum.commands.processes import ProcessExecutionError
from EukMetaSanity.src.tasks.task_class import Task, TaskList
from EukMetaSanity.src.utils.config_manager import ConfigManager

"""
Model the repeated regions of a FASTA sequence

"""


class RepeatsIter(TaskList):
    class Repeats(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def run(self) -> None:
            masked_db_path = os.path.join(self.wdir, self.record_id + "-mask_db")
            masked_fa_path = masked_db_path[:-3] + ".fna"
            self.output = {Data.OUT: [
                self.input[Data.IN][0],  # Taxonomic results for ab initio prediction
                masked_fa_path,  # Input FASTA file for ab initio
                masked_db_path,  # MMseqs database for use in metaeuk
            ]}
            super().run()

        def results(self) -> Dict[str, List[str]]:
            return super().results()

        def run_1(self):
            name = Data().repeat_modeling()[0]
            # Call protocol method
            getattr(self, self.cfg.config.get(name, ConfigManager.PROTOCOL))()

        # Simple repeat masking using mmseqs
        def simple(self):
            masked_db_path = self.output[Data.OUT][2]
            try:
                # Generate the masked sequence
                super().log_and_run(
                    self.program[
                        "masksequence",
                        self.input[Data.IN][2],
                        masked_db_path,
                        "--threads", self.threads,
                    ],
                    self.mode
                )
                # Output as FASTA file
                super().log_and_run(
                    self.program[
                        "convert2fasta",
                        masked_db_path,
                        self.output[Data.OUT][1],
                    ],
                    self.mode
                )
            except ProcessExecutionError as e:
                logging.info(e)

        # Complete masking using RepeatModeler/Masker
        def full(self):
            try:
                pass
            except ProcessExecutionError as e:
                logging.info(e)

    def __init__(self, input_paths: List[List[str]], cfg: ConfigManager, pm: PathManager,
                 record_ids: List[str], mode: int):
        # Get name of config section and required data
        name, ident, statement = Data().repeat_modeling()
        # Ensure protocol is passed and valid
        assert ConfigManager.PROTOCOL in cfg.config[name].keys()
        protocol = cfg.config.get(name, ConfigManager.PROTOCOL)
        workers = int(cfg.config.get(name, ConfigManager.WORKERS))
        super().__init__(
            # Call protocol function to generate list for superclass initialization
            [
                RepeatsIter.Repeats(
                    {Data.IN: input_path},
                    cfg,
                    pm,
                    record_id,
                    name,
                    mode,
                    required_data=[Data.IN],
                )
                for input_path, record_id in zip(input_paths, record_ids)
            ],
            # Remaining args as usual
            statement % (protocol + " search", workers, int(cfg.config.get(name, ConfigManager.THREADS))),
            workers,
            cfg,
            pm,
            mode
        )

    def run(self):
        super().run()

    def output(self) -> Tuple[List[List[str]], ConfigManager, PathManager, List[str], int]:
        return super().output()


if __name__ == "__main__":
    pass
