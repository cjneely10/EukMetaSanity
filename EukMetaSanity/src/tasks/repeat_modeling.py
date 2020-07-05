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
        def __init__(self, input_path_dict: Dict[str, List[str]], cfg: ConfigManager, pm: PathManager, record_id: str,
                     mode: int):
            # Needs an input file only
            super().__init__(input_path_dict, cfg, pm, record_id, Data().repeat_modeling()[0], [Data.IN], mode)

        def run(self) -> None:
            super().run()

        def results(self) -> Dict[str, List[str]]:
            return super().results()

        def run_1(self):
            name = Data().repeat_modeling()[0]
            # Call method
            mmseqs_mask_db, masked_fasta = getattr(self, self.cfg.config.get(name, ConfigManager.PROTOCOL))()
            self.output = {Data.OUT: [
                self.input[Data.IN][0],  # Taxonomic results for ab initio prediction
                masked_fasta,  # Input FASTA file for ab initio
                mmseqs_mask_db,  # MMseqs database for use in metaeuk
            ]}

        # Simple repeat masking using mmseqs
        def simple(self):
            masked_db_path = os.path.join(self.wdir, self.record_id + "-mask_db")
            masked_fa_path = masked_db_path[:-3] + ".fna"
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
                        masked_fa_path,
                    ],
                    self.mode
                )
                # Return masked-db and FASTA file paths
                return masked_db_path, masked_fa_path
            except ProcessExecutionError as e:
                logging.info(e)

        # Complete masking using RepeatModeler/Masker
        def full(self):
            pass

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
                    mode
                )
                for input_path, record_id in zip(input_paths, record_ids)
            ],
            # Remaining args as usual
            statement % (protocol, workers, int(cfg.config.get(name, ConfigManager.THREADS))),
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
