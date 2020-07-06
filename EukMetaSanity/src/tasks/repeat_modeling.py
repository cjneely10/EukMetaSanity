import os
import logging
from typing import List
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
            self.output = {Data.Type.OUT: [
                self.input[Data.Type.IN][0],  # Taxonomic results for ab initio prediction
                masked_fa_path,  # Input FASTA file for ab initio
                masked_db_path,  # MMseqs database for use in metaeuk
            ]}
            super().run()

        def run_1(self):
            name = Data().repeats()[0]
            # Call protocol method
            getattr(self, self.cfg.config.get(name, ConfigManager.PROTOCOL))()

        # Simple repeat masking using mmseqs
        def simple(self):
            masked_db_path = self.output[Data.Type.OUT][2]
            try:
                # Generate the masked sequence
                self.log_and_run(
                    self.program[
                        "masksequence",
                        self.input[Data.Type.IN][2],
                        masked_db_path,
                        "--threads", self.threads,
                    ]
                )
                # Output as FASTA file
                self.log_and_run(
                    self.program[
                        "convert2fasta",
                        masked_db_path,
                        self.output[Data.Type.OUT][1],
                    ]
                )
            except ProcessExecutionError as e:
                logging.info(e)

        # Complete masking using RepeatModeler/Masker
        def full(self):
            try:
                # Build database
                self.log_and_run(
                    self.program[
                        "-name", self.output[Data.Type.OUT][2],
                        (*self.added_flags),
                        self.input[Data.Type.IN][1],
                    ]
                )
                # Run RepeatModeler
                self.log_and_run(
                    self.program2[
                        "-pa", self.threads,
                        (*self.added_flags),
                        "-database", self.input[Data.Type.ACCESS][0],
                    ]
                )
                # Rename results
                if self.mode == 1:
                    os.replace(
                        os.path.join(self.wdir, self.record_id + "-families.fa"),
                        self.output[Data.Type.OUT][1],
                    )
            except ProcessExecutionError as e:
                logging.info(e)

    def __init__(self, input_paths: List[List[str]], cfg: ConfigManager, pm: PathManager,
                 record_ids: List[str], mode: int):
        dt = Data()
        protocol = cfg.config.get(dt.repeats()[0], ConfigManager.PROTOCOL)
        if protocol == "simple":
            super().__init__(RepeatsIter.Repeats, input_paths, record_ids, dt.repeats, cfg, pm, mode)
        else:
            super().__init__(RepeatsIter.Repeats, input_paths, record_ids, dt.repeats, cfg, pm, mode,
                             {Data.Type.ACCESS: [dt.repeats()[1]]})


if __name__ == "__main__":
    pass
