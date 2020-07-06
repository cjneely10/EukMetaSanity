import os
from typing import List, Dict
from EukMetaSanity.src.utils.data import Data
from EukMetaSanity.bin.fastagff3_to_gb import write_genbank
from EukMetaSanity.src.utils.path_manager import PathManager
from EukMetaSanity.src.utils.config_manager import ConfigManager
from EukMetaSanity.src.tasks.task_class import TaskList, Task, program_catch


class AbInitioIter(TaskList):
    class AbInitio(Task):
        def __init__(self, input_path_dict: Dict[Data.Type, List[str]], cfg: ConfigManager, pm: PathManager,
                     record_id: str, db_name: str, mode: int, required_data: List[str]):
            super().__init__(input_path_dict, cfg, pm, record_id, db_name, mode, required_data)

        def run(self) -> None:
            # Only looking for final trained ab initio prediction
            self.output = {Data.Type.OUT: [
                os.path.join(self.wdir, self.record_id + ".gff3"),  # Output gff3 ab initio predictions, final round
            ]}
            super().run()

        def run_1(self):
            name = Data().abinitio()[0]
            # Call protocol method
            getattr(
                self, self.cfg.config.get(name, ConfigManager.PROTOCOL)
            )(AbInitioIter.get_taxonomy(self.input[Data.Type.IN][0]))

        def augustus(self, tax_id: int):
            self._augustus(str(tax_id), 1)
            name = Data().abinitio()[0]
            for i in range(int(self.cfg.config.get(name, ConfigManager.ROUNDS)) - 1):
                self._augustus(self.record_id + str(i + 2), i + 2)

        @program_catch
        def _augustus(self, species: str, pos: int):
            out_gff = AbInitioIter.AbInitio._out_path(self.input[Data.Type.IN][1], ".%i.gff3" % pos)
            # Run prediction
            self.log_and_run(
                self.program[
                    "--codingseq=on",
                    "--stopCodonExcludedFromCDS=true",
                    "--species=%s" % species,
                    self.input[Data.Type.IN][1],
                    "--outfile", out_gff
                ]
            )
            # Parse to genbank
            out_gb = AbInitioIter.AbInitio._out_path(self.input[Data.Type.IN][1], ".%i.gb" % pos)
            write_genbank(
                self.input[Data.Type.IN][1],
                out_gff,
                out_gb
            )
            species_config_prefix = self.record_id + str(pos)
            # Write new species config file
            self.log_and_run(
                self.program2[
                    "--species=%s" % species_config_prefix,
                    out_gb
                ]
            )
            # Run training
            self.log_and_run(
                self.program3[
                    "--species=%s" % species_config_prefix,
                    out_gb
                ]
            )

        @program_catch
        def gmes(self, tax_id: int):
            pass

        @staticmethod
        def _out_path(_file_name: str, _ext: str) -> str:
            _file_name = _file_name.split(".")
            return ".".join(_file_name[:-1]) + _ext

    def __init__(self, input_paths: List[List[str]], cfg: ConfigManager, pm: PathManager, record_ids: List[str],
                 mode: int):
        dt = Data()
        super().__init__(AbInitioIter.AbInitio, input_paths, record_ids, dt.abinitio, cfg, pm, mode)

    @staticmethod
    def get_taxonomy(tax_results_file: str) -> int:
        _tax_results_file = open(tax_results_file, "r")
        # Get first line
        tax_id: int = 2759  # Default to Eukaryota if nothing better is found
        try:
            while True:
                line = next(_tax_results_file).rstrip("\r\n").split("\t")
                # Parse line for assignment
                _score, _tax_id, _assignment = float(line[0]), line[4], line[5].replace(" ", "")
                if _assignment in ("unclassified", "root"):
                    continue
                if _score < 80.0:
                    break
                # Keep new value if >= 80.0% of contigs map to the taxonomy
                else:
                    tax_id = _tax_id
        except StopIteration:
            return tax_id
