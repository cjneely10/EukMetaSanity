from io import TextIOWrapper
from typing import Tuple, List, Dict
from EukMetaSanity.src.utils.path_manager import PathManager
from EukMetaSanity.src.tasks.task_class import TaskList, Task
from EukMetaSanity.src.utils.config_manager import ConfigManager


class AbInitioIter(TaskList):
    class AbInitio(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def run(self) -> None:
            # Provide required output
            super().run()

    def __init__(self, new_task: type, input_paths: List[List[str]], record_ids: List[str], data_function: Callable,
                 cfg: ConfigManager, pm: PathManager, mode: int):
        super().__init__(AbInitioIter.AbInitio, input_paths, record_ids, data_function, cfg, pm, mode)


def get_taxonomy(tax_results_file: TextIOWrapper) -> int:
    # Get first line
    tax_id: int = 2759  # Default to Eukaryota if nothing better is found
    try:
        while True:
            line = next(tax_results_file).rstrip("\r\n").split("\t")
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