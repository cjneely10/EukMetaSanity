from typing import Dict, List, Type
from EukMetaSanity.tasks.base.task_class import TaskList
from EukMetaSanity.tasks.run.taxonomy import TaxonomyIter
from EukMetaSanity.tasks.run.ab_initio import AbInitioIter
from EukMetaSanity.tasks.run.repeat_modeling import RepeatsIter
from EukMetaSanity.tasks.run.initial_evidence import EvidenceIter


class TaskManager:
    def __init__(self):
        self._programs = {
            "run": [
                TaxonomyIter,
                RepeatsIter,
                AbInitioIter,
                EvidenceIter,
            ],
            "refine": [],
            "report": [],
        }

    @property
    def programs(self) -> Dict[str, List[Type[TaskList]]]:
        return self._programs
