from EukMetaSanity.src.tasks.taxonomy import TaxonomyIter
from EukMetaSanity.src.tasks.repeat_modeling import RepeatsIter


class TaskManager:
    def __init__(self):
        self.tasks = {
            "run": (
                TaxonomyIter,
                RepeatsIter,
            ),
            "refine": (),
            "report": (),
        }
