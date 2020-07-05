from EukMetaSanity.src.tasks.repeat_modeling import RepeatsIter
from EukMetaSanity.src.tasks.taxonomy import TaxonomyIter


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
