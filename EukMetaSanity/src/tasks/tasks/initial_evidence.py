from EukMetaSanity import Data, Task, TaskList, program_catch


class EvidenceIter(TaskList):
    class Evidence(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {Data.Type.OUT: []}

        def run(self) -> None:
            super().run()

        @program_catch
        def run_1(self):
            pass

    def __init__(self, *args, **kwargs):
        super().__init__(EvidenceIter.Evidence, "evidence", *args, **kwargs)
