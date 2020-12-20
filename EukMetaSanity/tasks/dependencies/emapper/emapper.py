import os
from EukMetaSanity import Task, TaskList, program_catch


class EMapperIter(TaskList):
    name = "emapper"
    requires = []
    depends = []

    class EMapper(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "emapper": os.path.join(self.wdir, self.record_id + ".emapper"),
            }

        @program_catch
        def run(self):
            self.parallel(
                self.local[self.config["python"]][
                    self.config["program"],
                    "-i", self.input["root"]["prot"],
                    "--output", os.path.join(self.wdir, self.record_id),
                    "--cpu", self.threads,
                    (*self.added_flags),
                ]
            )
            os.replace(
                os.path.join(self.wdir, self.record_id + ".emapper.annotations"),
                str(self.output["emapper"])
            )

    def __init__(self, *args, **kwargs):
        super().__init__(EMapperIter.EMapper, EMapperIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
