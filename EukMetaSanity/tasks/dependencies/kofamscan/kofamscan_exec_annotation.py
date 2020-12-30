import os
from EukMetaSanity import Task, TaskList, program_catch, set_complete


class ExecAnnotationIter(TaskList):
    name = "kofamscan.exec_annotation"
    requires = []
    depends = []

    class ExecAnnotation(Task):
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "kegg": os.path.join(self.wdir, self.record_id + ".kegg"),
            }

        @program_catch
        def run(self):
            self.parallel(
                self.program[
                    "--cpu", self.threads,
                    "--format", "detail",
                    (*self.added_flags),
                    "-o", self.output["kegg"],
                    "--tmp-dir", os.path.join(self.wdir, "tmp"),
                    self.dependency_input["prot"],
                ]
            )

    def __init__(self, *args, **kwargs):
        super().__init__(ExecAnnotationIter.ExecAnnotation, ExecAnnotationIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
