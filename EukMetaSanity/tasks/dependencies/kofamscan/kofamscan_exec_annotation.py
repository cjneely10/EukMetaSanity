"""
Module holds kofamscan.exec_annotation build functionality
"""
import os
from EukMetaSanity import Task, TaskList, program_catch, set_complete


class ExecAnnotationIter(TaskList):
    """ TaskList class iterates over kofamscan.exec_annotation tasks

    name: kofamscan.exec_annotation

    requires:

    depends:

    expects: prot[Path]

    output: kegg[Path]

    config:
        kofamscan.exec_annotation:
          program: /path/to/kofamscan/exec_annotation
          kolist: /path/to/kofam/ko_list
          profiles: /path/to/profiles/eukaryote.hal

    """
    name = "kofamscan.exec_annotation"
    requires = []
    depends = []

    class ExecAnnotation(Task):
        """
        Task class handles kofamscan.exec_annotation task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            self.output = {
                "kegg": os.path.join(self.wdir, self.record_id + ".kegg"),
            }

        @program_catch
        def run(self):
            """
            Run kofamscan.exec_annotation
            """
            self.parallel(
                self.program[
                    "--cpu", self.threads,
                    "--format", "detail",
                    "-o", self.output["kegg"],
                    "--tmp-dir", os.path.join(self.wdir, "tmp"),
                    self.dependency_input["prot"],
                    "-p", self.config["profiles"],
                    "-k", self.config["kolist"]
                ]
            )

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(ExecAnnotationIter.ExecAnnotation, ExecAnnotationIter.name, *args, **kwargs)
