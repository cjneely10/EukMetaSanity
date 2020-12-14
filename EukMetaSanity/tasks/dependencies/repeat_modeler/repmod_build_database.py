import os
from EukMetaSanity import Task, TaskList, program_catch


class BuildDatabaseIter(TaskList):
    name = "repmod.build_database"
    requires = []
    depends = []
    
    class BuildDatabase(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "db": os.path.join(self.wdir, self.record_id)
            }
            
        @program_catch
        def run(self):
            self.single(
                self.program[
                    "-name", self.output["db"],
                    str(self.input["root"]["fna"]),
                ]
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(BuildDatabaseIter.BuildDatabase, BuildDatabaseIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
