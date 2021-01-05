import os
from EukMetaSanity import Task, TaskList, program_catch, set_complete


class BuildDatabaseIter(TaskList):
    name = "repmod.build_database"
    requires = []
    depends = []
    
    class BuildDatabase(Task):
        @set_complete
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "db": [os.path.join(self.wdir, self.record_id)]
            }
            
        @program_catch
        def run(self):
            if len(os.listdir(self.wdir)) == 0:
                self.single(
                    self.program[
                        "-name", os.path.join(self.wdir, self.record_id),
                        str(self.dependency_input["fasta"]),
                    ]
                )
            
    def __init__(self, *args, **kwargs):
        super().__init__(BuildDatabaseIter.BuildDatabase, BuildDatabaseIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
