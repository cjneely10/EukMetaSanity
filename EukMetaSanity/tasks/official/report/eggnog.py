from EukMetaSanity import Task, TaskList, program_catch


class EggNogMapperIter(TaskList):
    """ Task runs eggnog-mapper on a set of gene calls

    Outputs: emapper
    Finalizes: emapper

    """
    name = "eggnog"
    requires = []
    depends = ["emapper"]

    class EggNogMapper(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "final": ["emapper.emapper"]
            }

        @program_catch
        def run(self):
            pass

    def __init__(self, *args, **kwargs):
        super().__init__(EggNogMapperIter.EggNogMapper, EggNogMapperIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
