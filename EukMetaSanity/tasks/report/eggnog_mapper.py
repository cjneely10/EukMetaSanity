import os
from EukMetaSanity import TaskList, Task, program_catch


class EggNOGMapper(TaskList):
    class EggNOG(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            pass

    def __init__(self, *args, **kwargs):
        super().__init__(EggNOGMapper.EggNOG, "eggnog", *args, **kwargs)


if __name__ == "__main__":
    pass
