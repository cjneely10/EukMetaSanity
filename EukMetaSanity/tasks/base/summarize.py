import os
from EukMetaSanity import Task, TaskList


class SummarizeIter(TaskList):
    class Summarize(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = self.input
            has_dict = False
            for _obj in self.output:
                if isinstance(_obj, dict):
                    has_dict = True
                    break
            if not has_dict:
                _data = {}
                for _val in self.output:
                    if isinstance(_val, str):
                        _prefix, _ext = os.path.splitext(_val)
                        _data.update({_ext[1:]: _val})
                self.output.append(_data)

        def run(self):
            super().run()

    def __init__(self, *args, **kwargs):
        super().__init__(SummarizeIter.Summarize, "summarize", *args, **kwargs)


if __name__ == "__main__":
    pass
