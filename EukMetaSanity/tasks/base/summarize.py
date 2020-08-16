import os
from EukMetaSanity import Task, TaskList

"""
Class automatically converts output to summary format
If a dict of data is present, its data is populated to be the final output
Otherwise, the data in the output is labelled using the extension on the file

"""


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

    def __init__(self, *args, **kwargs):
        super().__init__(SummarizeIter.Summarize, "summarize", *args, **kwargs)


if __name__ == "__main__":
    pass
