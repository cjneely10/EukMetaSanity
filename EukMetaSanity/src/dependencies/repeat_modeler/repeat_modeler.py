import glob
import os
import shutil
from typing import List, Union, Type

from plumbum import ProcessExecutionError

from yapim import Task, DependencyInput, touch


class RModRepeatModeler(Task):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output = {
            "model": os.path.join(self.wdir, "consensi.fa.classified")
        }

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        return []

    @staticmethod
    def depends() -> List[DependencyInput]:
        return [DependencyInput("RModBuildDatabase")]

    def run(self):
        """
        Run repmod.repeat_modeler
        """
        script = self.create_script(
            self.program[
                "-pa", str(int(self.threads) // 4 or 1),
                "-engine", "ncbi", "-database", str(self.input["RModBuildDatabase"]["db"][0])
            ],
            "modeler.sh"
        )
        try:
            self.parallel(script)
            _output = self.get_output_file()
            if len(_output) > 0:
                shutil.copyfile(_output[0], str(self.output["model"]))
            else:
                touch(str(self.output["model"]))
        except ProcessExecutionError:
            touch(str(self.output["model"]))

    def get_output_file(self):
        _output = glob.glob(os.path.join(self.wdir, "RM*", "consensi.fa.classified"))
        if len(_output) == 0:
            _output = glob.glob(os.path.join(self.wdir, "RM*", "consensi.fa"))
            if len(_output) == 0:
                _output = sorted(glob.glob(os.path.join(self.wdir, "RM*", "round*", "consensi.fa")), reverse=True)
                return _output
        return _output
