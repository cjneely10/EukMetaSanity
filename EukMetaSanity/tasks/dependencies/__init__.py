import sys
import inspect
from EukMetaSanity.tasks.dependencies.mmseqs import *
from EukMetaSanity.tasks.dependencies.repeat_modeler import *
from EukMetaSanity.tasks.dependencies.repeat_masker import *

dependencies = {dep.name: dep for name, dep in inspect.getmembers(sys.modules[__name__])}
