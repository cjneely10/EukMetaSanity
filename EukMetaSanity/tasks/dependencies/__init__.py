import sys
import inspect
from EukMetaSanity.tasks.dependencies.mmseqs import *
from EukMetaSanity.tasks.dependencies.repeat_modeler import *
from EukMetaSanity.tasks.dependencies.repeat_masker import *
from EukMetaSanity.tasks.dependencies.augustus import *
from EukMetaSanity.tasks.dependencies.genemark import *
from EukMetaSanity.tasks.dependencies.metaeuk import *

# Populate dependencies for easy loading
dependencies = {
    dep.name: dep for name, dep in inspect.getmembers(sys.modules[__name__])
    if isinstance(dep, type) and "EukMetaSanity" in repr(dep)
}
