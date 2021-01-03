import sys
import inspect
from EukMetaSanity.tasks.dependencies.mmseqs import *
from EukMetaSanity.tasks.dependencies.repeat_modeler import *
from EukMetaSanity.tasks.dependencies.repeat_masker import *
from EukMetaSanity.tasks.dependencies.augustus import *
from EukMetaSanity.tasks.dependencies.genemark import *
from EukMetaSanity.tasks.dependencies.metaeuk import *
from EukMetaSanity.tasks.dependencies.emapper import *
from EukMetaSanity.tasks.dependencies.kofamscan import *
from EukMetaSanity.tasks.dependencies.hisat2 import *
from EukMetaSanity.tasks.dependencies.gmap import *
from EukMetaSanity.tasks.dependencies.sambamba import *
from EukMetaSanity.tasks.dependencies.hisat2.hisat2 import Hisat2Iter
from EukMetaSanity.tasks.dependencies.hisat2.hisat2_build import Hisat2BuildIter
from EukMetaSanity.tasks.dependencies.gmap.gmap import GMapIter
from EukMetaSanity.tasks.dependencies.gmap.gmap_build import GMapBuildIter
from EukMetaSanity.tasks.dependencies.sambamba.sambamba_view import SambambaViewIter
from EukMetaSanity.tasks.dependencies.sambamba.sambamba_sort import SambambaSortIter

# Populate dependencies for easy loading
# TODO: Add `expects` to docstrings for all dependencies
dependencies = {
    dep.name: dep for name, dep in inspect.getmembers(sys.modules[__name__])
    if isinstance(dep, type) and "EukMetaSanity" in repr(dep)
}
