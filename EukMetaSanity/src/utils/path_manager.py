import os
from pathlib import Path
from plumbum import local
from numba import types, jit
from numba.experimental import jitclass

mkdir = local["mkdir"]


@jitclass((
        ("_base", types.string),
        ("_wdir", types.string),
))
class PathManager:

    @jit(types.string(types.string), nopython=True)
    def __init__(self, base_path):
        self._wdir = "wdir"
        self._base = str(Path(base_path).resolve())
        self._generate_directory_tree()

    # Get working dir
    @property
    def wdir(self):
        return os.path.join(self._base, self._wdir)

    # Get base dir
    @property
    def base(self):
        return str(self._base)

    # Add record directory to wdir
    @jit(types.void(types.string, types.List), nopython=True, cache=True)
    def add_dir(self, record_id, _subdirs=None):
        # Record base dir
        mkdir["-p", os.path.join(self.wdir, str(record_id))]()
        # Additional dirs, if needed
        if _subdirs is not None:
            all(mkdir["-p", os.path.join(self.wdir, str(record_id), _subd)]() for _subd in _subdirs)

    # Get record directory in wdir
    @jit(types.string(types.string), nopython=True, cache=True)
    def get_record_dir(self, record_id=None):
        loc = self.wdir
        if record_id is not None:
            loc = os.path.join(self.wdir, str(record_id))
        if os.path.exists(loc):
            return loc
        raise ValueError("%s dir does not exist in %s\n" % (record_id, self.wdir))

    # Generate initial directory tree
    @jit(types.void(), nopython=True, cache=True)
    def _generate_directory_tree(self):
        # Base output directory
        mkdir["-p", self._base]()
        # Working directory
        mkdir["-p", os.path.join(self._base, self._wdir)]()


if __name__ == "__main__":
    pass
