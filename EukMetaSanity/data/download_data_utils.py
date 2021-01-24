from typing import Iterator
from EukMetaSanity.data.data_utils import Merge, CreateTaxDB


def instructions(working_dir: str) -> Iterator:
    """ Run all data utilities needed for official EukMS pipelines

    :return: Iterator over all official EukMS pipeline database downloads
    """
    fxns = [
        Merge("odb-mmetsp_db", ["ortho_db", "mmetsp_db"]),
        CreateTaxDB("mmseqs.input", working_dir, ["ortho_db"])
    ]
    for fxn in fxns:
        yield fxn
