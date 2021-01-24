"""
Module houses functionality to use in EukMS installation to download all required data dependencies for
official EukMS pipelines
"""
from typing import Iterator
from EukMetaSanity.data.data_utils import Merge, CreateTaxDB


def instructions(working_dir: str) -> Iterator:
    """ Generate data utilities functors. Consumer should call each object in sequence.

    Current implementation:

    Download/generate mmseqs db for ODB and generate taxonomy database

    Download/generate mmseqs db for RFAM

    Download/generate mmseqs db for MMETSP

    Merge ODB and MMETSP databases into single database

    :return: Iterator over all official EukMS pipeline downloader classes
    """
    # TODO: generate input file for ODB parsing
    fxns = [
        CreateTaxDB("mmseqs.input", working_dir, ["ortho_db"]),
        Merge("odb-mmetsp_db", ["ortho_db", "mmetsp_db"]),
    ]
    for fxn in fxns:
        yield fxn
