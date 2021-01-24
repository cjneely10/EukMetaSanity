"""
Module houses functionality to use in EukMS installation to download all required data dependencies for
official EukMS pipelines
"""
from typing import Iterator, Generator
from EukMetaSanity.data.data import Fasta, MMSeqsDB
from EukMetaSanity.data.data_utils import Merge, CreateTaxDB


def download_data(wdir: str) -> Generator:
    database_downloads = [
        Fasta(
            config_identifier="ortho_db",
            data_url="https://v101.orthodb.org/download/odb10v1_all_og_fasta.tab.gz",
            wdir=wdir,
            expected="odb10v1_all_og_fasta.tab",
            unzip_command_args=["gunzip"],
        ),
        MMSeqsDB(
            config_identifier="mmetsp_db",
            data_url="https://wwwuser.gwdg.de/~compbiol/metaeuk/2019_11/TAX_DBs/MMETSP/TaxDB_MMETSP.tar.gz",
            wdir=wdir,
            expected="MMETSP",
            unzip_command_args=["tar", "-xzf", "-C", wdir],
        ),
    ]
    for db in database_downloads:
        yield db


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
