"""
Module houses functionality to use in EukMS installation to download all required data dependencies for
official EukMS pipelines
"""
import os
from pathlib import Path
from typing import Iterator, Generator
from EukMetaSanity.data.data_types import Fasta, MMSeqsDB
from EukMetaSanity.data.mmseqs_types import Merge, CreateTaxDB


def download_data(working_dir: str) -> Generator:
    """ Generate data download functors. Consumer should call each object in sequence.

    Current implementation:

    Download/generate mmseqs db for ODB

    Download/generate mmseqs db for RFAM

    Download/generate mmseqs db for MMETSP

    :param working_dir: Root directory for data downloads
    :return: Iterator over all official EukMS pipeline downloader functors
    """
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    database_downloads = (
        Fasta(
            config_identifier="ortho_db",
            data_url="https://v101.orthodb.org/download/odb10v1_all_og_fasta.tab.gz",
            wdir=working_dir,
            expected="odb10v1_all_og_fasta.tab",
            unzip_command_args=["gunzip"],
        ),
        MMSeqsDB(
            config_identifier="mmetsp_db",
            data_url="https://wwwuser.gwdg.de/~compbiol/metaeuk/2019_11/TAX_DBs/MMETSP/TaxDB_MMETSP.tar.gz",
            wdir=working_dir,
            expected="MMETSP",
            unzip_command_args=["tar", "-xzf", "-C", working_dir],
        ),
    )
    for db in database_downloads:
        yield db


def instructions(working_dir: str) -> Iterator:
    """ Generate data utilities functors. Consumer should call each object in sequence.

    Current implementation:

    Generate ODB taxonomy database

    Merge ODB and MMETSP databases into single database

    :return: Iterator over all official EukMS download utility functors
    """
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    fxns = (
        CreateTaxDB("mmseqs.input", working_dir, ["ortho_db"]),
        Merge("odb-mmetsp_db", ["ortho_db", "mmetsp_db"]),
    )
    for fxn in fxns:
        yield fxn
