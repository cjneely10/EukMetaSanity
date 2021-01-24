"""
Module holds functionality for downloading prerequisite data in official EukMetaSanity pipeline
"""
import os
from pathlib import Path
from typing import Iterable, Set, Optional, Sequence
from abc import abstractmethod
from plumbum import local
from plumbum.machines.local import LocalCommand


class DataUtil:
    """
    Parent class to all mmseqs data utilities. Defines functor that will handle specific utility operation.

    databases: set of database names under operation
    """
    def __init__(self, databases: Iterable[str]):
        """

        :param databases: Databases under operation
        """
        self._databases: Set[str] = set(databases)

    @abstractmethod
    def __call__(self) -> Optional[object]:
        """
        Child class implementation should perform given operation on stored database(s)
        """
        pass

    @property
    def databases(self) -> Set[str]:
        """ Get reference to set of stored databases

        :return: Set of stored database path/names
        """
        return self._databases

    @staticmethod
    def run(cmd: LocalCommand) -> Optional[object]:
        """ Print and run local command

        :param cmd: Local plumbum command
        """
        print(cmd)
        return cmd()


class Merge(DataUtil):
    """
    Combine two mmseqs2 databases into new database given path `new_db_name`
    """
    def __init__(self, new_db_name: str, databases: Iterable[str]):
        """ Merge databases into new_db_name

        :param new_db_name: Path/name for new database to generate
        :param databases: List of databases to merge
        """
        super().__init__(databases)
        self._new_db_name = new_db_name

    def __call__(self, *args, **kwargs):
        """
        Merge databases
        """
        return self.run(local["mmseqs"]["mergedbs", (*self.databases), self._new_db_name])


class CreateTaxDB(DataUtil):
    """
    Create taxonomy database using id mapping file
    """
    def __init__(self, id_mapping_files: Sequence[Path], wdir: str, databases: Sequence[str]):
        """ Create taxonomy database for each database in databases using files in id_mapping_files

        :param id_mapping_files: List of files to use in building taxonomy databases
        :param databases: List of databases to merge
        """
        super().__init__(databases)
        assert len(id_mapping_files) == len(self._databases)
        self._id_mapping_files = id_mapping_files
        self.wdir = wdir

    def __call__(self) -> Optional[object]:
        """
        Generate taxonomy database
        """
        # Download tax info
        self.run(local["wget"][
                     "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz",
                     "-O", os.path.join(self.wdir, "taxdump.tar.gz")
                 ])
        self.run(local["tar"][
                     "xzvf", os.path.join(self.wdir, "taxdump.tar.gz"), "-C", self.wdir
                 ])
        # Generate tax db
        for database, file in zip(self.databases, self._id_mapping_files):
            self.run(local["mmseqs"][
                            "createtaxdb",
                            database,
                            "tmp",
                            "--ncbi-tax-dump", self.wdir,
                            "--tax-mapping-file", database
                        ])
        return
