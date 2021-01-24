"""
Module holds functionality for downloading prerequisite data in official EukMetaSanity pipeline
"""
from typing import Iterable, Set, Optional
from abc import abstractmethod
from plumbum import local


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
        Functor will call merge operation
        """
        return self.merge()

    def merge(self):
        """
        Merge databases
        """
        pass
