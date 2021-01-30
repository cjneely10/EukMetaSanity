"""
Module holds logic to generate MMseqs database indices
"""
import os
from typing import Sequence
from abc import abstractmethod
from plumbum import local
from EukMetaSanity.api.data.data import DataUtil


class Index(DataUtil):
    """
    Index class should be extended to determine exact mmseqs command to use at runtime
    """
    @abstractmethod
    def command(self) -> str:
        """ MMseqs index creation subprogram name

        :return: Subprogram name
        """
        pass

    def __init__(self, threads: int, wdir: str, databases: Sequence[str]):
        """ Create mmseqs lin/index base class

        :param threads: Number of system threads
        :param wdir: Directory containing built databases
        :param split_mem_limit: Maximum allowed memory
        :param databases: List of databases for which to generate indices
        """
        super().__init__(databases)
        self.wdir = wdir
        self._threads = threads

    def __call__(self):
        """
        Run subprogram to create lin/index
        """
        for database in self.databases:
            self.run(local["mmseqs"][
                self.command(),
                os.path.join(self.wdir, database),
                "tmp",
                "--threads", self._threads,
                "--remove-tmp-files"
            ])


class CreateIndex(Index):
    """
    Child class object will create mmseqs index
    """
    def command(self) -> str:
        """ MMseqs index subcommand

        :return: MMseqs index subcommand
        """
        return "createindex"


class CreateLinIndex(Index):
    """
    Child class object will create mmseqs linear index
    """
    def command(self) -> str:
        """ MMseqs linear index subcommand

        :return: MMseqs index subcommand
        """
        return "createlinindex"
