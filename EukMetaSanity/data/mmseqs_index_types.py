"""
Module holds logic to generate MMseqs database indices
"""
import os
from abc import abstractmethod
from typing import Sequence

from plumbum import local

from EukMetaSanity.data.data import DataUtil


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

    def __init__(self, threads: int, wdir: str, databases: Sequence[str], *args):
        """ Create mmseqs lin/index base class

        :param threads: Number of system threads
        :param wdir: Directory containing built databases
        :param databases: List of databases for which to generate indices
        :param args: List of additional arguments to pass to index creation operation
        """
        super().__init__(databases)
        self.wdir = wdir
        self._threads = threads
        self.args = args

    def __call__(self):
        """
        Run subprogram to create lin/index
        """
        for database in self.databases:
            self.run(local["mmseqs"][
                self.command(),
                os.path.join(self.wdir, database),
                os.path.join(self.wdir, "tmp"),
                "--threads", self._threads,
                "--remove-tmp-files",
                (*self.args)
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
