"""
Module holds logic to generate MMseqs database indices
"""
# pylint: disable=useless-super-delegation
from typing import Sequence
from abc import abstractmethod
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

    def __init__(self, threads: int, split_mem_limit: str, generate: bool, databases: Sequence[str]):
        """ Create mmseqs lin/index base class

        :param threads: Number of system threads
        :param split_mem_limit: Maximum allowed memory
        :param generate: Flag determining if this index should be generated (or skip index generation)
        :param databases: List of databases for which to generate indices
        """
        super().__init__(databases)
        self._threads = threads
        self._mem_limit = split_mem_limit
        self._generate = generate

    def __call__(self):
        """
        Run subprogram to create lin/index
        """
        if not self._generate:
            return
        for database in self.databases:
            self.run(local["mmseqs"][
                self.command(),
                database,
                "tmp",
                "--threads", self._threads,
                "--split-memory-limit", self._mem_limit,
                "--remove-tmp-files"
            ])


class CreateIndex(Index):
    """
    Child class object will create mmseqs index
    """
    def __init__(self, *args, **kwargs):
        """
        Create index child class
        """
        super().__init__(*args, **kwargs)

    def command(self) -> str:
        """ MMseqs index subcommand

        :return: MMseqs index subcommand
        """
        return "createindex"


class CreateLinIndex(Index):
    """
    Child class object will create mmseqs linear index
    """
    def __init__(self, *args, **kwargs):
        """
        Create linindex child class
        """
        super().__init__(*args, **kwargs)

    def command(self) -> str:
        """ MMseqs linear index subcommand

        :return: MMseqs index subcommand
        """
        return "createlinindex"
