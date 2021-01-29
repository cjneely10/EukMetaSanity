"""
Module holds functionality for downloading various data types and converting them to MMseqs format
"""
from plumbum import local
from EukMetaSanity.api.data.data import Data


class Fasta(Data):
    """
    Class represents a FASTA file data type download
    """

    def __init__(self, *args, **kwargs):
        """ FASTA format will be input into mmseqs createdb

        :param fasta_file: Path to FASTA file downloaded
        :param args: Args to pass to superclass
        :param kwargs: kwargs to pass to superclass
        """
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs) -> str:
        """
        Generate mmseqs database from a FASTA file
        """
        self.run(local["mmseqs"]["createdb", self.data, self.db_name])
        super().__call__(*args, **kwargs)


class MMSeqsDB(Data):
    """
    Class represents an MMSEQs database download
    """

    def __init__(self, *args, **kwargs):
        """ Database will simply be extracted, default functor used

        :param args: Args to pass to superclass
        :param kwargs: kwargs to pass to superclass
        """
        super().__init__(*args, **kwargs)


class MSA(Data):
    """
    Class represents a MSA to convert to MMseqs profile format
    """
    def __init__(self, *args, **kwargs):
        """ Format is FASTA and will convert to MMseqs profile

        :param args: Args to pass to superclass
        :param kwargs: kwargs to pass to superclass
        """
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs) -> str:
        """
        Create profile format
        """
        self.run(local["mmseqs"]["convertmsa", self.data, self.db_name + "-msa"])
        self.run(local["mmseqs"]["msa2profile", self.db_name + "-msa", self.db_name])
        super().__call__(*args, **kwargs)
