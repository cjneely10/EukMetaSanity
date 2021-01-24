from plumbum import local
from EukMetaSanity.data.data import Data


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

    def __call__(self, *args, **kwargs):
        """
        Generate mmseqs database from a FASTA file
        """
        local["mmseqs"]["createdb", self.data, self.db_name]()


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

    def __call__(self, *args, **kwargs):
        """
        Overrides default call operator to create profile format
        """
        local["mmseqs"]["convertmsa", self.data, self.db_name + "-msa"]()
        local["mmseqs"]["msa2profile", self.db_name + "-msa", self.db_name]()
