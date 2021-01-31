"""
Module holds functionality for manipulating downloaded databases from initial format into mmseqs databases
"""
import os
from typing import Sequence, Callable, List
from plumbum import local
from EukMetaSanity.api.data.data import DataUtil
from EukMetaSanity.tasks.utils.helpers import prefix


class ConcatDBs(DataUtil):
    """
    Combine two mmseqs2 databases into new database

    new_db_name: database name to create using all in list
    """
    def __init__(self, wdir: str, new_db_name: str, databases: Sequence[str]):
        """ Merge databases into new_db_name

        :param wdir: Working directory containing databases
        :param new_db_name: Path/name for new database to generate
        :param databases: List of databases to merge
        """
        super().__init__(databases)
        self._new_db_name = new_db_name
        self.wdir = wdir

    def __call__(self):
        """
        Merge databases
        """
        databases = [os.path.join(self.wdir, database) for database in self.databases]
        self.run(
            local["mmseqs"]["concatdbs", databases[0], (*databases[1:], os.path.join(self.wdir, self._new_db_name))]
        )
        self.run(
            local["mmseqs"]["concatdbs",
                            (*[database + "_h" for database in databases]),
                            os.path.join(self.wdir, self._new_db_name + "_h")]
        )


class CreateMappingFiles(DataUtil):
    """
    Class will create mapping files from the index file generated by MMseqs. To be used to generate
    taxonomic database consisting of NCBI taxids mapping to each FASTA sequence identifier
    """
    def __init__(self, parsing_functions: Sequence[Callable], wdir: str, databases: Sequence[str]):
        """ Create a mapping file using the identifiers from an MMseqs database .lookup file

        Parsing functions provided must take only two positional parameters:
            1) Path to database lookup file
            2) Path to which to output mapping file (existing path will be overwritten)

        :param parsing_functions: List of functions to call on databases to generate a mapping file for each
        :param wdir: Working directory in which to generate mapping files
        :param databases: List of databases for which to generate mapping files
        """
        assert len(parsing_functions) == len(databases)
        super().__init__(databases)
        self.wdir = wdir
        self._parsing_functions = parsing_functions

    def __call__(self):
        """
        Call parsing function on each database's lookup file
        """
        for parsing_function, database in zip(self._parsing_functions, self.databases):
            parsing_function(database + ".lookup", os.path.join(self.wdir, prefix(database) + ".input"))


class CreateTaxDBs(DataUtil):
    """
    Create taxonomy database using id mapping file

    _id_mapping_file: Path to file with proper id to ncbitaxonomy mapping (see mmseqs documentation)
    wdir: Working directory for download options
    """
    def __init__(self, wdir: str, databases: Sequence[str]):
        """ Create taxonomy database for each database in databases

        :param wdir: Working directory into which to download current ncbi taxonomy dump
        :param databases: List of databases for which to create taxonomy databases using existing mapping files
        """
        super().__init__(databases)
        self.wdir = wdir

    def __call__(self):
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
        for database in self.databases:
            self.run(local["mmseqs"][
                            "createtaxdb",
                            os.path.join(self.wdir, database),
                            "tmp",
                            "--ncbi-tax-dump", self.wdir,
                            (*self._get_tax_mapping_file(database))
                        ])

    def _get_tax_mapping_file(self, database: str) -> List[str]:
        """ Get path to generated taxonomy mapping file

        :param database: Database name
        :return: List to expand into plumbum command
        """
        tax_mapping = []
        tax_file = os.path.join(self.wdir, prefix(database) + ".input")
        if os.path.exists(tax_file):
            tax_mapping = ["--tax-mapping-file", os.path.join(self.wdir, prefix(database) + ".input")]
        return tax_mapping