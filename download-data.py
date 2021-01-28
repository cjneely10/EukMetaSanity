#!/usr/bin/env python3
"""
Module downloads requisite data for official pipelines in EukMetaSanity
"""
import os
from pathlib import Path
from plumbum import cli
from EukMetaSanity.api.data.download_utils import download_data, manage_downloaded_data


class DataDownloader(cli.Application):
    _working_dir: str = str(Path(os.path.join(os.getcwd(), "data")).resolve())
    _index: bool = False
    _threads: int = 1
    _max_mem: str = "8G"

    @cli.switch(["-x", "--index"], help="Generate search index (optional, but takes a lot of space)")
    def set_index(self, _index):
        self._index = _index

    @cli.switch(["-t", "--threads"], int, help="Number of threads to use in database generation, default 1")
    def set_threads(self, threads):
        if threads < 1:
            raise ValueError("Must pass positive number of threads")
        self._threads = threads

    @cli.switch(["-m", "--max-mem"], str, help="Set max memory per split. E.g. 800B, 5K, 10M, 1G; default 8G")
    def set_max_memory(self, max_memory):
        if max_memory[-1] not in ("B", "K", "M", "G", "T"):
            raise ValueError("Memory string must be specific format - see help menu")
        self._max_mem = max_memory

    def main(self):
        # Generate working directory
        if not os.path.exists(self._working_dir):
            os.makedirs(self._working_dir)

        # Download data
        downloaded_databases = []
        for db_download in download_data(self._working_dir):
            downloaded_databases.append(db_download())

        # TODO
        # Create any needed lookup files

        # Run database utility protocols
        for util_instruction in manage_downloaded_data(self._working_dir, self._index, True, self._threads,
                                                       self._max_mem):
            util_instruction()

        # TODO
        # Generate base config files


def _odb_tax_parse(mmseqs_db_path: str, outfile: str):
    """ Generate NCBI taxonomy mappings from generated .lookup file

    :param mmseqs_db_path: Path to generated mmseqs database
    :param outfile: Output path for parsed orthodb mapping file
    """
    output_p = open(outfile, "w")
    mmseqs_input_fp = open(mmseqs_db_path + ".lookup", "r")
    for line in mmseqs_input_fp:
        line = line.split()
        output_p.write(line[1] + "\t" + line[1].split("_")[0] + "\n")
    output_p.close()


if __name__ == "__main__":
    DataDownloader.run()
