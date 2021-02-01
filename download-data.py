#!/usr/bin/env python3
"""
Module downloads requisite data for official pipelines in EukMetaSanity
"""
import os
from pathlib import Path
from plumbum import cli
from EukMetaSanity.tasks.official.config_generation_utils import update_config_files
from EukMetaSanity.tasks.official.download_utils import download_data, parsing_operations, manage_downloaded_data


class DataDownloader(cli.Application):
    _working_dir: str = str(Path(os.path.join(os.getcwd(), "data")).resolve())
    _index: bool = False
    _threads: int = 1

    @cli.switch(["-x", "--index"], help="Generate search index (optional, but takes a lot of space)")
    def set_index(self, _index):
        self._index = _index

    @cli.switch(["-t", "--threads"], int, help="Number of threads to use in database generation, default 1")
    def set_threads(self, threads):
        if threads < 1:
            raise ValueError("Must pass positive number of threads")
        self._threads = threads

    def main(self):
        # Generate working directory
        if not os.path.exists(self._working_dir):
            os.makedirs(self._working_dir)

        # Download data
        print("Downloading data")
        for db_download in download_data(self._working_dir):
            db_download()

        # Parse any required download data
        print("\nGenerating taxonomy lookup files")
        for parsing_operation in parsing_operations(self._working_dir):
            parsing_operation()

        # Run database utility protocols
        print("\nRunning MMseqs utility functions on data")
        for util_instruction in manage_downloaded_data(self._working_dir,
                                                       self._index,
                                                       True,
                                                       self._threads):
            util_instruction()

        # Generate system-specific config files
        update_config_files(self._working_dir)


if __name__ == "__main__":
    update_config_files(os.path.join(os.getcwd(), "data"))
    # DataDownloader.run()
