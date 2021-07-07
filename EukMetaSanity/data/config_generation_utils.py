"""
Module handles taking all config files in the official pipeline and generating
"""
import os
import glob
from pathlib import Path


def update_config_files(data_path: str):
    """ Get all config files matching pattern *-config.* and replace the "/path/to" strings with their actual paths

    """
    official_pipelines_dir = os.path.dirname(__file__)
    for config_file in glob.glob(os.path.join(official_pipelines_dir, "*", "*-config.*")):
        new_config_file = open(os.path.join(data_path, os.path.basename(config_file)), "w")
        original_file_ptr = open(config_file, "r")
        for line in original_file_ptr:
            line = line.replace("/path/to/data", str(Path(data_path).resolve()))
            new_config_file.write(line)
        new_config_file.close()
