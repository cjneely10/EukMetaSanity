"""
Module handles taking all config files in the official pipeline and generating
"""
import glob
import os
from pathlib import Path

from plumbum import local


def update_config_files(data_path: str, eukms_bin_dir: str):
    """ Get all config files matching pattern *-config.* and replace the "/path/to" strings with their actual paths

    """
    user_uid = local["whoami"]()
    for config_file in glob.glob(os.path.join(os.path.dirname(eukms_bin_dir), "*", "*", "*-config.*")):
        updated_config_dir = os.path.dirname(os.path.dirname(config_file))
        new_config_file = open(os.path.join(updated_config_dir, os.path.basename(config_file)), "w")
        original_file_ptr = open(config_file, "r")
        for line in original_file_ptr:
            line = line.replace("uid", user_uid).replace("/path/to", str(Path(data_path).resolve()))
            new_config_file.write(line)
        new_config_file.close()
