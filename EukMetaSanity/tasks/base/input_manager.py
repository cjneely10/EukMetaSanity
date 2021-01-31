"""
Module holds InputManager class
"""

import os
import sys
import pickle
from pathlib import Path
from typing import List, Dict, Optional
from Bio import SeqIO
# pylint: disable=no-member
from plumbum import colors
from EukMetaSanity import prefix
from EukMetaSanity.tasks.base.path_manager import PathManager
from EukMetaSanity.tasks.base.config_manager import ConfigManager


class InputManager:
    """ InputManager converts input passed at the user-level with the existing project structure
    and requested input at config level into single input dict for each task

    """

    def __init__(self, output_dir: str, input_dir: Optional[str], pm: PathManager, cfg: ConfigManager,
                 extension_list: List[str]):
        """ Generate object using existing pipeline directory structure and

        :param output_dir:
        :param input_dir:
        :param extension_list:
        """
        self.path_manager = pm
        self.path_manager.add_dirs(self.path_manager.MAGS)
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.extension_list = extension_list
        self.data: Dict[str, Dict[str, object]] = {}
        self.bases = cfg.config[ConfigManager.INPUT][ConfigManager.BASE].split(" ")
        for base in self.bases:
            if base == ConfigManager.ROOT:
                # Load input files from command line level
                self._get_input_files()
            else:
                # Load additional input based on existing data, if requested
                self._get_files_from_project_dir(os.path.join(
                    output_dir, ConfigManager.EXPECTED_RESULTS_DIR, base, base + ".pkl"
                ), cfg)
        self.record_ids = sorted(list(self.data.keys()))
        if len(self.record_ids) == 0:
            print(colors.bold & colors.warn | "No input files were found, exiting")
            sys.exit()

    @staticmethod
    def _simplify_fasta(fasta_file: str, record_id: str, storage_dir: str, extension: str) -> str:
        """ Simplify header in FASTA file using BioPython

        :param fasta_file: Path to FASTA file
        :param storage_dir: Directory to write simplified file
        :param extension: Extension to give file
        :return: Path to simplified FASTA file
        """
        # Get path to file
        out_file = str(Path(os.path.join(storage_dir, record_id + extension)).resolve())
        # Do not overwrite if file already present
        if os.path.exists(out_file):
            return out_file
        # Write simplified version
        records = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            record.description = ""
            record.name = ""
            records.append(record)
        SeqIO.write(records, out_file, "fasta")
        return out_file

    def _get_files_from_project_dir(self, summary_file: str, cfg: ConfigManager):
        """ Parse existing project directory, if it exists, for input data requested at the config
        file level

        :param summary_file: Path to summary file, if it exists
        :return:
        """
        if not os.path.exists(summary_file):
            print(colors.bold & colors.warn |
                  "Config file requested input files, but metadata file '%s' is not present" % summary_file)
            sys.exit()
        data: Dict[str, Dict[str, object]] = pickle.load(open(summary_file, "rb"))
        for record_id in data.keys():
            if record_id not in self.data.keys():
                self.data[record_id] = {}
            self.data[record_id].update(data[record_id])
            if "input" in cfg.config[ConfigManager.INPUT].keys():
                input_adjust_dict = cfg.config[ConfigManager.INPUT]["input"]
                self.data[record_id][list(input_adjust_dict.keys())[0]] = self.data[record_id][
                    list(input_adjust_dict.values())[0]]

    def _get_input_files(self):
        """ Load input files from passed input directory at command-line level

        :raises: AttributeError if identical record_ids are found
        """
        if self.input_dir is None:
            print(colors.bold & colors.warn |
                  "Config file requested input files at the command-line, but none were found")
            sys.exit()
        for file in os.listdir(self.input_dir):
            for ext in self.extension_list:
                if file.endswith(ext):
                    record_id = prefix(file)
                    if record_id in self.data.keys():
                        raise AttributeError("Input directory has multiple files with the same basename - exiting")
                    self.data[record_id] = {}
                    self.data[record_id]["fasta"] = InputManager._simplify_fasta(os.path.join(self.input_dir, file),
                                                                                 record_id,
                                                                                 self.path_manager.get_dir(
                                                                                     self.path_manager.MAGS),
                                                                                 ext)
                    self.path_manager.add_dirs(record_id)

    @property
    def input_prefixes(self) -> List[str]:
        """ List of record_ids identified after parsing input directory and existing results directory

        :return: List of record_ids
        """
        return self.record_ids

    @property
    def input_files(self) -> List[Dict[str, Dict[str, object]]]:
        """ List of input files (parsed into dictionary)

        :return: List of input files in order based on self.record_ids
        """
        return [{ConfigManager.ROOT: self.data[record_id]} for record_id in self.record_ids]

    def summarize_input(self) -> str:
        """ Method provides a summary of all input files that are provided to the run.

        Meant to be called at outer level as a means of displaying for user the files that will/may be
        involved in analysis

        :return: Str of input
        """
        input_keys = {key for data_dict in self.input_files for key in data_dict.keys() if key != ConfigManager.ROOT}
        return "\n".join((
            colors.bold & colors.blue | "\n".join((
                f"{record_id}: {len(self.input_files[i][ConfigManager.ROOT])} file(s) from command line, "
                f"{sum([len(self.input_files[i][key]) for key in input_keys])} file(s) from internal pipelines"
                for i, record_id in enumerate(self.input_prefixes)
            )),
            colors.reset | "",
            (colors.reset | "") +
            (colors.bold & colors.blue | f"\nNumber of unique records gathered: {len(self.input_prefixes)}"),
        ))
