import os
import pickle
from collections import defaultdict
from typing import List, Dict, Tuple, Optional
from Bio import SeqIO
from EukMetaSanity import prefix
from EukMetaSanity.tasks.base.path_manager import PathManager
from EukMetaSanity.tasks.base.config_manager import ConfigManager


# Get program-needed list of files for this step in pipeline
def _get_list_of_files(summary_file: str) -> Tuple[List[str], List[Dict[str, Dict[str, object]]]]:
    data = pickle.load(open(summary_file, "rb"))
    out_ids = sorted(list(data.keys()))
    out_dict_list = []
    for _id in out_ids:
        to_add = {"root": {}}
        for key, val in data[_id].items():
            if isinstance(val, dict):
                to_add["root"][key] = val
            else:
                to_add[key] = val
        out_dict_list.append(to_add)
    return out_ids, out_dict_list


# TODO: Handle loading from wdir all required input for specific pipeline
class InputManager:
    def __init__(self, output_dir: str, input_dir: Optional[str], pm: PathManager, cfg: ConfigManager,
                 extension_list: List[str]):
        """ Generate object using existing pipeline directory structure and

        :param output_dir:
        :param input_dir:
        :param extension_list:
        """
        self.pm = pm
        self.pm.add_dirs(self.pm.MAGS)
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.extension_list = extension_list
        self.data: Dict[str, Dict[str, object]] = {}
        self.get_input_files()
        self.base = cfg.config[ConfigManager.INPUT][ConfigManager.BASE]
        self.record_ids = sorted(list(self.data.keys()))

    @staticmethod
    def _simplify_fasta(fasta_file: str, record_id: str, storage_dir: str, extension: str) -> str:
        """ Simplify header in FASTA file using BioPython

        :param fasta_file: Path to FASTA file
        :param storage_dir: Directory to write simplified file
        :param extension: Extension to give file
        :return: Path to simplified FASTA file
        """
        # Get path to file
        out_file = os.path.join(storage_dir, record_id + extension)
        # Do not overwrite if file already present
        if os.path.exists(out_file):
            return out_file
        # Write simplified version
        SeqIO.write(SeqIO.parse(fasta_file, "fasta"), out_file, "fasta")
        return out_file

    def get_input_files(self):
        _input_files = []
        # TODO: This is still buggy
        for file in os.listdir(self.input_dir):
            for ext in self.extension_list:
                if file.endswith(ext):
                    record_id = prefix(file)
                    if record_id not in self.data.keys():
                        self.data[record_id] = {}
                    self.data[record_id][ext] = InputManager._simplify_fasta(file, record_id, self.input_dir, ext)
                    self.pm.add_dirs(record_id)

    def _parse_results_dir_for_needed_data(self):
        pass

    @property
    def input_files(self):
        return [{self.base: self.data[record_id]} for record_id in self.record_ids]

    @property
    def input_prefixes(self):
        return self.record_ids
