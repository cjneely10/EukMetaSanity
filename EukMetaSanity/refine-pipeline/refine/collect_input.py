from pathlib import Path
from typing import List, Union, Type, Tuple, Dict

from yapim import AggregateTask, DependencyInput


class CollectInput(AggregateTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def deaggregate(self) -> dict:
        input_data = {}
        for record_id, rnaseq_pair_list in CollectInput.get_rna_read_pairs(Path(self.config["rnaseq"])).items():
            if record_id not in input_data.keys():
                input_data[record_id] = {"rna_read_pairs": [], "transcripts": []}
            input_data[record_id]["rna_read_pairs"].extend(rnaseq_pair_list)
        for record_id, transcriptome_list in CollectInput.get_transcripts(Path(self.config["transcriptomes"])).items():
            if record_id not in input_data.keys():
                input_data[record_id] = {"rna_read_pairs": [], "transcripts": []}
            input_data[record_id]["transcripts"].extend(transcriptome_list)
        return input_data

    @staticmethod
    def requires() -> List[Union[str, Type]]:
        pass

    @staticmethod
    def depends() -> List[DependencyInput]:
        pass

    def run(self):
        pass

    @staticmethod
    def get_rna_read_pairs(rnaseq_mapping_file: Path) -> Dict[str, List[Tuple]]:
        if not rnaseq_mapping_file.exists():
            return {}
        file_ptr = open(rnaseq_mapping_file, "r")
        out = {}
        for line in file_ptr:
            line = line.rstrip("\r\n").split()
            pairs_string = line[1].split(";")
            if line[0] not in out.keys():
                out[line[0]] = []
            for pair in pairs_string:
                out[line[0]].append(tuple(pair.split(",")))
        file_ptr.close()
        return out

    @staticmethod
    def get_transcripts(transcripts_mapping_file: Path) -> Dict[str, List[str]]:
        if not transcripts_mapping_file.exists():
            return {}
        file_ptr = open(transcripts_mapping_file, "r")
        out = {}
        for line in file_ptr:
            line = line.rstrip("\r\n").split("\t")
            if line[0] not in out.keys():
                out[line[0]] = []
            out[line[0]].extend(line[1].split(","))
        return out
