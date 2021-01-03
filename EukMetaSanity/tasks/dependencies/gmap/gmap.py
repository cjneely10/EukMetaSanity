"""
Module holds gmap build functionality
"""

import os
from typing import List
from pathlib import Path
from EukMetaSanity import Task, TaskList, DependencyInput
from EukMetaSanity import MissingDataError
from EukMetaSanity import program_catch, prefix, set_complete


class GMapIter(TaskList):
    """ TaskList class iterates over gmap tasks

    name: gmap

    requires:

    depends: gmap.build

    output: sams

    final:

    """
    name = "gmap"
    requires = []
    depends = [DependencyInput("gmap.build")]

    class GMap(Task):
        """
        Task class handles gmap task
        """
        @set_complete
        def __init__(self, *args, **kwargs):
            """
            Instantiate class with given output
            """
            super().__init__(*args, **kwargs)
            if self.is_skip:
                sams = []
            else:
                sams = [os.path.join(self.wdir, prefix(transcript)) + ".sam" for transcript in self.get_transcripts()]
            self.output = {
                "sams": sams
            }

        @program_catch
        def run(self):
            """
            Run gmap
            """
            # Get transcripts
            for transcript, sam_file in zip(self.get_transcripts(), self.output["sams"]):
                # Generate genome index
                genome_idx = self.input["gmap.build"]["db"]
                _genome_dir = os.path.dirname(str(self.input["gmap.build"]["db"]))
                _genome_basename = os.path.basename(str(self.input["gmap.build"]["db"]))
                # Align
                self.parallel(
                    self.program[
                        "-D", _genome_dir, "-d", genome_idx,
                        "-t", self.threads,
                        transcript,
                        (*self.added_flags)
                    ] > str(sam_file)
                )

        def get_transcripts(self) -> List[str]:
            """ Parse transcriptome file into list of rna pairs to analyze

            :return: List of transcriptomes
            """
            if self.developer_mode:
                return []
            _path = str(Path(self.config["transcriptome"]).resolve())
            if not os.path.exists(_path):
                raise MissingDataError("Input transcriptome mapping file not found!")
            file_ptr = open(_path, "r")
            _id = self.record_id.replace("-mask", "")
            for line in file_ptr:
                if _id in line:
                    return line.rstrip("\r\n").split("\t")[1].split(",")
            return []

    def __init__(self, *args, **kwargs):
        """
        Instantiate TaskList
        """
        super().__init__(GMapIter.GMap, GMapIter.name, *args, **kwargs)
