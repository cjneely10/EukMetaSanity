import os
from collections import defaultdict
from typing import DefaultDict, List, Tuple
from EukMetaSanity.tasks.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch

"""
De novo assembly RNAseq

"""


class AssembleIter(TaskList):
    class Assemble(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # TODO: Ensure that this override step on wdir functions as expected
            self._wdir = os.path.join(os.path.dirname(self.wdir), "transcriptomes")
            files_dict = AssembleIter.Assemble.get_reads_files(self.rnaseq)
            built_files = []
            for files_tuple in files_dict[self.record_id]:
                if not all(os.path.exists(_f) for f in files_tuple for _f in f):
                    raise FileNotFoundError
                for file_tuple in files_tuple:
                    built_files.append(os.path.join(self.wdir, prefix(file_tuple[0]) + ".transcriptome.fna"))
            self.output = [
                *self.input,
                *built_files,  # Assembled transcriptomes
                built_files,  # List of files
            ]

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            # Get reads files
            files_dict = AssembleIter.Assemble.get_reads_files(self.rnaseq)
            for files_tuple in files_dict[self.record_id]:
                if not all(os.path.exists(_f) for f in files_tuple for _f in f):
                    raise FileNotFoundError
                for file_tuple in files_tuple:
                    _rna_basename = prefix(file_tuple[0])
                    out_file = os.path.join(self.wdir, _rna_basename + ".transcriptome.fna")
                    # Only assemble if not already done so
                    if os.path.exists(out_file):
                        continue
                    # Run Trinity
                    self.log_and_run(
                        self.program[
                            "--seqType", "fq",
                            (*self.added_flags),
                            "--left", file_tuple[0], "--right", file_tuple[1],
                            "--CPU", self.threads,
                            "--output", os.path.join(
                                self.wdir, "trinity" + "_" + _rna_basename
                            ),
                        ]
                    )
                    # Rename output file
                    os.replace(
                        os.path.join(
                            self.wdir,
                            "trinity" + "_" + _rna_basename,
                            "Trinity.fasta"
                        ),
                        out_file
                    )

        # Line format:
        # genome1\tSRR.1.fq,SRR.2.fq
        # genome2\tSRR.1.fq,SRR.2.fq
        # genome3\tSRR.1.fq,SRR.2.fq;SRR.3.fq,SRR.4.fq
        @staticmethod
        def get_reads_files(file_path: str) -> DefaultDict[str, List[Tuple[str, str]]]:
            reads_files = defaultdict(list)
            if not os.path.exists(file_path):
                return reads_files
            fp = open(file_path, "r")
            for line in fp:
                line = line.rstrip("\r\n").split("\t")
                reads_files[line[0]].append(
                    tuple([_v for _v in val.split(",") if _v != ""])
                    for val in line[1].split(";")
                    if val != ""
                )
            return reads_files

    def __init__(self, *args, **kwargs):
        super().__init__(AssembleIter.Assemble, "assemble", *args, **kwargs)


if __name__ == "__main__":
    pass
