import os
from pathlib import Path
from typing import List, Optional, Tuple
from EukMetaSanity.tasks.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch

"""
Map RNAseq to masked genome

"""


class RnaSeqIter(TaskList):
    class RnaSeq(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            out = []
            read_pairs = self.get_rna_read_pairs()
            for pair in read_pairs:
                out.append(os.path.join(self.wdir, prefix(pair[0])) + ".sorted.bam")
            self.output = [
                *self.input,  # Forward original data
                *out,  # Paths
                out,  # List of data
            ]

        @program_catch
        def run_1(self):
            # Get pairs to align
            read_pairs = self.get_rna_read_pairs()
            if len(read_pairs) > 0:
                # Generate genome index
                genome_prefix = os.path.join(self.wdir, self.record_id + "_db")
                self.log_and_run(
                    self.program_hisat2build[
                        self.input[0],
                        genome_prefix
                    ]
                )
                for pair in read_pairs:
                    out_prefix = os.path.join(self.wdir, prefix(pair[0]))
                    if len(pair) > 1:
                        _parse_args = ["-1", pair[0], "-2", pair[1]]
                    else:
                        _parse_args = ["-U", pair[0]]
                    # Align
                    self.log_and_run(
                        self.program_hisat2[
                            "-p", self.threads,
                            "-x", genome_prefix,
                            (*_parse_args),
                            "-S", out_prefix + ".sam",
                            "--novel-splicesite-outfile", out_prefix + ".splicesites",
                            (*self.added_flags),
                        ]
                    )
                    self.log_and_run(
                        self.local["awk"][
                            "{print $1\"\\tRNA_seq_junction\\tintron\\t\"$2\"\\t\"$3\"\\t500\\t\"$4\"\\t.\\t.\"}",
                            out_prefix + ".splicesites"
                        ] > out_prefix + ".hints.gff"
                    )
                    # Run sambamba
                    RnaSeqIter.RnaSeq.sambamba(self, out_prefix)

        def get_rna_read_pairs(self) -> Optional[List[Tuple[str, str]]]:
            _path = str(Path(self.rnaseq).resolve())
            if not os.path.exists(_path):
                return []
            fp = open(_path, "r")
            _id = self.record_id.replace("-mask", "")
            for line in fp:
                if line.startswith(_id):
                    pairs_string = line.rstrip("\r\n").split("\t")[1].split(";")
                    return [tuple(pair.split(",")) for pair in pairs_string]
            return []

        @staticmethod
        def sambamba(task_object: Task, out_prefix: str):
            # Convert to sorted bam
            task_object.log_and_run(
                task_object.program_sambamba[
                    "view",
                    "-S", out_prefix + ".sam",
                    "-f", "bam",
                    "-t", task_object.threads,
                    "-o", out_prefix + ".bam",
                ]
            )
            task_object.log_and_run(
                task_object.program_sambamba[
                    "sort",
                    "-t", task_object.threads,
                    "-m", task_object.sambamba_memlimit,
                    "-o", out_prefix + ".sorted.bam",
                    out_prefix + ".bam"
                ]
            )
            # Remove intermediary files
            task_object.local["rm"][out_prefix + ".sam", out_prefix + ".bam"]()

    def __init__(self, *args, **kwargs):
        super().__init__(RnaSeqIter.RnaSeq, "rnaseq", *args, **kwargs)


if __name__ == "__main__":
    pass
