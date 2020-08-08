import os
from EukMetaSanity.tasks.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch

"""
Map RNAseq to masked genome

"""


class RnaSeqIter(TaskList):
    class RnaSeq(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [*self.input]

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            # Get pairs to align
            read_pairs = self.get_rna_read_pairs()
            if read_pairs is None:
                self.output = []
                return
            # Generate genome index
            genome_prefix = os.path.join(self.wdir, prefix(self.input[0]) + "_db")
            self.log_and_run(
                self.program_hisat2build[
                    self.input[0],
                    genome_prefix
                ]
            )
            out = []
            for pair in read_pairs:
                out_prefix = os.path.join(self.wdir, prefix(pair[0]))
                # Align
                self.log_and_run(
                    self.program_hisat2[
                        "-p", self.threads,
                        "-x", genome_prefix,
                        "-1", pair[0],
                        "-2", pair[1],
                        "-S", out_prefix + ".sam",
                        (*self.added_flags),
                    ]
                )
                RnaSeqIter.RnaSeq.sambamba(self, out_prefix)
                # Store path to file in new output
                out.append(out_prefix + ".sorted.bam")
            self.output = [
                *self.output,  # Forward original data
                *out,  # Paths
                out,  # List of data
            ]

        def get_rna_read_pairs(self) -> [(str, str)]:
            if not os.path.exists(self.rnaseq):
                return
            fp = open(self.rnaseq, "r")
            for line in fp:
                if self.record_id in line:
                    pairs_string = line.rstrip("\r\n").split("\t")[1].split(";")
                    return [(p[0], p[1]) for pair in pairs_string for p in pair.split(",")]

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
            task_object.local["rm"][out_prefix + ".{sam,bam}"]()

    def __init__(self, *args, **kwargs):
        super().__init__(RnaSeqIter.RnaSeq, "rnaseq", *args, **kwargs)


if __name__ == "__main__":
    pass
