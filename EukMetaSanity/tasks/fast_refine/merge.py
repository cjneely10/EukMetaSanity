import os
from EukMetaSanity.tasks.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.run.initial_evidence import EvidenceIter


class MergeIter(TaskList):
    class Merge(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                *self.output,
                os.path.join(self.wdir, self.record_id) + ".nr.gff3"
            ]
        
        def run(self):
            super().run()
            
        @program_catch
        def run_1(self):
            gff3s = []
            for sorted_bam in self.input[-1]:
                out_prefix = os.path.join(self.wdir, prefix(sorted_bam))
                # Convert to BED
                self.log_and_run(
                    self.program_bedtools[
                        "bamtobed",
                        "-i", sorted_bam
                    ] > out_prefix + ".bed"
                )
                # Merge overlapping reads
                self.log_and_run(
                    self.program_bedtools[
                        "merge",
                        "-i", out_prefix + ".bed"
                    ] > out_prefix + ".merged.bed"
                )
                # Cluster to putative CDS regions
                self.log_and_run(
                    self.program_bedtools[
                        "cluster",
                        "-i", out_prefix + ".merged.bed",
                        (*self.added_flags),
                    ] > out_prefix + ".clustered.bed"
                )
                # Convert to gff3
                self.log_and_run(
                    self.local["bed-to-gff3.py"][
                        out_prefix + ".clustered.bed",
                        self.input[0],
                        "-o", out_prefix + ".gff3"
                    ]
                )
                gff3s.append(out_prefix + ".gff3")
            # Merge all results
            EvidenceIter.Evidence.merge(
                self,
                [*gff3s, self.input[1]],
                self.input[0],
                os.path.join(self.wdir, self.record_id),
            )
            
    def __init__(self, *args, **kwargs):
        super().__init__(MergeIter.Merge, "merge", *args, **kwargs)


if __name__ == "__main_":
    pass
