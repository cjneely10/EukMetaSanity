import os
from EukMetaSanity.tasks.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.run.initial_evidence import EvidenceIter


class MergeIter(TaskList):
    class Merge(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            _out = {
                "nr_gff3": os.path.join(self.wdir, self.record_id) + ".nr.gff3",
                "mask": self.input[0],
                "prot": os.path.join(self.wdir, self.record_id) + ".faa",
                "cds": os.path.join(self.wdir, self.record_id) + ".cds.fna",
            }
            self.output = [
                _out,  # Dictionary for accessing to write final summary
                *list(_out.values()),  # Regular list of values for final path checking
            ]
        
        def run(self):
            super().run()
            
        @program_catch
        def run_1(self):
            gff3s = []
            all_files = self.input[-2] + self.input[-1]
            for sorted_bam in all_files:
                print(sorted_bam)
                out_prefix = os.path.join(self.wdir, prefix(sorted_bam))
                if not os.path.exists(out_prefix + ".bed"):
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
                        "merge"
                        "-i", out_prefix + ".bed",
                        "-c", "1", "-o", "count"
                    ] > out_prefix + ".tmp.merged.bed"
                )
                # Remove regions that do not have enough coverage
                self.log_and_run(
                    self.local["awk"][
                        '$4 > %s' % self.min_depth,
                        out_prefix + ".tmp.merged.bed"
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
