import os
from Bio import SeqIO
from EukMetaSanity.tasks.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch


class FindRNAIter(TaskList):
    class FindRNA(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [
                *self.input,
                os.path.join(self.wdir, self.record_id + ".tblout"),
            ]
        
        def run(self):
            super().run()
            
        @program_catch
        def run_1(self):
            # Run cmscan
            self.log_and_run(
                self.program_cmscan[
                    "-Z", str(FindRNAIter.FindRNA.calculate_z(self.input[0])),
                    (*self.added_flags),
                    "--tblout", os.path.join(self.wdir, self.record_id + ".all.tblout"),
                    "--clanin", self.data_clanin,
                    "--cpu", self.threads,
                    self.data_cm,
                    self.input[0],
                ] > os.path.join(self.wdir, self.record_id + ".cmscan")
            )
            # Save non-overlapping hits
            self.log_and_run(
                self.local["grep"][
                    "-v", " = ",
                    os.path.join(self.wdir, self.record_id + ".all.tblout")
                ] > os.path.join(self.wdir, self.record_id + ".tblout")
            )

        @staticmethod
        def calculate_z(fasta_file: str) -> float:
            total: float = 0.0
            fp = SeqIO.parse(fasta_file, "fasta")
            for record in fp:
                total += float(len(record.seq))
            total *= 2.0
            return total / 1E6
            
    def __init__(self, *args, **kwargs):
        super().__init__(FindRNAIter.FindRNA, "rfam", *args, **kwargs)


if __name__ == "__main_":
    pass
