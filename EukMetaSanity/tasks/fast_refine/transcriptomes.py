import os
from typing import List, Optional

from EukMetaSanity.tasks.fast_refine.rnaseq import RnaSeqIter
from EukMetaSanity.tasks.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch


class TranscriptomesIter(TaskList):
    class Transcriptomes(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = [*self.input]
        
        def run(self):
            super().run()
            
        @program_catch
        def run_1(self):
            # Get transcripts
            transcripts = self.get_transcripts()
            if transcripts is None:
                self.output = [*self.input]
                return
            # Generate genome index
            genome_idx = prefix(self.input[0]) + "_db"
            _genome_dir = os.path.dirname(self.input[0])
            _genome_basename = os.path.basename(self.input[0])
            self.log_and_run(
                self.program_gmapbuild[
                    "-d", genome_idx,
                    "-D", _genome_dir, _genome_basename
                ]
            )
            out = []
            for transcript in transcripts:
                out_prefix = os.path.join(self.wdir, prefix(transcript))
                # Align
                self.log_and_run(
                    self.program_gmap[
                        "-D", _genome_dir, "-d", genome_idx,
                        "-t", self.threads,
                        "-f", "samse", transcript
                    ] > out_prefix + ".sam"
                )
                # Run sambamba
                RnaSeqIter.RnaSeq.sambamba(self, out_prefix)
                # Store path to file in new output
                out.append(out_prefix + ".sorted.bam")
            self.output = [
                *self.output,  # Forward original data
                *out,  # Paths
                out,  # List of data
            ]

        def get_transcripts(self) -> Optional[List[str]]:
            if not os.path.exists(self.transcriptomes):
                return
            fp = open(self.transcriptomes, "r")
            for line in fp:
                if self.record_id in line:
                    return [p for p in line.rstrip("\r\n").split("\t")[1].split(";") if p != ""]
            
    def __init__(self, *args, **kwargs):
        super().__init__(TranscriptomesIter.Transcriptomes, "transcriptomes", *args, **kwargs)


if __name__ == "__main_":
    pass
