import os
from pathlib import Path
from typing import List, Optional
from EukMetaSanity.tasks.official.refine.rnaseq import RnaSeqIter
from EukMetaSanity.tasks.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch


class TranscriptomesIter(TaskList):
    class Transcriptomes(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            out = []
            transcripts = self.get_transcripts()
            for transcript in transcripts:
                out.append(os.path.join(self.wdir, prefix(transcript)) + ".sorted.bam")
            self.output = [
                *self.input[:-1],  # Forward original data
                *out,  # Paths
                self.input[-1],  # aligned reads
                out,  # List of data
            ]
            
        @program_catch
        def run_1(self):
            # Get transcripts
            transcripts = self.get_transcripts()
            if len(transcripts) > 0:
                # Generate genome index
                genome_idx = self.record_id + "_db"
                _genome_dir = os.path.dirname(self.input[0])
                _genome_basename = os.path.basename(self.input[0])
                self.log_and_run(
                    self.program_gmapbuild[
                        "-d", genome_idx,
                        "-D", _genome_dir, _genome_basename
                    ]
                )
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

        def get_transcripts(self) -> Optional[List[str]]:
            _path = str(Path(self.transcriptomes).resolve())
            if not os.path.exists(_path):
                return []
            fp = open(_path, "r")
            _id = self.record_id.replace("-mask", "")
            for line in fp:
                if _id in line:
                    return line.rstrip("\r\n").split("\t")[1].split(",")
            return []
            
    def __init__(self, *args, **kwargs):
        super().__init__(TranscriptomesIter.Transcriptomes, "transcriptomes", *args, **kwargs)


if __name__ == "__main_":
    pass
