import os
import shutil
from Bio import SeqIO
from pathlib import Path
from collections import Counter
from plumbum import ProcessExecutionError
from EukMetaSanity import InvalidProtocolError
from EukMetaSanity import Task, TaskList, program_catch, prefix
from EukMetaSanity.tasks.official.run.helpers.taxonomy import get_taxonomy
from EukMetaSanity.tasks.official.run.helpers.abinitio import augustus_taxon_ids


class AbInitioIter(TaskList):
    """ Task runs either Augustus or GeneMark ab initio prediction pipelines on genome file

    Outputs: ab-gff3
    Finalizes: ab-gff3

    """
    name = "abinitio"
    requires = ["taxonomy", "repeats"]
    
    class AbInitio(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.output = {
                "ab-gff3": os.path.join(self.wdir, self.record_id + ".gff3"),
                "final": ["ab-gff3"]
            }
            self.rounds = int(self.rounds)
            
        @program_catch
        def run(self):
            # Use augustus predictor
            if self.protocol == "augustus":
                self.augustus()
            # Use GeneMark predictor
            elif self.protocol == "gmes":
                self.gmes()
            else:
                raise InvalidProtocolError

        @program_catch
        def gmes(self):
            tax = get_taxonomy(str(self.input["taxonomy"]["tax-report"]), 2.0, "order")
            subset_db_outpath = os.path.join(self.wdir, self.record_id + "-tax-prots_%s" % prefix(self.data))
            _fasta_output = os.path.join(self.wdir, self.record_id + ".faa")
            self.parallel(
                self.program_mmseqs[
                    "filtertaxseqdb",
                    self.data,
                    subset_db_outpath,
                    "--taxon-list", tax[1],
                    "--threads", self.threads,
                ],
                "1:00:00"
            )
            # Output as FASTA file
            self.single(
                self.program_mmseqs[
                    "convert2fasta",
                    subset_db_outpath,
                    _fasta_output,
                ]
            )
            try:
                # Run prothint
                self.parallel(
                    self.program_prothint[
                        str(self.input["root"]["fna"]),
                        _fasta_output,
                        "--workdir", self.wdir,
                        "--threads", self.threads,
                    ]
                )
            except ProcessExecutionError as e:
                pass
            ev_vals = ["--ES"]
            if os.path.exists(os.path.join(self.wdir, "prothint.gff")):
                ev_vals = ["--EP", os.path.join(self.wdir, "prothint.gff"),
                           "--evidence", os.path.join(self.wdir, "evidence.gff")]
            script = self.create_script(
                self.program_gmes[
                    "--sequence", str(self.input["root"]["fna"]),
                    (*ev_vals),
                    "--cores", self.threads, (*self.added_flags),
                    ("--fungus" if "fungi" == tax[0] else "")
                ],
                "abinitio.sh"
            )
            # Run script
            self.parallel(script)
            # Move program to match required output name
            self.single(
                self.program_gffread[
                    os.path.join(self.wdir, "genemark.gtf"), "-G",
                    "-o", str(self.output["ab-gff3"])
                ]
            )
            self.single(
                self.local["sed"][
                    "-i", "s/GeneMark.hmm/ab-initio/g",
                    str(self.output["ab-gff3"])
                ]
            )
            if os.path.exists(_fasta_output):
                self.local["rm"][_fasta_output]()
            for file in os.listdir(self.wdir):
                file = os.path.join(self.wdir, file)
                if subset_db_outpath in file:
                    self.local["rm"][file]()
            
    def __init__(self, *args, **kwargs):
        super().__init__(AbInitioIter.AbInitio, AbInitioIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
