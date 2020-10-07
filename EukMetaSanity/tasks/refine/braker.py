import os
from EukMetaSanity.tasks.utils.helpers import prefix
from EukMetaSanity import Task, TaskList, program_catch
from EukMetaSanity.tasks.run.taxonomy import TaxonomyIter
from EukMetaSanity.tasks.run.initial_evidence import EvidenceIter


class BrakerIter(TaskList):
    class Braker(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            bams = ""
            if len(self.input) > 2:
                bams = (",".join((*self.input[-1], *self.input[-2])))
            if bams != "":
                bams = ["--bam=%s" % (",".join((*self.input[-1], *self.input[-2])))]
            else:
                bams = []
            _record_id = self.record_id.replace("-mask", "")
            if len(bams) > 0:
                if os.path.exists(os.path.join(self.wdir, _record_id + ".gff3")):
                    os.remove(os.path.join(self.wdir, _record_id + ".gff3"))
                out = {
                    "cds": os.path.join(self.wdir, _record_id + ".cds.fna"),
                    "prot": os.path.join(self.wdir, _record_id + ".faa"),
                    "nr_gff3": os.path.join(self.wdir, _record_id + ".gff3"),
                    "fna": self.input[2],
                }
                self.output = [
                    out,  # Paths
                    *list(out.values()),  # List of data
                ]
            else:
                self.output = []

        @program_catch
        def run_1(self):
            bams = ""
            if len(self.input) > 2:
                bams = (",".join((*self.input[-1], *self.input[-2])))
            if bams != "":
                bams = ["--bam=%s" % (",".join((*self.input[-1], *self.input[-2])))]
            else:
                bams = []
            if len(bams) > 0:
                tax = TaxonomyIter.Taxonomy.get_taxonomy(self.input[1], 0, self.level)
                _tax = []
                if "fungi" in tax[0]:
                    _tax = ["--fungus"]
                subset_db_outpath = os.path.join(self.wdir, self.record_id + "-tax-prots_%s" % prefix(self.data))
                _fasta_output = os.path.join(self.wdir, self.record_id + ".faa")
                if not os.path.exists(_fasta_output):
                    self.log_and_run(
                        self.program_mmseqs[
                            "filtertaxseqdb",
                            self.data,
                            subset_db_outpath,
                            "--taxon-list", tax[1],
                            "--threads", self.threads,
                        ]
                    )
                    # Output as FASTA file
                    self.log_and_run(
                        self.program_mmseqs[
                            "convert2fasta",
                            subset_db_outpath,
                            _fasta_output,
                        ]
                    )
                _added = []
                if "--skipGeneMark-ES" in self.added_flags:
                    # Output CDS and prot fasta
                    _gtf = os.path.join(self.wdir, self.record_id + ".gtf")
                    self.log_and_run(
                        self.program_gffread[self.input[3]] > _gtf
                    )
                    _added = ["--geneMarkGtf=%s" % _gtf]
                if not os.path.exists(os.path.join(self.wdir, "GeneMark-ET", "genemark.gtf")) or \
                        not os.path.exists(os.path.join(self.wdir, "augustus.hints.gtf")):
                    self.log_and_run(
                        self.program_braker[
                            "--cores=%s" % str(self.threads),
                            "--genome=%s" % self.input[0],
                            (*bams),
                            "--prot_seq=%s" % _fasta_output,
                            "--workingdir=%s" % self.wdir,
                            (*_tax),
                            "--species=%s" % self.record_id,
                            (*self.added_flags),
                            "--prg=exonerate",
                            (*_added)
                        ]
                    )
                _record_id = self.record_id.replace("-mask", "")
                _files = [
                    os.path.join(self.wdir, "GeneMark-ET", "genemark.gtf"),
                    os.path.join(self.wdir, "augustus.hints.gtf"),
                ]
                _out = os.path.join(self.wdir, "braker.gtf")
                if "--gff3" in self.added_flags:
                    _files = [
                        os.path.join(self.wdir, "GeneMark-ET", "genemark.gff3"),
                        os.path.join(self.wdir, "augustus.hints.gff3"),
                    ]
                    _out = os.path.join(self.wdir, "braker.gff3"),
                # Merge final results if final output failed
                if not os.path.exists(_out):
                    self.log_and_run(
                        self.program_gffread[
                            (*_files),
                            "--merge", "-G", "-S",
                            "-o", os.path.join(self.wdir, _record_id + ".gff3"),
                            "-g", self.input[2],
                            "-x", os.path.join(self.wdir, _record_id + ".cds.fna"),
                            "-y", os.path.join(self.wdir, _record_id + ".faa"),
                        ]
                    )
                # Else just keep
                else:
                    self.output[0]["braker"] = _out

    def __init__(self, *args, **kwargs):
        super().__init__(BrakerIter.Braker, "braker", *args, **kwargs)


if __name__ == "__main_":
    pass
