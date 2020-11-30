import os
import shutil
from Bio import SeqIO
from pathlib import Path
from collections import Counter
from EukMetaSanity import InvalidProtocolError
from EukMetaSanity import Task, TaskList, program_catch, prefix
from EukMetaSanity.tasks.test_tasks.run.helpers.taxonomy import get_taxonomy
from EukMetaSanity.tasks.test_tasks.run.helpers.abinitio import augustus_taxon_ids


class AbInitioIter(TaskList):
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

        def augustus(self):
            # Initial training based on best species from taxonomy search
            out_gff = self._augustus(self._augustus_tax_ident(), 1, str(self.input["root"]["fna"]))
            self._train_augustus(1, str(self.input["root"]["fna"]), out_gff)
            # Remaining rounds of re-training on generated predictions
            for i in range(self.rounds):
                _last = False
                if i == self.rounds - 1:
                    _last = True
                out_gff = self._augustus(self.record_id + str(i + 1), i + 2, str(self.input["root"]["fna"]), _last)
                if i != self.rounds - 1:
                    self._train_augustus(i + 2, str(self.input["root"]["fna"]), out_gff)
            # Move any augustus-generated config stuff
            self._handle_config_output()
            # Rename final file
            os.replace(out_gff, str(self.output["ab-gff3"]))

        def _augustus(self, species: str, _round: int, _file: str, _last: bool = False):
            out_gff = os.path.join(
                self.wdir, AbInitioIter.AbInitio._out_path(str(self.input["root"]["fna"]), ".%i.gb" % _round)
            )
            # Chunk file predictions
            record_p = SeqIO.parse(_file, "fasta")
            progs = []
            out_files = []
            out_gffs = []
            for record in record_p:
                out_file_path = os.path.join(self.wdir, str(record.id) + ".fna")
                _out_gff = out_gff + "-" + str(record.id)
                out_files.append(out_file_path)
                out_gffs.append(_out_gff)
                SeqIO.write([record], out_file_path, "fasta")
                # Run prediction
                progs.append(
                    self.program_augustus[
                        "--codingseq=on",
                        "--stopCodonExcludedFromCDS=false",
                        "--species=%s" % species,
                        "--outfile=%s" % _out_gff,
                        ("--gff3=on" if _last else "--gff3=off"),
                        out_file_path,
                    ]
                )
            self.batch(progs)
            # Combine files
            (
                self.local["cat"][out_gffs] |
                self.local["gffread"]["-o", out_gff + ".tmp", "-F", "-G", "--keep-comments"]
            )()
            # Make ids unique
            self._make_unique(out_gff)
            # Remove intermediary files
            all([os.remove(_file) for _file in out_files])
            all([os.remove(_file) for _file in out_gffs])
            return out_gff

        def _handle_config_output(self):
            # Move the augustus training folders to their wdir folders
            config_dir = os.path.join(
                os.path.dirname(os.path.dirname(Path(str(self.program_augustus)).resolve())),
                "config", "species"
            )
            for i in range(1, int(self.rounds)):
                shutil.move(
                    os.path.join(config_dir, self.record_id + str(i + 1)),
                    self.wdir
                )

        @program_catch
        def _augustus_tax_ident(self) -> str:
            tax_db = os.path.join(self.wdir, self.record_id + "-augustus_db")
            if not os.path.exists(tax_db + ".m8"):
                # Run taxonomy search
                self.parallel(
                    self.program_mmseqs[
                        "linsearch",
                        str(self.input["taxonomy"]["seq_db"]),  # Input FASTA sequence db
                        self.data,  # Input augustus-db
                        tax_db,  # Output tax db
                        os.path.join(self.wdir, "tmp"),
                        "-c", "0.6",
                        "--cov-mode", "0",
                        "-e", "0.01",
                        "--remove-tmp-files",
                        "--threads", self.threads,
                    ]
                )
                # Output results
                self.parallel(
                    self.program_mmseqs[
                        "convertalis",
                        str(self.input["taxonomy"]["seq_db"]),  # Input FASTA sequence db
                        self.data,  # Input augustus-db
                        tax_db,  # Input tax db
                        tax_db + ".m8",  # Output results file
                        "--threads", self.threads,
                        "--format-output", "query,target,pident,taxid,taxname,taxlineage",
                    ],
                    "2:00:00"
                )
            # Return optimal taxonomy
            augustus_ids_dict = augustus_taxon_ids()
            found_taxa = Counter()
            with open(tax_db + ".m8", "r") as R:
                for line in R:
                    line = line.rstrip("\r\n").split()
                    # Count those that pass the user-defined cutoff value
                    if line[3] in augustus_ids_dict.keys() and float(line[2]) >= (float(self.cutoff) / 100.):
                        found_taxa[line[3]] += 1
            return augustus_ids_dict[found_taxa.most_common()[0][0]]

        @program_catch
        def _train_augustus(self, _round: int, _file: str, out_gff: str):
            # Remove old training directory, if needed
            config_dir = os.path.join(
                os.path.dirname(os.path.dirname(Path(str(self.program_augustus)).resolve())),
                "config", "species", self.record_id + str(_round)
            )
            if os.path.exists(config_dir):
                shutil.rmtree(config_dir)
            # Parse to genbank
            out_gb = os.path.join(self.wdir, AbInitioIter.AbInitio._out_path(_file, ".%i.gb" % _round))
            self.local["gff2gbSmallDNA.pl"][
                out_gff,
                _file,
                "1000",
                out_gb
            ]()

            species_config_prefix = self.record_id + str(_round)
            # Write new species config file
            self.program_new_species_pl[
                "--species=%s" % species_config_prefix,
                out_gb
            ]()
            # Run training
            self.program_etraining[
                "--species=%s" % species_config_prefix,
                out_gb
            ]()
            return out_gff

        @staticmethod
        def _make_unique(out_gff):
            gff_fp = open(out_gff + ".tmp", "r")
            out_fp = open(out_gff, "w")
            i = 1
            line = next(gff_fp)
            while True:
                if line.startswith("#"):
                    out_fp.write(line)
                else:
                    line = line.split("\t")
                    if line[2] == "transcript":
                        out_fp.write("\t".join((
                            line[0],
                            "ab-initio",
                            *line[2:-1],
                            "ID=gene%i\n" % i
                        )))
                        try:
                            line = next(gff_fp).split("\t")
                        except StopIteration:
                            break
                        while line[2] != "transcript":
                            out_fp.write("\t".join((
                                line[0],
                                "ab-initio",
                                *line[2:-1],
                                "Parent=gene%i\n" % i
                            )))
                            try:
                                line = next(gff_fp).split("\t")
                            except StopIteration:
                                break
                        i += 1
                try:
                    line = next(gff_fp)
                except StopIteration:
                    break
            out_fp.close()

        @staticmethod
        def _out_path(_file_name: str, _ext: str) -> str:
            return os.path.basename(os.path.splitext(_file_name)[0]) + _ext

        @program_catch
        def gmes(self):
            tax = get_taxonomy(str(self.input["root"]["fna"]), 2.0, "species")
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
            # Run prothint
            self.parallel(
                self.program_prothint[
                    str(self.input["root"]["fna"]),
                    _fasta_output,
                    "--workdir", self.wdir,
                    "--threads", self.threads,
                ]
            )
            script = self.create_script(
                self.program_gmes[
                    "--sequence", str(self.input["root"]["fna"]),
                    "--EP", os.path.join(self.wdir, "prothint.gff"),
                    "--evidence", os.path.join(self.wdir, "evidence.gff"),
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
            self.local["rm"][_fasta_output]()
            for file in os.listdir(self.wdir):
                file = os.path.join(self.wdir, file)
                if subset_db_outpath in file:
                    self.local["rm"][file]()
            
    def __init__(self, *args, **kwargs):
        super().__init__(AbInitioIter.AbInitio, AbInitioIter.name, *args, **kwargs)


if __name__ == "__main_":
    pass
