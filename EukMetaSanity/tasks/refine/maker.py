import os
from EukMetaSanity import Task, TaskList, program_catch


class MakerIter(TaskList):
    class Maker(Task):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def run(self):
            super().run()

        @program_catch
        def run_1(self):
            self.reformat_repeats_gff3()
            self.generate_ctl_file()

        def generate_ctl_file(self):
            # Create base file and move to wdir
            self.log_and_run(self.program["-CTL"])
            self.log_and_run(self.local["cp"]["*.ctl", self.wdir])
            opts_file = os.path.join(self.wdir, "maker_opts.ctl")
            # Add proper paths to file
            # Edit path to genome
            self.log_and_run(
                self.local["sed"][
                    "-i", "s/%s/%s/g" % ("genome=", "genome=%s" % self.input[4]),
                    opts_file,
                ]
            )
            # Incorporate search for simple repeats with MAKER, use existing models from Run
            self.log_and_run(
                self.local["sed"][
                    "-i", "s/%s/%s/g" % ("model_org=", "model_org=simple"),
                    opts_file,
                ]
            )
            self.log_and_run(
                self.local["sed"][
                    "-i",
                    "s/%s/%s/g" % (
                        "rm_gff=", "rm_gff=%s" % os.path.join(self.wdir, self.record_id + ".mask.complex.reformat.gff3")
                    ),
                    opts_file,
                ]
            )
            # Add metaeuk evidence
            self.log_and_run(
                self.local["sed"][
                    "-i", "s/%s/%s/g" % ("other_gff=", "other_gff=%s" % self.input[5]),
                    opts_file,
                ]
            )
            # Parse additional evidence as needed

        def reformat_repeats_gff3(self):
            # Isolate complex repeats
            self.log_and_run(
                self.local["grep"][
                    "-v", "-e", "Satellite", "-e", ")n", "-e", "-rich",
                    self.input[3]
                ] > os.path.join(self.wdir, self.record_id + ".mask.complex.gff3")
            )
            # Reformat to work with MAKER
            self.log_and_run(
                self.local["cat"][os.path.join(self.wdir, self.record_id + ".mask.complex.gff3")] |
                self.local["perl"][
                    "-ane",
                    r'$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", '
                    r'@F)."\n"} print $_',
                ] > os.path.join(self.wdir, self.record_id + ".mask.complex.reformat.gff3")
            )

        @staticmethod
        def get_proteins_file(file_path: str) -> str:
            pass

        @staticmethod
        def get_est_file(file_path: str) -> str:
            pass

    def __init__(self, *args, **kwargs):
        super().__init__(MakerIter.Maker, "maker", *args, **kwargs)


if __name__ == "__main__":
    pass
