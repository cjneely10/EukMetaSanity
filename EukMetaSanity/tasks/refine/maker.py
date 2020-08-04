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
            if not os.path.exists("maker_opts.ctl"):
                self.log_and_run(self.program["-CTL"])
                self.log_and_run(self.local["cp"]["*.ctl", self.wdir])
            opts_file = os.path.join(self.wdir, "maker_opts.ctl")
            # Add proper paths to file
            # Edit path to genome
            self.log_and_run(
                self.local["sed"][
                    "-i", "s/%s/%s/" % ("genome=", "genome=%s" % self.input[4]),
                    opts_file,
                ]
            )
            # Incorporate search for simple repeats with MAKER, use existing models from Run
            self.log_and_run(
                self.local["sed"][
                    "-i", "s/%s/%s/" % ("model_org=", "model_org=simple"),
                    opts_file,
                ]
            )
            self.log_and_run(
                self.local["sed"][
                    "-i",
                    "s/%s/%s/" % (
                        "rm_gff=", "rm_gff=%s" % os.path.join(self.wdir, self.record_id + ".mask.complex.reformat.gff3")
                    ),
                    opts_file,
                ]
            )
            # Add metaeuk evidence
            self.log_and_run(
                self.local["sed"][
                    "-i", "s/%s/%s/" % ("other_gff=", "other_gff=%s" % self.input[5]),
                    opts_file,
                ]
            )
            # Parse additional evidence as needed
            prot_file = self.get_data_file(self.proteins)
            if prot_file is not None:
                self.log_and_run(
                    self.local["sed"][
                        "-i", "s/%s/%s/" % ("protein=", "protein=%s" % prot_file),
                        opts_file
                    ]
                )
                self.log_and_run(
                    self.local["sed"][
                        "-i", "s/protein2genome=0/protein2genome=1/",
                        opts_file
                    ]
                )
            est_file = self.get_data_file(self.est)
            if est_file is not None:
                self.log_and_run(
                    self.local["sed"][
                        "-i", "s/%s/%s/g" % ("est=", "est=%s" % est_file),
                        opts_file
                    ]
                )
                self.log_and_run(
                    self.local["sed"][
                        "-i", "s/est2genome=0/est2genome=1/",
                        opts_file
                    ]
                )
            # Parse user args into config file

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

        def get_data_file(self, file_path: str):
            if not os.path.exists(file_path):
                return
            with open(file_path, "r") as w:
                for line in w:
                    if line.startswith(self.record_id):
                        return line.rstrip("\r\n").split("\t")[1]

    def __init__(self, *args, **kwargs):
        super().__init__(MakerIter.Maker, "maker", *args, **kwargs)


if __name__ == "__main__":
    pass
