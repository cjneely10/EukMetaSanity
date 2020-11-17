import os
from typing import List, Dict
from plumbum.machines.local import LocalCommand, LocalMachine


class SLURMCaller:

    OUTPUT_SCRIPTS = "slurm-runner.sh"

    def __init__(self, wdir: str, cmd: LocalCommand, config_data: Dict[str, str], _local: LocalMachine, slurm_flags: List[str]):
        """

        :param wdir: Working directory to write slurm script
        :param cmd: Command to run within script
        :param config_data: Dict of data to parse for sbatch script header
        :param _local: local object from which to get sbatch command object
        :param slurm_flags: Added flags to include in script, if any
        """
        self.cmd = cmd
        self.wdir = wdir
        self.config = config_data
        self.local = _local
        self.added_flags = slurm_flags

        # Generated job id
        self.job_id = -1
        # Path of script that is running
        self.script = os.path.join(self.wdir, SLURMCaller.OUTPUT_SCRIPTS)
        # Create slurm script in working directory
        self.generate_script()

    def __str__(self):
        return open(os.path.join(self.wdir, SLURMCaller.OUTPUT_SCRIPTS)).read()

    def __repr__(self):
        return self.__str__()

    # Run script written
    def launch_script(self) -> bool:
        pass

    def is_running(self) -> bool:
        pass

    def has_launched(self) -> bool:
        pass

    def has_ended(self) -> bool:
        pass

    def parse_stdout_for(self) -> str:
        return ""

    def generate_script(self):
        fp = open(self.script, "w")
        fp.write("#!/bin/bash\n")

        fp.close()

    @staticmethod
    def create_header_line(param: str, arg: str) -> str:
        return "#SBATCH %s=%s" % (param, arg)

    # Call will run script using sbatch and periodically check for when to stop running
    def __call__(self, *args, **kwargs):
        pass
