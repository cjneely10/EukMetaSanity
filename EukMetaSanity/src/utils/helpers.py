import logging

"""
General functions

"""


# Function logs and runs dask command
def log_and_run(cmd, test: int):
    logging.info(str(cmd))
    if test == 1:
        cmd()


def touch(_path: str):
    open(_path, "w").close()
