import logging

"""
General functions

"""


# Function logs and runs dask command
def log_and_run(cmd):
    logging.info(str(cmd))
    cmd()
