import logging


def log_and_run(cmd):
    logging.info(str(cmd))
    cmd()
