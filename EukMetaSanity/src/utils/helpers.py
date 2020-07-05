"""
General functions

"""


# Mimic touch command from linux
def touch(_path: str):
    open(_path, "w").close()


# Get prefix of path - e.g. for /path/to/file_1.ext, return file_1
def prefix(_path: str) -> str:
    return ".".join(_path.split("/")[-1].split(".")[:-1])
