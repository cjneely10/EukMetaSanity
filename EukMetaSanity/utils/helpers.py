import os
from typing import Tuple, List

"""
General functions

"""


# Mimic touch command from linux
def touch(_path: str):
    open(_path, "w").close()


# Get prefix of path - e.g. for /path/to/file_1.ext, return file_1
def prefix(_path: str) -> str:
    return os.path.basename(os.path.splitext(_path)[0])


# All available AUGUSTUS taxon ids
def augustus_taxon_ids() -> Tuple[List[str], List[str]]:
    tax_ids = """9606
7227
3702
6279
7159
7070
6183
5911
130081
4577
5811
6239
746128
162425
5062
33178
40559
5476
4929
5482
38033
5501
5346
49451
37769
40410
5207
4959
6035
33169
5518
5037
28985
29883
7757
5689
36914
148305
5141
5306
4924
64495
4932
4896
119072
6334
5270
4952
7425
4081
3055
400682
42068
4565
9031
7955
562
1280
3067""".split()
    identifiers =
