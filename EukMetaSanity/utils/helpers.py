import os
from typing import Dict

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
def augustus_taxon_ids() -> Dict[str, str]:
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
    identifiers = """human
fly
arabidopsis
brugia
aedes
tribolium
schistosoma
tetrahymena
galdieria
maize
toxoplasma
caenorhabditis
aspergillus_fumigatus
aspergillus_nidulans
aspergillus_oryzae
aspergillus_terreus
botrytis_cinerea
candida_albicans
candida_guilliermondii
candida_tropicalis
chaetomium_globosum
coccidioides_immitis
coprinus
coyote_tobacco
cryptococcus_neoformans_gattii
cryptococcus_neoformans_neoformans_B
(cryptococcus)
debaryomyces_hansenii
encephalitozoon_cuniculi_GB
eremothecium_gossypii
fusarium_graminearum
histoplasma_capsulatum
kluyveromyces_lactis
laccaria_bicolor
lamprey
leishmania_tarentolae
lodderomyces_elongisporus
magnaporthe_grisea
neurospora_crassa
phanerochaete_chrysosporium
pichia_stipitis
rhizopus_oryzae
saccharomyces_cerevisiae_S288C
schizosaccharomyces_pombe
thermoanaerobacter_tengcongensis
trichinella
ustilago_maydis
yarrowia_lipolytica
nasonia
tomato
chlamydomonas
amphimedon
pneumocystis
wheat
chicken
zebrafish
E_coli_K12
s_aureus
volvox
""".split()
    return {key: val for key, val in zip(tax_ids, identifiers)}
