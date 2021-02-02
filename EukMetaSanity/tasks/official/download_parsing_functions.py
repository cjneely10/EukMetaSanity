"""
Module consists of functions to parse databases that are associated with official EukMS pipelines
"""


def odb_tax_parse(mmseqs_db_lookup_path: str, outfile: str):
    """ Generate NCBI taxonomy mappings from generated .lookup file

    :param mmseqs_db_lookup_path: Path to generated mmseqs database
    :param outfile: Output path for parsed orthodb mapping file
    """
    mmseqs_input_fp = open(mmseqs_db_lookup_path, "r")
    output_p = open(outfile, "w")
    for line in mmseqs_input_fp:
        line = line.split()
        output_p.write(line[1] + "\t" + line[1].split("_")[0] + "\n")
    output_p.close()
