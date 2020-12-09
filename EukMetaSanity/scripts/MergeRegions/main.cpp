#include "Gene.h"
#include "Parser.h"
#include "exceptions.h"
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 5) {
        cout << "Usage: MergeRegions <fasta-file> <gff3-file> <gff3-file> <out-prefix>" << endl;
        return -1;
    }
    string fasta_file = argv[1];
    string gff3_file = argv[2];
    string evidence_file = argv[3];
    string output_prefix = argv[4];

    try {
        Parser gff(gff3_file, fasta_file, MergeType::EXTEND);
        gff.read_into(evidence_file, MergeType::CONFIRM);
        gff.write(output_prefix);
    } catch (runtime_error& e) {
        cout << e.what() << endl;
        return 1;
    }
    return 0;
}
