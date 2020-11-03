#include "argparse.h"
#include "exceptions.h"
#include "gene.h"
#include "parser.h"
#include <fstream>
#include <string>

using namespace std;

/// Validators
struct FileExists {
    bool operator()(string& ref) {
        fstream fs;
        fs.open(ref, fstream::in);
        bool opened = fs.is_open();
        if (opened)
            fs.close();
        return opened;
    }
    const char* ERR = "Unable to locate file";
} file_exists;

/// Main logic
int main(int argc, char* argv[]) {
    string fasta_file;
    string gff3_file;
    string evidence_file;
    string output_prefix = "out";
    {
        ArgParse::ArgParse ap(argc, argv, "Finalize output of EukMetaSanity");
        ap.refer(output_prefix, {"-o", "--output-prefix"}, "Output prefix for files", ArgParse::Type::FLAG, false);
        ap.refer(gff3_file, {"gff3-file"}, "First GFF3 file path", ArgParse::Type::DEFAULT, true, file_exists);
        ap.refer(evidence_file, {"gff3-file2"}, "Evidence GFF3 file path", ArgParse::Type::DEFAULT, true, file_exists);
        ap.refer(fasta_file, {"fasta-file"}, "FASTA file path", ArgParse::Type::DEFAULT, true, file_exists);
        ap.parse_or_exit();
    }
    try {
        GffParser gff(gff3_file, fasta_file, MergeType::EXTEND);
        gff.read_into(evidence_file, MergeType::CONFIRM);
        gff.write(output_prefix);
    } catch (runtime_error& e) {
        cout << e.what() << endl;
        return 1;
    }
    return 0;
}