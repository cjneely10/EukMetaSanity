#include "../lib/parser-test.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    MergeParserTest test(argc, argv);
    test.make_test("FASTA file load of size 369", [&test]() {
        auto mm = ParserTester::load_fasta_file(test.fasta_file);
        return NTest::expect(mm.size() == 369);
    });
    test.make_test("Load GFF3 file without failing", [&test]() {
        GffParser parser(test.gff3_file, test.fasta_file, MergeType::EXTEND);
        size_t contigs = parser.size();
        size_t genes_found = 0;
        cout << "Num contigs: " << contigs << endl;
        GffParser::iterator iter = parser.begin();
        while (iter != parser.end()) {
            cout << "Contig id: " << iter->first << endl;
            genes_found += iter->second.size();
            cout << "Num genes: " << iter->second.size() << endl;
            for (auto gene : iter->second)
                cout << gene << endl;
            ++iter;
        }
        return NTest::expect(genes_found == 16);
    });
    test.test();
}