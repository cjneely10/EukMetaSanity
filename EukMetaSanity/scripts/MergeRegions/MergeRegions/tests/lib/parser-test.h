#ifndef MERGE_PARSER_TEST_H
#define MERGE_PARSER_TEST_H

#include <string>
#include "ntest.h"
#include "../../lib/parser.h"

class ParserTester: protected GffParser {
    public:
        ParserTester(std::string& f, std::string& g) : GffParser(g, f) { }
        static std::unordered_map<std::string, std::vector<char>> load_fasta_file(std::string f) {
            return GffParser::load_fasta_file(f);
        }
};

class MergeParserTest: public NTest {
    public:
        MergeParserTest(int argc, char** argv): NTest() {
            if (argc != 3) return;
            fasta_file = argv[1];
            gff3_file = argv[2];
        }
        std::string fasta_file;
        std::string gff3_file;
    private:
};

#endif