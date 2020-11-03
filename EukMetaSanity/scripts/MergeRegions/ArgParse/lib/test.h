#ifndef ARGPARSE_TEST_H
#define ARGPARSE_TEST_H

#include "argparse.h"
#include "ntest.h"
#include <iostream>
#include <string>
#include <vector>

using namespace ArgParse;

struct StringLengthValidator {
    bool operator()(std::string& val) { return val.size() > 4; }
    const char* ERR = "Too small of a string";
} length_validator;

class ArgParseTest : public NTest {
public:
    ArgParseTest(int& _argc, char** _argv) : NTest() {
        bool val = false;
        std::string data = "old";
        std::string data2 = "old2";
        size_t value = 1;
        {
            ArgParse::ArgParse ap(_argc, _argv, "ArgParse test program");

            ap.refer(val, {"--store-true"}, "make true", Type::PARSE_TRUE);
            ap.refer(value, {"--flagged-arg"}, "store flag", Type::FLAG);
            ap.refer(data, {"data-to-store"}, "store item", Type::DEFAULT, false);
            ap.refer(data2, {"data2-to-store"}, "store 2nd item", Type::DEFAULT, true, length_validator);
            ap.parse_or_exit();
        }
        std::cout << val << std::endl;
        std::cout << value << std::endl;
        std::cout << data << std::endl;
        std::cout << data2 << std::endl;
    }

private:
};

#endif
