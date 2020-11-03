#ifndef GENOME_TEST_H
#define GENOME_TEST_H

#include "../../lib/record.h"
#include "../../lib/region.h"
#include "ntest.h"
#include <sstream>
#include <string>
#include <vector>

class GeneTest : public NTest {
public:
    GeneTest() : NTest() { create_test_records(); }
    std::vector<Record> records;

private:
    void create_test_records();
};

void GeneTest::create_test_records() {
    records.push_back(Record(100, 600, 1, 1));
    records.push_back(Record(150, 650, 1, 1));
    records.push_back(Record(200, 450, 1, 1));
    records.push_back(Record(500, 1000, 1, 1));
    records.push_back(Record(2000, 10000, 1, 1));
}

#endif