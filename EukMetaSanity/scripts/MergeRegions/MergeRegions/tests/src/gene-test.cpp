#include "../lib/gene-test.h"
#include "../../lib/gene.h"

using namespace std;

int main(int argc, char* argv[]) {
    GeneTest test;
    for (auto region: test.records) {
        cout << region << endl;
    }
    test.make_test(
        "Merge overlapping regions",
        [&test](){
            ssize_t id = 0;
            Gene gene(++id, 1, MergeType::EXTEND);
            for (auto region: test.records) {
                gene.insert(region);
            }
            cout << gene << endl;
            size_t sz = gene.size();
            gene.set_merge(MergeType::CONFIRM);
            gene.insert(Record(2500, 12000, 1, 1));
            cout << gene << endl;
            size_t sz2 = gene.size();
            gene.set_merge(MergeType::EXTEND);
            gene.insert(Record(15000, 17000, 1, 1));
            cout << gene << endl;
            size_t sz3 = gene.size();
            return NTest::expect(sz == sz2 && sz3 - 1 == sz);
        }
    );
    test.make_test(
        "Find, get, etc.",
        [&test](){
            ssize_t id = 0;
            Gene gene(++id, 1, MergeType::EXTEND);
            for (auto region: test.records) {
                gene.insert(region);
            }
            cout << gene << endl;
            RegionSet::iterator it = gene.find(Record(1000, 12000, 1, 1));
            if (it != gene.end()) cout << *it << endl;
            return NTest::expect(it == gene.begin());
        }
    );
    test.make_test(
        "Compare gene locations",
        [&test](){
            ssize_t id = 0;
            Gene gene(++id, 1, MergeType::EXTEND);
            for (auto region: test.records) {
                gene.insert(region);
            }
            Gene gene2(++id, 1, MergeType::EXTEND);
            gene2.insert(Record(11000, 12000, 1, 1));
            return NTest::expect(gene2 < gene);
        }
    );
    test.make_test(
        "Filter regions",
        [&test](){
            ssize_t id = 0;
            Gene gene(++id, 1, MergeType::EXTEND);
            for (auto region: test.records) {
                gene.insert(region);
            }
            cout << gene << endl;
            gene.set_merge(MergeType::CONFIRM);
            gene.insert(Record(1000, 12000, 1, 1));
            gene.insert(Record(15000, 17000, 1, 1));
            cout << gene << endl;
            return NTest::expect(gene.size() == 3);
        }
    );
    test.make_test(
        "Test merging",
        [&test](){
            // Initial gene
            ssize_t id = 0;
            Gene gene(++id, 1, MergeType::EXTEND);
            // Exon 1
            gene.insert(Record(100, 500, 1, 1));
            // Exon 2
            gene.insert(Record(1000, 5000, 1, 1));
            cout << gene << endl;
            // Switch to confirmation
            gene.set_merge(MergeType::CONFIRM);
            // Single long mapped protein
            gene.insert(Record(100, 5000, 1, 1));
            auto result = gene.get(2);
            for (auto g: result) cout << g << endl;
            return true;
        }
    );
    test.test();
}