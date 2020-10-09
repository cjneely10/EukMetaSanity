# CPP Library for Bioinformatics Research

## Installation
Requires c++17.

```
git clone https://github.com/cjneely10/lib.git
cd lib
```

Adjust the `INSPATH` variable to match your installation location

```
CC = g++
INSPATH = /home/Chris/lib
```

Install using `make`. Run tests if desired

```
make
make tests
make run-tests
```

All header files will now be in `$(INSPATH)`

## NTest

Test your code by subclassing NTest. Example code:

```
//// test.h
#include "ntest.h"
#include <string>
#include <vector>
#include <map>

#ifndef TESTER_H
#define TESTER_H
class NTestHash: public NTest {
    public:
        NTestHash(size_t size) : NTest() {
            for (size_t i = 0; i < size; i++) {
                ids.push_back(random_string(50));
            }
            make_test(
                "Create", 
                [this]() {
                    std::cerr << "Creating mm" << std::endl;
                    double mm_time = NTest::experiment(
                        [this](){
                            std::map<std::string, size_t> mm;
                            for (size_t i = 0; i < this->size(); i++) {
                                mm[this->get(i)] = i;
                            }
                        }
                    );
                    std::cerr << mm_time << std::endl;
                    return NTest::expect(2 > mm_time);
                }
            );
        }
    private:
        const std::string get(const size_t& i) const { return ids[i]; }
        const size_t size() const { return ids.size(); }
        std::vector<std::string> ids;
};

#endif



//// main.cpp
#include "test.h"


using namespace std;

int main() {
    NTestHash test(1000000);
    test.test();
}
```

## HashTable
Achieve lookup-time that is, on average, twice as fast as `std::set`. Example code:

```
#include "hashtable.h"


using namespace std;

int main() {
    HashTable<std::string, size_t> ht;
    ht.insert("a", 1);
    ht["a"] == 1;  // Does nothing, but is true
}

```

