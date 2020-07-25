#include "repeats.h"

RepeatsLocation::RepeatsLocation(std::istream* infile) {
    file = infile;
    repeats = new std::vector<std::tuple<int, int>>;
}

RepeatsLocation::~RepeatsLocation() {
    delete repeats;
}

void RepeatsLocation::read_next() {

}