#include <fstream>
#include <string>
#include <tuple>
#include <set>
#include <iostream>

typedef unsigned long long int GenomeCoord;

#ifndef REPEATS_H
#define REPEATS_H
class RepeatsLocation {
    public:
        RepeatsLocation(std::istream*);
        ~RepeatsLocation();
        void read_next();  // Read in next record's repeat data
        bool is_in_repeat_region(GenomeCoord*);  // Check if genomic coord is a repeat
        std::string id;
    private:
        std::set<GenomeCoord>* repeats;  // Repeats in contig
        std::istream *file;  // File pointer
};
#endif