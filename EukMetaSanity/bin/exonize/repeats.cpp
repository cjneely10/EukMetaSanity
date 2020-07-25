#include "repeats.h"

RepeatsLocation::RepeatsLocation(std::istream* infile) {
    file = infile;
    repeats = NULL;
}

RepeatsLocation::~RepeatsLocation() {
    clear_repeats();
}

/// Private

// Clear repeats data struct
void RepeatsLocation::clear_repeats() {
    if (repeats) delete repeats;
}


/// Public

// Read in contig data
void RepeatsLocation::read_next() {
    clear_repeats();
    repeats = new std::set<GenomeCoord>;
    std::string line;
    getline(*file, line);
    // Read in id
    id = line.substr(
        1,  // Start after the FASTA delimiter
        line.find(' ') - 1  // Read length
    );
    char val;
    *file >> val;
    GenomeCoord num_endlines = 0;
    // Read up until the next record
    GenomeCoord pos = 1;
    while (val != '>') {
        // Read repeats structure into internal data struct
        if (val == 'N') repeats->insert(pos - num_endlines);
        if (val == '\n') num_endlines += 1;
        else pos += 1;
        *file >> val;
    }
}

// Search for coord in set of coordinates for contig
bool RepeatsLocation::is_in_repeat_region(GenomeCoord* pos) {
    if (repeats->find(*pos) != repeats->end()) return true;
    return false;
}

// Get stored id
std::string RepeatsLocation::get_id() {
    return id;
}