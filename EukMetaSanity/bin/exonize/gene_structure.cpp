#include "gene_structure.h"
#include "repeats.cpp"

GeneStructure::GeneStructure(std::istream* _fasta, std::istream* _gff3,
        std::ostream* _out) {
    fasta_file = _fasta;
    gff3_file = _gff3;
    out_file = _out;
    _region_start = 0;
    _region_end = 0;
    repeats = new RepeatsLocation(fasta_file);
}

GeneStructure::~GeneStructure() {
    delete repeats;
}

/// Private


/// Public

// Get left-most coordinate of current region
GenomeCoord GeneStructure::get_region_start() {
    return _region_start;
}

// Get right-most coordinate of current region
GenomeCoord GeneStructure::get_region_end() {
    return _region_end;
}

// Load in next CDS grouping
std::string GeneStructure::next() {
    return "";
}