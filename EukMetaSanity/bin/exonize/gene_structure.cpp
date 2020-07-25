#include "gene_structure.h"
#include "repeats.cpp"

GeneStructure::GeneStructure(std::istream* _fasta, std::istream* _gff3,
        std::ostream* _out) {
    fasta_file = _fasta;
    gff3_file = _gff3;
    out_file = _out;
    gene_structure = NULL;
    optimal_path = NULL;
    _region_start = 0;
    _region_end = 0;
    repeats = new RepeatsLocation(_fasta);
    repeats->read_next();
}

GeneStructure::~GeneStructure() {
    delete repeats;
}

/// Private

// Create node using 1-based coordinate
Node* GeneStructure::new_node(GenomeCoord pos, bool is_best, Type val_type) {
    Node* node = new Node;
    node->is_best = is_best;
    node->pos = pos;
    node->next = new Node*;
    node->pos_type = val_type;
    return node;
}


/// Public

// Get left-most coordinate of current region
GenomeCoord GeneStructure::get_region_start() {
    return _region_start;
}

// Get right-most coordinate of current region
GenomeCoord GeneStructure::get_region_end() {
    return _region_end;
}