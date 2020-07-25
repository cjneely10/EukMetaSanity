#include <fstream>
#include <tuple>
#include "repeats.h"

// Store if pointing to next part of exon structure,
// Or if to skip an intron to the next exon
enum Type {
    exon,
    intron,
};

// Position of gene structure boundary,
struct Node {
    GenomeCoord pos;  // Position of feature
    Type pos_type;  // Type of node stored
    Node** next;  // Pointer to next positions
    bool is_best;  // Track if part of optimal exon structure
};

#ifndef GENESTRUCTURE_H
#define GENESTRUCTURE_H
class GeneStructure {
    public:
        GeneStructure(std::istream*, std::istream*, std::ostream*);
        ~GeneStructure();
        GenomeCoord get_region_start();  // Current region's limit start
        GenomeCoord get_region_end();  // Current region limit end
        void load_next();  // Read in region chunk
        void write();  // Write and clear current region chunk
    private:
        RepeatsLocation* repeats;
        std::istream* fasta_file;
        std::istream* gff3_file;
        std::ostream* out_file;
        GenomeCoord _region_start;
        GenomeCoord _region_end;
        Node* gene_structure;  // Internal exon/intron structure
        Node* optimal_path;  // Store best path for ease in writing
        Node* new_node(GenomeCoord, bool, Type);
        void clear_node(Node*);  // Erase internal structure
        void add_node(GenomeCoord, Type);  // Add a node
        Node* find_prev(GenomeCoord*);  // Find node with closest value
        void insert_after(Node*);  // Place node after a found node
        void find_best_path(Node*, Node*);  // Load optimal path
};
#endif