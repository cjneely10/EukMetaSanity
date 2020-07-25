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
    long long int pos;  // Position of feature
    std::tuple<Type, Node*> next;  // Pointer to next position and type
    bool is_best;  // Track if part of optimal exon structure
};

// Path for traceback
struct Path {
    Node* cur;
    Node* next;
};

class GeneStructure {
    public:
        GeneStructure(std::istream*, std::istream*, std::ostream*);
        ~GeneStructure();
        long long int region_start;  // Current region's limit start
        long long int region_end;  // Current region limit end
        void load_next();  // Read in region chunk
        void write();  // Write and clear current region chunk
    private:
        Node* gene_structure;  // Internal exon/intron structure
        Path* optimal_path;  // Store best path for ease in writing
        void clear_gene(Node*);  // Erase internal structure
        void add_node(long long, Type);  // Add a node
        Node* find_prev(long long*);  // Find node with closest value
        void insert_after(Node*);  // Place node after a found node
        void find_best_path(Node*, Path*);  // Load optimal path
};
