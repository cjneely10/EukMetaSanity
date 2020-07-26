#include <fstream>
#include <tuple>
#include <string>
#include "repeats.h"

#ifndef GENESTRUCTURE_H
#define GENESTRUCTURE_H
class GeneStructure {
    public:
        GeneStructure(std::istream*, std::istream*, std::ostream*);
        ~GeneStructure();
        GenomeCoord get_region_start();  // Current region's limit start
        GenomeCoord get_region_end();  // Current region limit end
        std::string next();  // Return corrected GFF3 line
    private:
        RepeatsLocation* repeats;
        std::istream* fasta_file;
        std::istream* gff3_file;
        std::ostream* out_file;
        GenomeCoord _region_start;
        GenomeCoord _region_end;
};
#endif