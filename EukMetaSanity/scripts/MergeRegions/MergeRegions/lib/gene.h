#ifndef GENE_H
#define GENE_H

#include "record.h"
#include "region.h"
#include <iostream>
#include <set>
#include <stack>
#include <vector>

/* Gene contains a set of Regions of exons/CDS */

typedef std::vector<char> DnaString;
typedef std::set<Region> RegionSet;

enum MergeType { CONFIRM, EXTEND };

class Gene {
public:
    Gene(const size_t& i, const ssize_t& s, const MergeType& mt)
            : _id(i), _strand(s), merge_type(mt), min_pos(1000000000), max_pos(0) {}
    const size_t& id() const { return _id; }               // ID of gene
    const ssize_t& strand() const { return _strand; }      // Strand
    size_t size() const { return regions.size(); }         // Number of regions stored
    DnaString cds(const DnaString& dna) const;             // Get CDS string
    DnaString protein(const DnaString& dna) const;         // Get protein string
    size_t range() const;                                  // End_end - Start_start
    size_t range_start() const;                            // Start_start
    size_t range_end() const;                              // End_end
    std::vector<Region> get(const size_t& min_val) const;  // Get gene in ordered list
    void insert(const Record& record);                     // Insert record in sorted position
    void set_merge(const MergeType& mt);                   // Set merge method for future inserts

    // Find nearest region to record
    RegionSet::iterator find(const Record& record);
    // Aliases
    RegionSet::iterator begin();
    RegionSet::reverse_iterator rbegin();
    RegionSet::iterator end();
    RegionSet::reverse_iterator rend();

    // Stream operations
    friend std::ostream& operator<<(std::ostream& o, const Gene& rhs);
    // Comparison to order genes in contig - normal order based on position/strand
    bool operator<(const Gene& gene) { return gene.range_end() < range_start(); }
    // Friend comparisons
    friend bool operator<(const Gene& gene, const Gene& gene2) { return gene2.range_end() < gene.range_start(); }
    Gene& operator+=(const Gene& gene) {
        std::set<Region>::const_iterator it = gene.regions.cbegin();
        while (it != gene.regions.cend()) {
            regions.insert(*it);
            ++it;
        }
        return *this;
    }
    // Gene& operator()(const MergeType& min_val) { set_merge(min_val); }
private:
    size_t _id;            // Internal id storage
    ssize_t _strand;       // +/-
    RegionSet regions;     // Stored regions corresponding to genes
    MergeType merge_type;  // Compare or extend
    size_t min_pos;
    size_t max_pos;
};

#endif