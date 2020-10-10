#ifndef REGION_H
#define REGION_H

#include <math.h>
#include <string>
#include <cstdlib>
#include <iostream>
#include "exceptions.h"

/* Region is exon/CDS location in a genome */


enum Strand {
    POS,
    NEG
};

struct Region {
    Region(const size_t& s, const size_t& e, const ssize_t& st, const size_t& off, const size_t& ct = 1)
        : is_terminal(false)
        , start(s)
        , end(e)
        , strand(st == 1 ? Strand::POS : Strand::NEG)
        , offset(off)
        , _count(ct) {
            if (start >= end) throw RangeException(
                std::string("Invalid region formed from " + s) + std::string(" to " + e)
            );
        }
    bool is_terminal;  // Mark as first/last exon
    size_t start;  // 0-based coordinate of start
    size_t end;  // 0-based coordinate of end
    Strand strand;  // Strand
    size_t offset; // Offset
    mutable size_t _count;  // Internal tracking for filter function
    std::string to_strand() const {
        switch (strand) {
            case Strand::POS:
                return "+";
            default:
                return "-";
        }
    }
    friend std::ostream& operator<<(std::ostream& o, const Region& rhs) {
        // o << "Start " << rhs.start << ", End: " << rhs.end;
        for (size_t i = 0; i < log10(rhs.end - rhs.start); i++) o << '|';
        return o;
    }
    // Define in terms of what is `not` part of the region
    bool operator<(const Region& rhs) const { return rhs.end < start; }
};
#endif