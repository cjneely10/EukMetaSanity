#ifndef RECORD_H
#define RECORD_H

#include "exceptions.h"
#include <string>
#include <iostream>

/* Record is line of GFF3 file */

struct Record {
    Record(const size_t& st, const size_t& en, const size_t& off, const ssize_t& str)
            : start(st), end(en), offset(off), strand(str) {
        if (start >= end)
            throw RangeException(std::string("Invalid region formed from " + st) + std::string(" to " + en));
    }
    size_t start;    // Start coordinate 0-based, inclusive
    size_t end;      // End coordinate, non-inclusive
    size_t offset;   // Offset for ORF
    ssize_t strand;  // Strand
    friend std::ostream& operator<<(std::ostream& o, const Record& record) {
        o << "Start: " << record.start + 1;
        o << " End: " << record.end;
        o << " Offset: " << record.offset;
        o << " Strand: " << record.strand;
        return o;
    }
};

#endif