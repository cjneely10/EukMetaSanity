#include "gene.h"

#include <algorithm>

size_t Gene::range() const {
    if (regions.size() == 0)
        throw RangeException("Empty Gene");
    // Entire region spanned by gene
    return range_end() - range_start();
}

size_t Gene::range_start() const {
    if (regions.size() == 0)
        throw RangeException("Empty Gene");
    return min_pos;
}

size_t Gene::range_end() const {
    if (regions.size() == 0)
        throw RangeException("Empty Gene");
    return max_pos;
}

std::vector<Region> Gene::get(const size_t& min_val) const {
    // Push values to stack as means of sorting in forward order
    std::vector<Region> tmp;
    RegionSet::const_iterator it = regions.cbegin();
    while (it != regions.cend()) {
        tmp.push_back(*it);
        ++it;
    }
    std::sort(tmp.begin(), tmp.end(), std::less<Region>());
    std::reverse(tmp.begin(), tmp.end());
    tmp.begin()->is_terminal = true;
    tmp.rbegin()->is_terminal = true;

    std::vector<Region> out;
    std::vector<Region>::const_iterator tt = tmp.cbegin();
    while (tt != tmp.cend()) {
        if (tt->is_terminal || tt->_count >= min_val)
            out.push_back(*tt);
        ++tt;
    }
    return out;
}

void Gene::insert(const Record& record) {
    // Create Region object
    Region region = Region(record.start, record.end, record.strand, record.offset);
    // Search for existing record
    RegionSet::iterator it = regions.find(region);
    // Insert if not present
    if (it == regions.end())
        regions.insert(region);
    // Update with new bounds
    else {
        if (merge_type == MergeType::EXTEND) {
            // Extend forward within given direction
            if (region.offset == it->offset) {
                region.end = record.end < it->end ? it->end : record.end;
                region.start = record.start < it->start ? record.start : it->start;
            }
            // Update region and re-insert
            region._count = it->_count + 1;
            regions.erase(it);
            regions.insert(region);
            // Track start/end boundaries
            if (region.start < min_pos)
                min_pos = it->start;
            if (region.end > max_pos)
                max_pos = it->end;
        } else {
            // Update all overlapping regions starting at find spot
            RegionSet::iterator matched_regions = it;
            while (matched_regions != regions.end()) {
                // For ovelapping regions
                if (std::max(matched_regions->start, region.start) <= std::min(region.end, matched_regions->end)) {
                    matched_regions->_count += 1;
                    ++matched_regions;
                } else
                    break;
            }
        }
        return;
    }
    if (region.start < min_pos)
        min_pos = region.start;
    if (region.end > max_pos)
        max_pos = region.end;
}

void Gene::set_merge(const MergeType& mt) { merge_type = mt; }

RegionSet::iterator Gene::find(const Record& record) {
    Region region = Region(record.start, record.end, record.strand, record.offset);
    return regions.find(region);
}

RegionSet::iterator Gene::begin() { return regions.begin(); }

RegionSet::reverse_iterator Gene::rbegin() { return regions.rbegin(); }

RegionSet::iterator Gene::end() { return regions.end(); }

RegionSet::reverse_iterator Gene::rend() { return regions.rend(); }

Gene& Gene::operator+=(const Gene& gene) {
    std::set<Region>::const_iterator it = gene.regions.cbegin();
    while (it != gene.regions.cend()) {
        regions.insert(*it);
        ++it;
    }
    return *this;
}

//// Friend functions

std::ostream& operator<<(std::ostream& o, const Gene& rhs) {
    std::vector<Region> tmp = rhs.get(1);
    std::vector<Region>::iterator it = tmp.begin();
    o << "gene" << rhs.id() << " " << rhs.strand() << std::endl;
    size_t end = it->end;
    while (it != tmp.end()) {
        // o << *it << " ";
        o << *it;
        ++it;
        if (it != tmp.end()) {
            size_t bigger = it->start > end ? it->start : end;
            size_t smaller = it->start < end ? it->start : end;
            for (size_t i = 0; i < log10(bigger - smaller); i++)
                o << '-';
        }
        end = it->end;
    }
    return o;
}

std::ostream& operator<<(std::ostream& o, const std::vector<Region>& rhs) {


    return o;
}