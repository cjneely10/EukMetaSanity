#include "gene.h"
#include <algorithm>

size_t Gene::range() const {
    if (regions.size() == 0) throw RangeException("Empty Gene");
    // Entire region spanned by gene
    return range_end() - range_start();
}

size_t Gene::range_start() const {
    if (regions.size() == 0) throw RangeException("Empty Gene");
    return min_pos;
}

size_t Gene::range_end() const {
    if (regions.size() == 0) throw RangeException("Empty Gene");
    return max_pos;
}

std::vector<Region> Gene::get(const size_t& min_val) const {
    // Push values to stack as means of sorting in forward order
    std::vector<Region> tmp;
    RegionSet::const_iterator it = regions.cbegin();
    while (it != regions.cend()) {
        if (it->_count >= min_val) tmp.push_back(*it);
        ++it;
    }
    std::sort(tmp.begin(), tmp.end());
    return tmp;
}

void Gene::insert(const Record& record) {
    if (_strand != record.strand) return;
    // Create Region object
    Region region = Region(record.start, record.end, record.strand, record.offset);
    // Search for existing record
    RegionSet::iterator it = regions.find(region);
    // Insert if not present
    if (it == regions.end()) regions.insert(region);
    // Update with new bounds
    else {
        if (merge_type == MergeType::EXTEND) {
            // Extend forward within given direction
            switch (_strand) {
                case 1:
                    region.end = record.end < it->end ? it->end : record.end;
                    // region.start = region.start > it->start ? it->start : record.start;
                    break;
                case -1:
                    // region.end = record.end > it->end ? it->end : record.end;
                    region.start = region.start < it->start ? it->start : record.start;
                    break;
                default:
                    break;
            }
            region._count = it->_count + 1;
            regions.erase(it);
            regions.insert(region);
            if (region.start < min_pos) min_pos = region.start;
            if (region.end > max_pos) max_pos = region.end;
            return;
        } else {
            it->_count = it->_count + 1;
            (--it)->_count += 1;
        }
    }
    if (region.start < min_pos) min_pos = region.start;
    if (region.end > max_pos) max_pos = region.end;
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

//// Friend functions

std::ostream& operator<<(std::ostream& o, const Gene& rhs) {
    std::vector<Region> tmp = rhs.get(1);
    std::vector<Region>::iterator it = tmp.begin();
    size_t end = it->end;
    while (it != tmp.end()) {
        // o << *it << " ";
        o << *it;
        ++it;
        size_t bigger = it->start > end ? it->start : end;
        size_t smaller = it->start < end ? it->start : end;
        if (it != tmp.end())
            for (size_t i = 0; i < log10(bigger - smaller); i++) o << '-';
        end = it->end;
    }
    return o;
}