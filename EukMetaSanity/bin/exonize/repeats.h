#include <fstream>
#include <string>
#include <tuple>
#include <vector>

class RepeatsLocation {
    public:
        RepeatsLocation(std::istream*);
        ~RepeatsLocation();
        void read_next();  // Read in next record's repeat data
        bool is_in_repeat_region(int*);  // Check if genomic coord is a repeat
    private:
        std::vector<std::tuple<int, int>>* repeats;  // Repeats in contig
        std::istream *file;  // File pointer
};