#include <tuple>
#include <utility>
#include <iostream>
#include <algorithm>
#include "parser.h"

inline void complement(char& cds) {
    switch (cds) {
        case 'T':
            cds = 'A';
            break;
        case 'A':
            cds = 'T';
            break;
        case 'C':
            cds = 'G';
            break;
        default:
            cds = 'C';
            break;
    }
}

void GffParser::write(const std::string& output_file) const {
    // Gff3 write
    std::fstream gff3_p;
    std::fstream fasta_p;
    std::fstream prot_p;
    std::vector<std::tuple<std::string, std::vector<char>>> fasta_list;
    gff3_p.open(output_file + ".nr.gff3", std::fstream::out);
    fasta_p.open(output_file + ".cds.fna", std::fstream::out);
    prot_p.open(output_file + ".faa", std::fstream::out);

    FastaMap fmap = load_fasta_file(_fasta_file);
    TranslationTable tmap = load_translation_table();
    write_gff3(gff3_p, fasta_list, fmap);
    // FASTA write
    write_fasta(fasta_p, fasta_list);
    // FASTA prot write
    write_protein(prot_p, fasta_list, tmap);
    // Close pointers
    gff3_p.close();
    fasta_p.close();
    prot_p.close();
}

// Read gff3 file into data
void GffParser::read_into(const std::string& gff3_file, const MergeType& mt) {
    // Open and read gff3 file
    std::fstream fp;
    fp.open(gff3_file, std::fstream::in);
    parse_gff3(fp, mt);
    fp.close();
}

void GffParser::parse_gff3(std::istream& fp, const MergeType& mt) {
    std::string line, id, source, _type;
    std::string vv, ee;
    size_t start, end, offset;
    char strand;
    size_t gene_count = 0;
    // Get initial line
    std::getline(fp, line);
    while (!fp.eof()) {
        if (line.size() > 0 && line.at(0) != '#') {
            std::stringstream ss;
            ss << line;
            ss >> id >> source >> _type >> start >> end >> ee >> strand;
            // Gene feature identified
            if (GffParser::at_gene_start(_type)) {
                // Get next line
                std::getline(fp, line);
                if (fp.eof()) return;
                std::stringstream ss;
                ss << line;
                ss >> id >> source >> _type >> start >> end >> ee >> strand >> offset;
                // Check if sequence id present, if not add
                Genome::iterator iter = find(id);
                // Sequence not present - insert contig object
                if (iter == this->end()) {
                    this->insert(std::make_pair(id, Contig()));
                    iter = find(id);
                }
                // Get reference to Contig object
                Contig& contig = iter->second;
                // Generate Record of GFF3 file
                Record record(start - 1, end, offset, GffParser::to_strand(strand));
                // Create gene object
                Gene gene(++gene_count, GffParser::to_strand(strand), mt);
                gene.insert(record);
                // Add regions to gene
                while (!GffParser::at_gene_start(_type)) {
                    std::getline(fp, line);
                    if (fp.eof()) break;
                    std::stringstream ss;
                    ss << line;
                    ss >> id >> source >> _type >> start >> end >> ee >> strand >> offset;
                    if (GffParser::at_gene_start(_type)) break;
                    if (_type == "CDS") {
                        // Line's record
                        Record record(start - 1, end, offset, GffParser::to_strand(strand));   
                        gene.insert(record);
                    }
                }
                // Find old gene
                Contig::iterator loc = contig.find(gene);
                // Combine together data
                if (loc != contig.end()) {
                    gene += *loc;
                    contig.erase(loc);
                }
                // Set in position
                contig.insert(gene);
            } else std::getline(fp, line);
        } else std::getline(fp, line);
    }
}

void GffParser::write_gff3(std::fstream& fp, FastaList& fasta_list, const FastaMap& fmap) const {
    fp << "##gff3-format" << std::endl;
    GffParser::const_iterator iter = this->cbegin();
    size_t i = 0;
    while (iter != this->cend()) {
        // Write sequence
        const std::set<Gene>& contigs = iter->second;
        std::set<Gene>::const_reverse_iterator genes = contigs.crbegin();
        FastaMap::const_iterator record = fmap.find(iter->first);
        fp << "# Begin region " << iter->first << std::endl;
        while (genes != contigs.crend()) {
            // Write gene info
            std::vector<std::vector<char>> list_of_cds_lists;
            const std::vector<Region> tmp = genes->get(1);
            if (tmp.size() == 0) {
                ++genes;
                continue;
            }
            char str = genes->strand() == 1 ? '+' : '-';
            fp << iter->first << '\t' << "EukMS" << '\t' << "gene" << '\t';
            fp << genes->range_start() + 1 << '\t' << genes->range_end() << '\t';
            fp << ".\t" << str << "\t.\tID=gene" << ++i << std::endl;
            std::vector<Region>::const_reverse_iterator reg = tmp.crbegin();
            size_t count = 0;
            while (reg != tmp.crend()) {
                size_t ss = 0;
                size_t ee = 0;
                if (reg->strand == 1) ss += reg->offset;
                else ee += reg->offset;
                // Write regions in gene
                fp << iter->first << '\t' << "EukMS" << '\t' << "CDS" << '\t';
                fp << reg->start + 1 << '\t' << reg->end << '\t';
                fp << "0" << '\t' << reg->to_strand() << '\t';
                fp << reg->offset << '\t';
                fp << "ID=gene" << i << "-cds" << ++count << ";Parent=gene" << i;
                fp << std::endl;
                std::vector<char> cds;
                for (size_t k = reg->start + ss; k < reg->end - ee; k++) {
                    cds.push_back(record->second.at(k));
                }
                list_of_cds_lists.push_back(cds);
                ++reg;
            }
            std::vector<char> gene_string;
            // Add to gene string
            for (size_t _c = 0; _c < list_of_cds_lists.size(); _c++) {
                std::vector<char>& cds = list_of_cds_lists.at(_c);
                for (size_t _i = 0; _i < cds.size(); _i++) {
                    gene_string.push_back(cds.at(_i));
                }
            }
            if (genes->strand() == -1) {
                std::reverse(gene_string.begin(), gene_string.end());
                for (size_t _i = 0; _i < gene_string.size(); _i++) 
                    complement(gene_string.at(_i));
            }
            std::string __id("gene");
            std::string num_string = std::to_string(i);
            for (size_t ns = 0; ns < num_string.size(); ns++) __id.push_back(num_string.at(ns));
            __id += " strand";
            __id += str;
            fasta_list.push_back(std::make_tuple(__id, gene_string));
            ++genes;
        }
        ++iter;
    }
}

void GffParser::write_fasta(std::fstream& fp, const FastaList& fasta_list) const {
    std::vector<std::tuple<std::string, std::vector<char>>>::const_iterator record = fasta_list.begin();
    while (record != fasta_list.cend()) {
        fp << '>' << std::get<0>(*record) << std::endl;
        const std::vector<char>& rec = std::get<1>(*record);
        size_t i;
        for (i = 0; i < rec.size(); i++) {
            if (i > 0 && i % 80 == 0) fp << std::endl;
            fp << rec.at(i);
        }
        fp << std::endl;
        ++record;
    }
}

void GffParser::write_protein(std::fstream& fp, const FastaList& fasta_list, const TranslationTable& tmap) const {
    std::vector<std::tuple<std::string, std::vector<char>>>::const_iterator record = fasta_list.begin();
    while (record != fasta_list.cend()) {
        fp << '>' << std::get<0>(*record) << std::endl;
        const std::vector<char>& rec = std::get<1>(*record);
        size_t i;
        for (i = 0; i < rec.size(); i += 3) {
            if (i + 2 >= rec.size()) break;
            if (i > 0 && i % 80 == 0) fp << std::endl;
            std::string val;
            val.push_back(toupper(rec.at(i)));
            val.push_back(toupper(rec.at(i + 1)));
            val.push_back(toupper(rec.at(i + 2)));
            fp << tmap.find(val)->second;
        }
        fp << std::endl;
        ++record;
    }
}

TranslationTable GffParser::load_translation_table() {
    return {
        {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'T'}, {"TTG", 'L'},
        {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
        {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '*'}, {"TAG", '*'},
        {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '*'}, {"TGG", 'W'},
        {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
        {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
        {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
        {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
        {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'},
        {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
        {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
        {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
        {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
        {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
        {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
        {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
    };
}

bool GffParser::at_gene_start(const std::string type) {
    return type == "gene" || type == "transcript";
}

ssize_t GffParser::to_strand(const char& strand) {
    return strand == '+' ? 1 : (strand == '-' ? -1 : throw StrandComparisonException("Invalid strand provided"));
}

FastaMap GffParser::load_fasta_file(const std::string& f) {
    std::fstream fs;
    std::unordered_map<std::string, std::vector<char>> out;    
    std::unordered_map<std::string, std::vector<char>>::iterator it;
    fs.open(f, std::fstream::in);
    if (!fs.is_open()) {
        throw FileException("Unable to open fasta file");
    }
    std::string line;
    std::getline(fs, line);
    while (!fs.fail()) {
        // FASTA record
        if (line.size() > 0 && line.at(0) == '>') {
            std::string id = line.substr(1, line.find(' '));
            // Add data to set
            out.insert(std::make_pair(id, std::vector<char>()));
            it = out.find(id);
            // Read in first line of data
            std::getline(fs, line);
            // Store sequence
            while (line.size() > 0 && line.at(0) != '>') {
                for (auto c: line) it->second.push_back(c);
                std::getline(fs, line);
            }
        } else std::getline(fs, line);
    }
    return out;
}