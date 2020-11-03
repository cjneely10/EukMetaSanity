#ifndef GFFPARSER_H
#define GFFPARSER_H

#include "exceptions.h"
#include "gene.h"
#include "record.h"
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

typedef std::unordered_map<std::string, DnaString> FastaMap;
typedef std::vector<std::tuple<std::string, std::vector<char>>> FastaList;
typedef std::unordered_map<std::string, char> TranslationTable;
typedef std::set<Gene> Contig;
typedef std::unordered_map<std::string, Contig> Genome;

class GffParser : public Genome {
public:
    // Constructor generate translation table and fasta file map
    GffParser(const std::string& gff3_file, const std::string& fasta_file, const MergeType& mt = MergeType::CONFIRM)
            : Genome(), _fasta_file(fasta_file) {
        read_into(gff3_file, mt);
    }
    void write(const std::string& output_file) const;
    // Read gff3 file into Genome structure
    void read_into(const std::string& gff3_file, const MergeType& mt);

protected:
    // Parse GFF3 file
    void parse_gff3(std::istream& fp, const MergeType& mt);
    // Write inner results to file
    void write_gff3(std::fstream& fp, FastaList& fasta_list, const FastaMap& fmap) const;
    // Convert regions to fasta and write
    void write_fasta(std::fstream& fp, const FastaList& fasta_list) const;
    void write_protein(std::fstream& fp, const FastaList& fasta_list, const TranslationTable& tmap) const;
    static FastaMap load_fasta_file(const std::string& f);

private:
    std::string _fasta_file;
    // Check if beginning of gene region
    static bool at_gene_start(const std::string type);
    // Convert char strand to internal strand identifier
    static ssize_t to_strand(const char& strand);
    static TranslationTable load_translation_table();
};

#endif