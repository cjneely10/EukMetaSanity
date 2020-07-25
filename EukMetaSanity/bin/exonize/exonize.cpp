#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include "gene_structure.h"

/// Helper function declarations
bool get_fasta_file(int*, char*[], std::string*);
bool get_gff3_file(int*, char*[], std::string*);
bool print_help(int*, char*[]);
void get_output_file(int*, char*[], std::string*);


/// Main program logic
int main(int argc, char* argv[]) {
    // Variable declarations
    std::string fasta, gff3, output;
    // Parse input args
    if (print_help(&argc, argv)) return 0;
    if (!get_fasta_file(&argc, argv, &fasta)) return 1;
    if (!get_gff3_file(&argc, argv, &gff3)) return 2;
    get_output_file(&argc, argv, &output);
    // Confirm passed values are valid
    std::ifstream fasta_file(fasta);
    if (fasta_file.fail()) {
        std::cout << "FASTA file does not exist" << std::endl;
        return 1;
    }
    std::ifstream gff3_file(gff3);
    if (gff3_file.fail()) {
        std::cout << "GFF3 file does not exist" << std::endl;
        return 2;
    }
    std::ofstream outfile(output);
    if (outfile.fail()) {
        std::cout << "Unable to open output file" << std::endl;
        return 3;
    }
    return 0;
}

/// Helper functions
bool print_help(int* _argc, char* _argv[]) {
    for (int i = 1; i < *_argc; i++) {
        if (strcmp(_argv[i], "-h") == 0 || strcmp(_argv[i], "--help") == 0) {
            std::cout << "Usage: exonize <-f/--fasta-file file.mask.fna> " <<
                "<-g/--gff3-file file.gff3> [-o/--output-file file.gff3]" << std::endl;
            return true;
        }
    }
    return false;
}

bool get_fasta_file(int* _argc, char*_argv[], std::string* file_ptr) {
    for (int i = 1; i < *_argc; i++) {
        if (strcmp(_argv[i], "-f") == 0 || strcmp(_argv[i], "--fasta-file") == 0) {
            *file_ptr = _argv[i + 1];
            break;
        }
    }
    if (*file_ptr == "") {
        std::cout << "FASTA file not provided" << std::endl;
        return false;
    }
    return true;
}

bool get_gff3_file(int* _argc, char*_argv[], std::string* file_ptr) {
    for (int i = 1; i < *_argc; i++) {
        if (strcmp(_argv[i], "-g") == 0 || strcmp(_argv[i], "--gff3-file") == 0) {
            *file_ptr = _argv[i + 1];
            break;
        }
    }
    if (*file_ptr == "") {
        std::cout << "GFF3 file not provided" << std::endl;
        return false;
    }
    return true;
}

void get_output_file(int* _argc, char*_argv[], std::string* file_ptr) {
    for (int i = 1; i < *_argc; i++) {
        if (strcmp(_argv[i], "-o") == 0 || strcmp(_argv[i], "--output-file") == 0) {
            *file_ptr = _argv[i + 1];
            return;
        }
    }
    *file_ptr = "/dev/stdout";
}
