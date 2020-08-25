extern crate argparse;
extern crate bio;
extern crate bio_types;

// IO operations
use std::fs::File;
use std::path::{Path, PathBuf};
use argparse::{ArgumentParser, Store};

// Bio operations
use bio::io::{gff, fasta};
use bio::data_structures::annot_map::AnnotMap;
use bio_types::annot::contig::Contig;
use bio_types::strand::ReqStrand;

fn main() {
    let mut fasta_file = "?".to_string();
    let mut gff3_file = "?".to_string();
    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
            "Write final merged results from EukMetaSanity"
        );
        ap.refer(&mut fasta_file)
            .add_option(&["-f", "--fasta-file"], Store, "Path to FASTA file")
            .required();
        ap.refer(&mut gff3_file)
            .add_option(&["-g", "--gff3-file"], Store, "Path to GFF3 file")
            .required();
        ap.parse_args_or_exit();
    }
    let fasta_file = Path::new(&fasta_file);
    let gff3_file = Path::new(&gff3_file);
    let (mut f_iter, mut g_iter, mut g_writer) = initialize_readers_writer(
        &fasta_file, &gff3_file
    );
    println!(
        "{:?}",
        g_iter
        .records()
        .filter_map(|x| x.ok())
        .take(8).collect::<Vec<gff::Record>>()
    );
}

/// Create reader/writer objects with valid paths
fn initialize_readers_writer(fasta_file: &Path, gff3_file: &Path) 
            -> (fasta::Reader<File>, gff::Reader<File>, gff::Writer<File>) {
    let freader = fasta::Reader::new(
        File::open(fasta_file).expect("Unable to open FASTA file!")
    );
    let greader = gff::Reader::new(
        File::open(gff3_file).expect("Unable to open GFF3 file!"), 
        gff::GffType::GFF3
    );

    let mut buf = PathBuf::from(gff3_file.parent().unwrap());
    buf.push(gff3_file.file_stem().unwrap());
    buf.set_extension("nr.gff3");

    let gwriter = gff::Writer::new(
        File::create(buf).expect("Unable to open output file!"), 
        gff::GffType::GFF3
    );
    (freader, greader, gwriter)
}

/// Read data from GFF3 into AnnotMap
fn load_to_map() {

}