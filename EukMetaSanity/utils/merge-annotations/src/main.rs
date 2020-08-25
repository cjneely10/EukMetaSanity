extern crate argparse;
extern crate bio;

use argparse::{ArgumentParser, Store};

fn main() {
    let mut fasta_file = "?".to_string();
    let mut gff3_file = "?".to_string();
    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
            "Write final merged results from EukMetaSanity"
        );
        ap.refer(&mut fasta_file)
            .add_option(&["-f", "--fasta-file"], Store, "Path to FASTA file").required();
        ap.refer(&mut gff3_file)
            .add_option(&["-g", "--gff3-file"], Store, "Path to GFF3 file").required();
        ap.parse_args_or_exit();
    }
    println!("{} {}", fasta_file, gff3_file);
}
