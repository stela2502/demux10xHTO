use clap::Parser;
//use flate2::Compression;
//use flate2::write::GzEncoder;
// use flate2::read::GzDecoder;

use needletail::parse_fastx_file;
//use std::collections::HashSet;
//use std::io::BufWriter;
// use std::io::BufReader;
// use std::str;

mod cellids;
use crate::cellids::CellIds;
//use std::collections::BTreeMap;

mod geneids;
use crate::geneids::GeneIds;

//use std::path::PathBuf;
//use std::fs::File;
use std::fs;

//use std::io::prelude::*;
use std::time::SystemTime;

//use crate::glob;
//use glob::glob; // to search for files
//use regex::Regex;

// first, reproduce the appproach from
// https://github.com/jeremymsimon/SPLITseq/blob/main/Preprocess_SPLITseq_collapse_bcSharing.pl

/// Split a pair of BD rhapsody fastq files (R1 and R2) into sample specific fastq pairs
#[derive(Parser)]
#[clap(version = "0.1.0", author = "Stefan L. <stefan.lang@med.lu.se>, Rob P. <rob@cs.umd.edu>")]
struct Opts {
    /// the input R1 reads file
    #[clap(short, long)]
    reads: String,
    /// the input R2 samples file
    #[clap(short, long)]
    file: String,
    /// the barcodes table name<tab>bc
    #[clap(short, long)]
    bc: String,
    /// the outpath
    #[clap(short, long)]
    outpath: String,
}


#[derive(Debug, Deserialize)]
struct BCMapRecord {
    name: String,
    bc: String,
}

/// parse_bc_map will read in the bc map two column table and create the lookup samples 
/// class with this data
fn parse_bc_map(bc_map: &str, sub_len: usize ) -> &SampleIds {

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(true)
        .from_path(bc_map)
        .expect("cannot open barcode file");

    let sub_len = 9;
    let mut samples = SampleIds::new( sub_len );

    for result in rdr.deserialize() {
        // Notice that we need to provide a type hint for automatic
        // deserialization.
        let record: BCMapRecord = result.expect("could not deserialize barcode map record.");

        // now I need to get the seq into a [&u8] vector....
        // from https://gist.github.com/jimmychu0807/strinf-conversion.rs 
        // let byte1: Vec<u8> = src1.iter().map(|c| *c as u8).collect::<Vec<_>>();
        let seq: Vec<u8> =  record.bc.chars().map(|c| *c as u8).collect::<Vec<_>>();
        samples.add( seq, record.name );
    }

    samples 
}


fn main() {
    // parse the options

    let now = SystemTime::now();
    
    let opts: Opts = Opts::parse();

    match fs::create_dir_all(&opts.outpath){
        Ok(_) => (),
        Err(e) => panic!("I could not create the outpath: {}", e)
    };

    println!("starting to collect the cell ids per sample - can take a LONG time");

    let sub_len = 9;
    let mut samples:GeneIds = parse_bc_map( &opts.bc, sub_len );// = Vec::with_capacity(12);


    //  now we need to get a CellIDs object, too
    let mut cells = CellIds::new();

    let mut unknown = 0;

    {
        // need better error handling here too    
        // for now, we're assuming FASTQ and not FASTA.
        let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
        let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");

        while let Some(record2) = readefile.next() {
            if let Some(record1) = readereads.next() {
                let seqrec = record2.expect("invalid record");
                let seqrec1 = record1.expect("invalid record");

                // totally unusable sequence
                if seqrec1.seq().len() < 20 { // 16 + some UMI
                    unknown +=1;
                    continue;
                }
                //let seq = seqrec.seq().into_owned();
                for nuc in seqrec1.seq() {  
                    if *nuc ==b'N'{
                        unknown +=1;
                        continue;
                    }
                }
            
                match samples.get( &seqrec.seq() ){
                    Ok(id) => {
                        match cells.get( &seqrec1.seq() ){
                            Ok(gene_id) => {

                                cell2sample[id as usize].insert( val ); // will never insert one element twice. Great!
                            },
                            Err(_err) => {
                                //println!("{}",err);
                                continue
                            }, //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
                        };
                    },
                    Err(_err) => {
                        //println!("{}",err);
                        //seqrec1.write(&mut file1_ambig_out, None)?;
                        //seqrec.write(&mut file2_ambig_out, None)?;
                        //seqrec1.write(&mut outBuff1, None)?;
                        //seqrec.write(&mut outBuff2, None)?;
                        unknown +=1;
                    }
                }
            } else {
                println!("file 2 had reads remaining, but file 1 ran out of reads!");
            }
        }
        
        println!( "collected sample info:");
        /* for i in 0..cell2sample.len(){
            if cell2sample[i].len() > 0 {
        
                match samples.read_get( i as u32 )  {
                    Some(val) => println!( "    sample {}: {} reads and {} cells", i+1, val.total, cell2sample[i].len() ),
                    None => println!( "    sample {}: {} reads and {} cells", i+1, "na", cell2sample[i].len() ),
                };
        
                let fp1 = PathBuf::from(opts.reads.clone());
                //println!( "this is a the filename of the fastq file I'll use {}", fp1.file_name().unwrap().to_str().unwrap() );
                let file_path = PathBuf::from(&opts.outpath).join(
                    format!("{}.sample{}.ints.txt", fp1.file_name().unwrap().to_str().unwrap(), i+1)
                );
                let file = match File::create( file_path ){
                    Ok(file) => file,
                    Err(err) => {
                        eprintln!("Error: {:#?}", err);
                        std::process::exit(1)
                    }
                };
                let mut writer = BufWriter::new(&file);
                for int in &cell2sample[i]{
                    //file.write( )
                    //BigEndian::write_u32(file, int).unwrap();
                    writeln!(writer, "{}", int).unwrap();
                    //println!( "        with sample {}",int );
                }
            }
        }
        */
        println!(     "genomic reads: {} reads", unknown );
    }

}

