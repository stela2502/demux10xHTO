use clap::Parser;
use serde::Deserialize;
//use flate2::Compression;
//use flate2::write::GzEncoder;
// use flate2::read::GzDecoder;

use needletail::parse_fastx_file;
use kmers::naive_impl::Kmer;
//use std::collections::HashSet;
use std::io::BufWriter;
// use std::io::BufReader;
// use std::str;

mod cellids;
use crate::cellids::CellIds10x;
//use std::collections::BTreeMap;

mod geneids;
use crate::geneids::GeneIds;

use std::path::PathBuf;
use std::fs::File;
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
    /// the input R2 genes file
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

/// parse_bc_map will read in the bc map two column table and create the lookup genes 
/// class with this data
fn parse_bc_map(bc_map: &str, sub_len: usize ) -> GeneIds {

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(true)
        .from_path(bc_map)
        .expect("cannot open barcode file");

    //let sub_len = 9;
    let mut genes = GeneIds::new( sub_len );

    for result in rdr.deserialize() {
        // Notice that we need to provide a type hint for automatic
        // deserialization.
        let record: BCMapRecord = result.expect("could not deserialize barcode map record.");

        // now I need to get the seq into a [&u8] vector....
        // from https://gist.github.com/jimmychu0807/strinf-conversion.rs 
        // let byte1: Vec<u8> = src1.iter().map(|c| *c as u8).collect::<Vec<_>>();
        //println!("what do I have here: {}", record.bc);
        let seq: Vec<u8> =  record.bc.chars().map(|c| c as u8).collect::<Vec<_>>();
        genes.add( &seq, record.name );
    }

    genes 
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
    let mut genes = parse_bc_map( &opts.bc, sub_len );// = Vec::with_capacity(12);


    //  now we need to get a CellIDs object, too
    let mut cells = CellIds10x::new(sub_len);

    let mut unknown:u32 = 0;
    let mut ok:u32 = 0;

    {
        // need better error handling here too    
        // for now, we're assuming FASTQ and not FASTA.
        let mut readereads = parse_fastx_file(&opts.reads).expect("valid path/file");
        let mut readefile = parse_fastx_file(&opts.file).expect("valid path/file");

        'main: while let Some(record2) = readefile.next() {
            if let Some(record1) = readereads.next() {
                let seqrec = record2.expect("invalid record");
                let seqrec1 = record1.expect("invalid record");

                // totally unusable sequence
                if seqrec1.seq().len() < 20 { // 16 + some UMI
                    unknown +=1;
                    continue 'main;
                }
                let seq = seqrec1.seq().into_owned(); // here the cells are coded
                for nuc in &seq {  
                    if *nuc ==b'N'{
                        unknown +=1;
                        continue 'main;
                    }
                }

                //let kmer = needletail::kmer::Kmers::new(&seq[16..seq.len()], (seq.len()-16) as u8 ).next();
                for kmer in needletail::kmer::Kmers::new(&seq[15..seq.len()], (seq.len()-16) as u8 ){ // should get me exactly 1
                    // println!("this will break - trying to add an entry for \nr1:{:?}\nr2:{:?}", 
                    //     std::str::from_utf8( &seqrec1.seq().into_owned()), 
                    //     std::str::from_utf8( &seqrec.seq().into_owned()) 
                    // );
                    let umi:u64  =  Kmer::from( kmer ).into_u64();
                    let cell:u64  = Kmer::from( &seq[0..16] ).into_u64();
                    let cell_name = match std::str::from_utf8( &seq[0..16] ){
                        Ok(t) => t,
                        Err(err) => panic!("the cellID could not be converted to string! {}",err)
                    };

                    match genes.get( &seqrec.seq().into_owned() ){
                        Ok(gene_id) => {
                            //println!("Gene has matched");
                            match cells.get( cell, cell_name.to_string() ){
                                Ok(cell_data) => {
                                    ok += 1;
                                    cell_data.add(gene_id.id, umi);
                                },
                                Err(_err) => {
                                    //eprintln!("But somehow we could not get the entry added!?! {} ", err);
                                    unknown +=1;
                                    continue 'main;
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
                }
            } else {
                println!("file 2 had reads remaining, but file 1 ran out of reads!");
            }
        }
        
        println!( "collected sample info:");
       
        let fp1 = PathBuf::from(opts.reads.clone());
        //println!( "this is a the filename of the fastq file I'll use {}", fp1.file_name().unwrap().to_str().unwrap() );
        let file_path = PathBuf::from(&opts.outpath).join(
            format!("Cell2Sample.{}.tsv", fp1.file_name().unwrap().to_str().unwrap() )
        );
        let file = match File::create( file_path ){
            Ok(file) => file,
            Err(err) => {
                eprintln!("Error: {:#?}", err);
                std::process::exit(1)
            }
        };
        let writer = BufWriter::new(&file);

        match cells.write ( writer, genes ) {
            Ok(_) => (),
            Err(err) => panic!("Error in the data write: {}", err)
        };
        
        println!("usable: {} reads", ok );
        println!("bad   : {} reads", unknown );
    }

    match now.elapsed() {
       Ok(elapsed) => {println!("demux10x finished in {} sec", elapsed.as_secs());},
       Err(e) => {println!("Error: {e:?}");}
    }


}

pub fn fill_kmer_vec<'a>( seq: needletail::kmer::Kmers<'a>, kmer_vec: &mut Vec<u64>) {
   kmer_vec.clear();
   let mut bad = 0;
   for km in seq {
        // I would like to add a try catch here - possibly the '?' works?
        // if this can not be converted it does not even make sense to keep this info
        for nuc in km{
            if *nuc ==b'N'{
                bad = 1;
            }
        }
        if bad == 0{
            // let s = str::from_utf8(km);
            // println!( "this is the lib: {:?}",  s );
            kmer_vec.push(Kmer::from(km).into_u64());
        }
   }
}
