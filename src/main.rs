use clap::Parser;
use flate2::Compression;
use flate2::write::GzEncoder;
// use flate2::read::GzDecoder;

use needletail::parse_fastx_file;
use std::collections::HashSet;
use std::io::BufWriter;
// use std::io::BufReader;
// use std::str;

mod cellids;
use crate::cellids::CellIds;
use std::collections::BTreeMap;

mod sampleids;
use crate::sampleids::SampleIds;

use std::path::PathBuf;
use std::fs::File;
use std::fs;

use std::io::prelude::*;
use std::time::SystemTime;

//use crate::glob;
use glob::glob; // to search for files
use regex::Regex;

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

struct Ofiles {
    pub count: u32,
    //pub file1: GzEncoder<File>,
    //pub file2: GzEncoder<File>,
    pub buff1: BufWriter<GzEncoder<File>>,
    pub buff2: BufWriter<GzEncoder<File>>
}

impl Ofiles{
    pub fn new(id: usize, opts: &Opts )->Self {

        let fp1 = PathBuf::from(opts.reads.clone());
        let fp2 = PathBuf::from(opts.file.clone());
        let file1_path = PathBuf::from(&opts.outpath).join(format!("{}.{}.{}", opts.mode, id, fp1.file_name().unwrap().to_str().unwrap() ));
        let file2_path = PathBuf::from(&opts.outpath).join(format!("{}.{}.{}", opts.mode, id, fp2.file_name().unwrap().to_str().unwrap() ));
        
        // need better error handling here too
        // println!( "why does this file break? {}", file1_path.display() );
        let f1 = match File::create(file1_path){
            Ok(file) => file,
            Err(err) => panic!("The file {} cound not be created: {}", format!("{}.{}.{}", opts.mode, id, fp1.file_name().unwrap().to_str().unwrap() ), err)
        };
        let f2 = match File::create(file2_path){
            Ok(file) => file,
            Err(err) => panic!("The file {} cound not be created: {}", format!("{}.{}.{}", opts.mode, id, fp2.file_name().unwrap().to_str().unwrap() ), err)
        };
        
        let file1 = GzEncoder::new(f1, Compression::default());
        let file2 = GzEncoder::new(f2, Compression::default());

        let buff1 = BufWriter::new( file1 );
        let buff2 = BufWriter::new( file2 );

        let count:u32 = 0;
        Self{
            count,
            //file1,
            //file2,
            buff1,
            buff2
        }
    }
    pub fn close( &mut self ){

        match self.buff1.flush(){
            Ok(_) => (),
            Err(e) => eprintln!("Could not flush R1: {}",e),
        };
        match self.buff2.flush(){
            Ok(_) => (),
            Err(e) => eprintln!("Could not flush R2: {}",e),
        };
    }
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
        let seq: Vec<u8> =  =  result.bc.iter().map(|c| *c as u8).collect::<Vec<_>>();
        samples.add( seq, result.name );
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
    let mut samples:SampleIds = parse_bc_map( &opts.bc, sub_len );// = Vec::with_capacity(12);


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
                if seqrec1.seq().len() < 53 {
                    unknown +=1;
                    continue;
                }
                //let seq = seqrec.seq().into_owned();

                match samples.get( &seqrec.seq(), 9, 10 ){
                    Ok(id) => {
                        match cells.to_cellid( &seqrec1.seq(), vec![0,9], vec![21,30], vec![43,52]){
                            Ok(val) => {
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
        for i in 0..cell2sample.len(){
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
        println!(     "genomic reads: {} reads", unknown );
    }

}

