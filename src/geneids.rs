/// Geneids is a class that holds a set of kmers and allows to match a longer sequence to the set of kmers.
/// The id matching most kmers is returned.

use std::collections::BTreeMap;
use std::collections::HashSet;

use kmers::naive_impl::Kmer;
use crate::fill_kmer_vec;

//mod cellIDsError;
//use crate::cellids::cellIDsError::NnuclError;

#[derive(Debug)]
pub struct Info {
//    pub cells: CellIds10x,
    pub id: u64,
    pub name: std::string::String
}


impl  Info{
        pub fn new( id: u64, name: std::string::String )-> Self {
            let loc_name = name.clone();
            Self {
                id,
                name: loc_name,
            }
        }
}

/// GeneIds harbors the antibody tags
/// these sequences have (to my knowmledge) 15 bp os length
/// but that can be different, too. 
/// Hence we store here 
/// kmers       : the search object
/// seq_len     : the length of the sequences (10x oversequences them)
/// kmer_size   : the length of the kmers
/// names       : a hashset for the gene names
/// bad_entries : a hash to save bad entries (repetetive ones)
pub struct GeneIds{    
    pub kmers: BTreeMap<u64, Info>,
    pub seq_len: usize,
    kmer_size: usize,
    pub names: BTreeMap<std::string::String, usize>,
    bad_entries: HashSet<u64>
}

// here the functions
impl GeneIds{
    pub fn new(kmer_size: usize )-> Self {
        let kmers = BTreeMap::<u64, Info>::new();
        let names = BTreeMap::<std::string::String, usize>::new();
        let bad_entries = HashSet::<u64>::new();
        let seq_len = 0;
        Self {
            kmers,
            seq_len,
            kmer_size: kmer_size,
            names,
            bad_entries
        }
    }

    pub fn add(&mut self, seq: &[u8], name: std::string::String ){
        
        if seq.len() > self.seq_len{
            self.seq_len = seq.len() 
        }
        for kmer in needletail::kmer::Kmers::new(seq, self.kmer_size as u8 ) {
            for nuc in kmer {
                if *nuc ==b'N'{
                    continue;
                }
            }
            //println!("Adding a gene id os length {} with seq {:?}", self.kmer_size, std::str::from_utf8(kmer) );
            // if  id == 1 { let s = str::from_utf8(kmer); println!( "this is the lib: {:?}",  s )};
            let km = Kmer::from(kmer).into_u64();
            if self.kmers.contains_key ( &km ){
                self.bad_entries.insert( km.clone() );
                self.kmers.remove( &km );
            }else {
                let info = Info::new(km, name.clone() );
                self.kmers.insert(km, info);
                self.names.insert( name.clone(), self.names.len() );
            }
        }
    }

    pub fn get(&mut self, seq: &[u8] ) -> Option< u64 >{
        
        // let min_value = 2;
        // let min_z = 1;
        // let mut max_value = 0;
        // let mut ret:u32 = 0;
        let kmers = needletail::kmer::Kmers::new(seq, self.kmer_size as u8 );
        let mut kmer_vec = Vec::<u64>::with_capacity(60);

        fill_kmer_vec(kmers, &mut kmer_vec);

        let mut ret:Option<u64> = None;

        if kmer_vec.len() == 0 {
            eprintln!( "Seq length: {};  -> problematic sequence: {:?}", self.kmer_size, std::str::from_utf8( seq ) );
            return ret
        }  

        for km in kmer_vec{
            ret = match self.kmers.get(&km){
                Some(c1) => Some(c1.id), 
                None => continue
            };
        }
        return ret
    }

    // pub fn to_ids( &self,  ret:&mut Vec<Info> )  {
    //     ret.clear();
    //     for (_i, obj) in &self.kmers {
    //         ret.push(*obj );
    //     }
    // }

    pub fn to_header( &self ) -> std::string::String {
        let mut ret= Vec::<std::string::String>::with_capacity( self.names.len() +2 );
        //println!( "I get try to push into a {} sized vector", self.names.len());
        for (obj, _id) in &self.names {
            //println!( "Pushing {} -> {}", obj, *id-1);
            ret.push( format!( "{}", obj)) ;
        }
        ret.push("Faction total".to_string());
        ret.push("Most likely name".to_string());
        return "CellID\t".to_owned()+&ret.join("\t")
    }

}



