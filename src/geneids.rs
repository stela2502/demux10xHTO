/// Geneids is a class that holds a set of kmers and allows to match a longer sequence to the set of kmers.
/// The id matching most kmers is returned.

use std::collections::BTreeMap;
//use std::collections::HashSet;

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

// and here the data
pub struct GeneIds{    
    pub kmers: BTreeMap<u64, Info>,
    kmer_size: usize,
    max_value: u32,
}

// here the functions
impl GeneIds{
    pub fn new(kmer_size: usize)-> Self {
        let kmers = BTreeMap::<u64, Info>::new();
        let max_value:u32 = 0;
        Self {
            kmers,
            kmer_size: kmer_size,
            max_value,
        }
    }

    pub fn add(&mut self, seq: &[u8], name: std::string::String ){
        for kmer in needletail::kmer::Kmers::new(seq, self.kmer_size as u8 ) {
            for nuc in kmer {
                if *nuc ==b'N'{
                    continue;
                }
            }
            //println!("Adding a gene id os length {} with seq {:?}", self.kmer_size, std::str::from_utf8(kmer) );
            // if  id == 1 { let s = str::from_utf8(kmer); println!( "this is the lib: {:?}",  s )};
            let km = Kmer::from(kmer).into_u64();
            let info = Info::new(km, name.clone() );
            self.kmers.insert(km, info);
        }
        self.max_value += 1;
    }


    pub fn get(&mut self, seq: &[u8] ) -> Result< &mut Info, &str>{
        
        // let min_value = 2;
        // let min_z = 1;
        // let mut max_value = 0;
        // let mut ret:u32 = 0;
        let kmers = needletail::kmer::Kmers::new(seq, self.kmer_size as u8 );
        let mut kmer_vec = Vec::<u64>::with_capacity(60);

        fill_kmer_vec(kmers, &mut kmer_vec);

        let id = 0;

        if kmer_vec.len() == 0 {
            eprintln!( "problematic sequence: {:?}", std::str::from_utf8( seq ) );
            return Err::<&mut Info, &str>( "Genes NoMatch")
        }  
        let km = kmer_vec[id];

        //println!("GeneIds::get - checking this sequence: {} and the at pos {}", km, id );
        let ret = match self.kmers.get_mut(&km){
            Some(c1) => c1, 
            None => return Err::<&mut Info, &str>( "Genes NoMatch"), 
        };
        Ok( ret )
    }

    pub fn to_ids( &self,  ret:&mut Vec<u64> )  {
        ret.clear();
        for (i, _obj) in &self.kmers {
            ret.push(*i );
        }
    }

    pub fn to_header( &self ) -> std::string::String {
        let mut ret: Vec<std::string::String>;
        ret = vec![" ".to_string(); self.kmers.len()];
        let mut id = 0;
        for (_, obj) in &self.kmers {
            ret[id] = obj.name.clone();
            id += 1;
        }
        return ret.join("\t")
    }

}




