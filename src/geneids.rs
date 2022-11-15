/// Geneids is a class that holds a set of kmers and allows to match a longer sequence to the set of kmers.
/// The id matching most kmers is returned.

use std::collections::BTreeMap;
use std::collections::HashSet;

use kmers::naive_impl::Kmer;


//mod cellIDsError;
//use crate::cellids::cellIDsError::NnuclError;

#[derive(Debug)]
pub struct Info {
        pub cells: CellIds10x,
        pub name: &str
}


impl Info{
        pub fn new( name: &str )-> Self {
            name = name.clone();
            Self {
                name,
            }
        }
        pub fn add_umi(&mut self, seq: u64 ) {
            self.umi.insert( seq );
        }
        pub fn to_str<'life>( &mut self ) -> &'life str{
            let ret = format("{}\t{}",self.name,  self.umi.len() );
            return (ret)
        }
}

// and here the data
pub struct GeneIds{    
    kmers: BTreeMap<u64, u32>,
    kmer_size: usize,
    max_value: u32,
    pub info: BTreeMap<u32, Info>
}

// here the functions
impl GeneIds{
    pub fn new(kmer_size: usize)-> Self {
        let kmers = BTreeMap::<u64, u32>::new();
        let max_value:u32 = 0;
        let info = BTreeMap::<u32, Info>::new();
        Self {
            kmers,
            kmer_size: kmer_size,
            max_value,
            info
        }
    }

    pub fn add(&mut self, seq: &[u8], name: &str ){
        for kmer in needletail::kmer::Kmers::new(seq, self.kmer_size as u8 ) {
            // if  id == 1 { let s = str::from_utf8(kmer); println!( "this is the lib: {:?}",  s )};
            let km = Kmer::from(kmer).into_u64();
            let info = Info::new( name );
            self.kmers.insert(km, info);
        }
        self.max_value += 1;
    }


    pub fn get(&mut self, seq: &[u8] ) -> Result< u32, &str>{
        
        let min_value = 2;
        let min_z = 1;
        let mut max_value = 0;
        let mut ret:u32 = 0;
        let kmers = needletail::kmer::Kmers::new(seq, self.kmer_size as u8 );
        let mut kmer_vec = Vec::<u64>::with_capacity(60);

        fill_kmer_vec(kmers, &mut kmer_vec);

        let mut id = 0;

        let km = kmer_vec[id];
        //println!("GeneIds::get - checking this sequence: {} and the at pos {}", km, id );
        let ret = match self.kmers.get_mut(&km){
            Some(c1) => c1, 
            None => Err::<u32, &str>( "Genes NoMatch"), 
        }
        OK( ret )
    }
}


fn fill_kmer_vec<'a>(seq: needletail::kmer::Kmers<'a>, kmer_vec: &mut Vec<u64>) {
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

