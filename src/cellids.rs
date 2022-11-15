/// Cellids is a class that should implement the Rhapsody_cell_keys.r.txt functionality
/// Mainly to translate between R1 read and cell ID.
/// It should also fix some read errors in the cell ids. But that will come a little later

use std::collections::BTreeMap;
use kmers::naive_impl::Kmer;
use std::collections::HashSet;

//use std::thread;

//mod cellIDsError;
//use crate::cellids::cellIDsError::NnuclError;


/// GeneCount is an entry in the CellData and therefore does not need to collect it's own id
pub GeneCount{
    pub umi: HashSet<u64>
}

impl GeneCount{
    pub fn new() ->Self{
        umi = HashSet::new();
        Self{
            umi
        }
    }

    pub fn add( &mut self, umi:u64 ){
        self.umi.insert( umi );
    }

}

/// CellData here is a storage for the total UMIs. UMIs will be checked per cell
/// But I do not correct the UMIs here - even with sequencing errors 
/// we should get a relatively correct picture as these folow normal distribution.
pub struct CellData{
    pub gene: BTreeMap<u64, GeneCount>
}

impl CellData{
    pub fn new( ) -> Self{
        gene =  BTreeMap::new(); // to collect the sample counts
        Self{
            gene
        }
    }
    pub fn add(&mut self, geneid: u64, umi:u64 ){
        match self.gene.get( geneid) {
            Ok(mut gene) => gene.add( umi );
            Err(_) => {
                let gc = GeneCount::new();
                gc.insert( umi );
                gene.insert( geneid, gc );
        }
    }
    fn filleVec( &self, samples:Vec<u64>, ret:&Vec<u32>) {
        for i in 0..samples.len() {
            ret[i] = match self.gene.get( sample[id] ){
                Ok(id) => id,
                Err(_) => 0
            }
        }
    }
    pub fn to_str<'live>(&self, samples:Vec<u64> ) -> &'live str {
        let data = Vec<u32>;
        data = vec!(0, samples.len());
        fillVec( samples, data );
        let ret = "\t".join( data ); 
        return ret
    }
}

// and here the data
pub struct CellIds10x<'a>{    
    //kmer_len: usize,
    cells = BTreeMap<u64, CellData>
}


// here the functions
impl CellIds10x<'_>{

    pub fn new()-> Self {

        let mut cells = BTreeMap::<u64, CellData::new(10)>::new();
        let kmer_len = 5;
        // TACAGAACA

        Self {
            //kmer_len,
            cells
        }
    }

    pub fn add (&mut self, r1: &[u8], gene_id: u64  )-> Result< u32, &CellData>{
        let mut cell_id:CellData;

        let kmer =  &r1[0..16];
        for nuc in km{  
            if *nuc ==b'N'{
                //let tmp = std::str::from_utf8(km)?;
                return Err::<u32, &str>( "NnuclError");
                //Err::<i32, NnuclError<'_>>( NnuclError::new( &tmp ));
            }
        }
        let km = Kmer::from(kmer).into_u64();
        let umi = Kmer::from( &r1[16..r1.len()] ).into_u64();
        cell_id += match self.cells.get( km ){
            Some(c1) => {
                c1.add( gene_id, umi );
            },
            None => {
                let data = CellData::new();
                data.add( gene_id, umi );
                self.cells.insert( &km, data );
            }
        }
    }

}


// macro_rules! to_cellid {
//    ($r1: expr, $r2: expr) => {
//       //1-9,22-30,44-52
//       $r1.to_cellid($r1, $r2, [0,8], [21,29], [43,51]  )
//    };
// }

#[cfg(test)]
mod tests {
    #[test]
    fn getsamples() {
        let cells = crate::cell_ids::new();

        let mut primer = b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCT";
        let mut id:u32 = 1;
        let mut exp= ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
        match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
            Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
            Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
        };
        
        let mut exp2 = vec![b"GTCGCTATA", b"TACAGGATA", b"AAGCCTTCT"];
        assert_eq!( cells.to_sequence( exp ), exp2 );
        // 3, 3, 3
        primer = b"CTTCACATANNNNNNNNNNNNTGTGAAGAANNNNNNNNNNNNNCACAAGTAT";
        id = 3;
        exp = ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
        exp2 = vec![b"CTTCACATA", b"TGTGAAGAA", b"CACAAGTAT"];
        match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
            Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
            Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
        };
        //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
        assert_eq!( cells.to_sequence( exp ), exp2 );

        // and the last one
        primer = b"TGCGATCTANNNNNNNNNNNNCAACAACGGNNNNNNNNNNNNNCATAGGTCA";
        id = 96;
        exp = ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
        exp2 = vec![b"TGCGATCTA", b"CAACAACGG", b"CATAGGTCA"];
        assert_eq!( 884735+1 , exp);
        match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
            Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
            Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
        };
        //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
        assert_eq!( cells.to_sequence( exp ), exp2 );        
    }
}


