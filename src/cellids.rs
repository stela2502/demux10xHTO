/// Cellids is a class that should implement the Rhapsody_cell_keys.r.txt functionality
/// Mainly to translate between R1 read and cell ID.
/// It should also fix some read errors in the cell ids. But that will come a little later

use std::collections::BTreeMap;
use kmers::naive_impl::Kmer;

use crate::sampleids::SampleIds;

//use std::thread;

//mod cellIDsError;
//use crate::cellids::cellIDsError::NnuclError;


pub struct CellData{
    pub values = Vec<u32>
}


impl CellData{
    pub fn new(id:&[u8; 16], lenth: usize) -> Self{
        values = vec![ 0; length ];
        Self{
            id,
            values
        }
    }
    pub fn add1(&mut self, pos: usize ){
        self.values[pos] += 1;
    }
    pub fn to_string(&self, ) -> &str {

        ret:&mut str = std::str::from_utf8(self.id );
        for dat in self.values{
            ret = format!( "{} {}",ret, dat ); 
        }
        return ret;
    }
}

// and here the data
pub struct CellIds<'a>{    
    //kmer_len: usize,
    cells = BTreeMap<u64, CellData>
}


// here the functions
impl CellIds<'_>{

    pub fn new()-> Self {

        let mut cells = BTreeMap::<u64, CellData::new(10)>::new();
        let kmer_len = 5;
        // TACAGAACA

        Self {
            //kmer_len,
            cells
        }
    }

    pub fn to_cellid (&mut self, r1: &[u8], c1: Vec<usize>, c2: Vec<usize>, c3: Vec<usize>  )-> Result< u32, &CellData>{
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
        cell_id += match self.cells.get( &km ){
            Some(c1) => {
                //println!("to_cellid the c1 {}", c1 );
                c1
            },
            None => {
                self.cells.insert( &km, CellData::new(10) );
                Some( match self.cells.get( &km ) )
            }
        }
        return cell_id;
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


