/// Cellids is a class that should implement the Rhapsody_cell_keys.r.txt functionality
/// Mainly to translate between R1 read and cell ID.
/// It should also fix some read errors in the cell ids. But that will come a little later

use std::collections::BTreeMap;
//use kmers::naive_impl::Kmer;
use std::collections::HashSet;

//use std::thread;
use crate::geneids::GeneIds;
//use crate::geneids::Info;

//mod cellIDsError;
//use crate::cellids::cellIDsError::NnuclError;
use std::io::BufWriter;
use std::fs::File;
use std::io::Write;

use std::path::PathBuf;

/// GeneCount is an entry in the CellData and therefore does not need to collect it's own id
pub struct GeneCount{
    pub umi: HashSet<u64>
}

impl GeneCount{
    pub fn new() ->Self{
        let umi = HashSet::new();
        Self{
            umi
        }
    }

    pub fn insert( &mut self, umi:u64 ){
        self.umi.insert( umi );
    }

}

/// CellData here is a storage for the total UMIs. UMIs will be checked per cell
/// But I do not correct the UMIs here - even with sequencing errors 
/// we should get a relatively correct picture as these folow normal distribution.
pub struct CellData{
    pub kmer_size: usize,
    pub name: std::string::String,
    pub gene: BTreeMap<u64, GeneCount>
}

impl CellData{
    pub fn new( kmer_size:usize, name: std::string::String ) -> Self{
        let gene =  BTreeMap::new(); // to collect the sample counts
        let loc_name = name.clone();
        Self{
            kmer_size,
            name: loc_name,
            gene
        }
    }


    // pub fn get(&mut self, seq: &[u8] ) -> Result< &GeneCount, &str>{
        
    //     //let min_value = 2;
    //     //let min_z = 1;
    //     //let mut max_value = 0;
    //     //let mut ret:u32 = 0;
    //     let kmers = needletail::kmer::Kmers::new(seq, self.kmer_size as u8 );
    //     let mut kmer_vec = Vec::<u64>::with_capacity(60);

    //     self.fill_kmer_vec(kmers, &mut kmer_vec);

    //     let id = 0;

    //     let km = kmer_vec[id];
    //     //println!("GeneIds::get - checking this sequence: {} and the at pos {}", km, id );
    //     let ret = match self.gene.get_mut(&km){
    //         Some(c1) => c1, 
    //         None => return Err::<&GeneCount, &str>( "Genes NoMatch"), 
    //     };
    //     Ok( ret )
    // }


    pub fn add(&mut self, geneid: u64, umi:u64 ){
        match self.gene.get_mut( &geneid) {
            Some( gene) => gene.insert( umi ),
            None => {
                let mut gc = GeneCount::new();
                gc.insert( umi );
                self.gene.insert( geneid, gc );
            }
        }
    }

    // fn fill_kmer_vec<'a>(&self, seq: needletail::kmer::Kmers<'a>, kmer_vec: &mut Vec<u64>) {
    //    kmer_vec.clear();
    //    let mut bad = 0;
    //    for km in seq {
    //         // I would like to add a try catch here - possibly the '?' works?
    //         // if this can not be converted it does not even make sense to keep this info
    //         for nuc in km{
    //             if *nuc ==b'N'{
    //                 bad = 1;
    //             }
    //         }
    //         if bad == 0{
    //             // let s = str::from_utf8(km);
    //             // println!( "this is the lib: {:?}",  s );
    //             kmer_vec.push(Kmer::from(km).into_u64());
    //         }
    //    }
    // }
    
    pub fn to_str<'live>(&mut self, genes:&GeneIds ) -> std::string::String {

        let mut data = Vec::<std::string::String>::with_capacity( genes.names.len()+3 );
        data.push(self.name.clone());

        let mut names = Vec::<std::string::String>::with_capacity( genes.names.len() );
        for (name, _id) in &genes.names{
            names.push( name.clone() )
        }
        let mut nums:Vec<u32>;
        let mut real_id:usize;

        nums = vec!(0 as u32; genes.names.len());

        for (id , info ) in &genes.kmers {
            real_id = match &genes.names.get( &info.name ){
                Some(num) => **num,
                None => panic!("I can not find the name {} in the genes object", info.name )
            };
            nums[ real_id-1 ] += match self.gene.get( id ){
                Some(gene_count) => gene_count.umi.len() as u32,
                None => 0 as u32
            }
        }

        //println!("I have a data vector wiith capacity {} and a nums vector with {}", genes.names.len()+1, nums.len());
        
        // get some statis and a little bit more
        let mut total:u32 = 0;
        let mut max:u32 = 0;

        for i in 0..nums.len() {
            total += nums[i];
            if  nums[i] > max{
                max = nums[i]
            }
            data.push( nums[i].to_string())
        }
        for i in 0..nums.len() {
            if nums[i] == max {
                data.push( names[i].clone());
                break; 
            }
        }
        data.push( (max as f32 / total as f32 ).to_string());

        let ret = data.join( "\t" );
        format!( "{}",ret)
    }
}



// This CellIds10x needs to copy some of the logics from split2samples - no it actually is totally dufferent
// Here we look for new sample ids and each sample id needs to be a total match to the previousely identified sample id
// Of cause I could also implement something with a whitelist. But that is for the future.
pub struct CellIds10x{    
    kmer_size: usize,
    //kmers: BTreeMap<u64, u32>,
    cells: BTreeMap<u64, CellData>
}


// here the functions
impl <'a> CellIds10x{

    pub fn new(kmer_size:usize )-> Self {

        let cells = BTreeMap::new();
        // let kmer_len = 5;
        // TACAGAACA

        Self {
            kmer_size,
            cells
        }
    }

    /// here the get checks for a complete match of the cell ID
    /// and if taht fails we need to add
    pub fn get(&mut self, cell_id: u64, name: std::string::String ) -> Result< &mut CellData, &str>{
        
        //println!("CellIDs::get cell_id: {}", cell_id );
        if ! self.cells.contains_key( &cell_id ){
            let data = CellData::new(self.kmer_size, name );
            self.cells.insert( cell_id, data );
        }

        let ret = match self.cells.get_mut(&cell_id){
            Some(c1) => c1, 
            None => return Err::< &mut CellData, &str>("BTreeMap Upstream error")
        };
        Ok( ret )
    }

    // pub fn add (&mut self, cell_id: u64, gene_id: u64, umi: u64, name: std::string::String  )-> Result< (), &str>{
    //     //let mut cell_id:CellData;

    //     match self.cells.get_mut( &cell_id ){
    //         Some(c1) => {
    //             c1.add( gene_id, umi );
    //         },
    //         None => {
    //             let mut data = CellData::new(self.kmer_size, name );
    //             data.add( gene_id, umi );
    //             self.cells.insert( cell_id, data );
    //         }
    //     };
    //     Ok( () )
    // }

    pub fn write (&mut self, file_path: PathBuf, genes: &GeneIds) -> Result< (), &str>{
        
        let file = match File::create( file_path ){
            Ok(file) => file,
            Err(err) => {
                panic!("Error: {:#?}", err);
            }
        };
        let mut writer = BufWriter::new(&file);

        match writeln!( writer, "{}", genes.to_header() ){
            Ok(_) => (),
            Err(err) => {
                eprintln!("write error: {}", err);
                return Err::<(), &str>("Header could not be written")
            }
        };

        for ( _id,  cell_obj ) in &mut self.cells {

            //println!( "get something here?: {}", cell_obj.to_str( &gene_ids ) );

            match writeln!( writer, "{}", cell_obj.to_str( genes )){
            // the compiler thought this might be more correct...
            //match writeln!( writer, "{}", cell_obj.to_str( <Vec<u64> as Borrow<Borrowed>>::borrow(&gene_ids).clone() ) ){
             Ok(_) => (),
             Err(err) => {
                 eprintln!("write error: {}", err);
                 return Err::<(), &str>("cell data could not be written")   
             }
            }
        }
        Ok( () )
    }
}


// macro_rules! to_cellid {
//    ($r1: expr, $r2: expr) => {
//       //1-9,22-30,44-52
//       $r1.to_cellid($r1, $r2, [0,8], [21,29], [43,51]  )
//    };
// }

// #[cfg(test)]
// mod tests {
//     #[test]
//     fn getsamples() {
//         let cells = crate::cell_ids::new();

//         let mut primer = b"GTCGCTATANNNNNNNNNNNNTACAGGATANNNNNNNNNNNNNAAGCCTTCT";
//         let mut id:u32 = 1;
//         let mut exp= ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
//         match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
//             Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
//             Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
//         };
        
//         let mut exp2 = vec![b"GTCGCTATA", b"TACAGGATA", b"AAGCCTTCT"];
//         assert_eq!( cells.to_sequence( exp ), exp2 );
//         // 3, 3, 3
//         primer = b"CTTCACATANNNNNNNNNNNNTGTGAAGAANNNNNNNNNNNNNCACAAGTAT";
//         id = 3;
//         exp = ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
//         exp2 = vec![b"CTTCACATA", b"TGTGAAGAA", b"CACAAGTAT"];
//         match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
//             Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
//             Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
//         };
//         //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
//         assert_eq!( cells.to_sequence( exp ), exp2 );

//         // and the last one
//         primer = b"TGCGATCTANNNNNNNNNNNNCAACAACGGNNNNNNNNNNNNNCATAGGTCA";
//         id = 96;
//         exp = ((id-1)* 96 * 96 + (id-1) * 96 + (id-1) +1) as u32;
//         exp2 = vec![b"TGCGATCTA", b"CAACAACGG", b"CATAGGTCA"];
//         assert_eq!( 884735+1 , exp);
//         match cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52]){
//             Ok(val) => assert_eq!( val , exp ), // will never insert one element twice. Great!
//             Err(err) => (), //we mainly need to collect cellids here and it does not make sense to think about anything else right now.
//         };
//         //assert_eq!( cells.to_cellid( primer, vec![0,9], vec![21,30], vec![43,52])? , exp );
//         assert_eq!( cells.to_sequence( exp ), exp2 );        
//     }
// }


