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

/// CellData here is a storage for the total UMIs. UMIs will be checked per cell
/// But I do not correct the UMIs here - even with sequencing errors 
/// we should get a relatively correct picture as these folow normal distribution.
pub struct CellData{
    pub kmer_size: usize,
    pub name: std::string::String,
    pub genes: BTreeMap<usize, HashSet<u64>>
}

impl CellData{
    pub fn new( kmer_size:usize, name: std::string::String ) -> Self{
        let genes =  BTreeMap::new(); // to collect the sample counts
        let loc_name = name.clone();
        Self{
            kmer_size,
            name: loc_name,
            genes
        }
    }

    pub fn add(&mut self, geneid: usize, umi:u64 ){
        //println!("adding gene id {}", geneid );
        match self.genes.get_mut( &geneid ) {
            Some( gene ) => {
                gene.insert( umi ); // the gene has already been added - check if umi matters
                }, 
            None => {
                let mut gc:HashSet<u64> = HashSet::new(); //to store the umis
                gc.insert( umi );
                self.genes.insert( geneid, gc );
            }
        }
    }
    
    pub fn to_str<'live>(&mut self, gene_info:&GeneIds ) -> std::string::String {

        let mut data = Vec::<std::string::String>::with_capacity( gene_info.names.len()+3 );
        data.push(self.name.clone());

        // here our internal data already should be stored with the same ids as the gene names.
        let mut total = 0;
        let mut max = 0;
        let mut max_name:std::string::String = "na".to_string();

        for (name, id) in &gene_info.names {
            //println!("I collect expression for id {}", id);
            match self.genes.get( id  ){
                Some(hash) => {
                    let n = hash.len();
                    total += n;
                    if  n > max{
                        max = n;
                        max_name = name.clone();
                    }
                    data.push( n.to_string() )
                },
                None => {
                    data.push( 0.to_string() )
                }
            }
        }
        data.push( max_name.clone() ); // max expressing gene (or sample id in an HTO analysis)
        data.push( (max as f32 / total as f32 ).to_string()); // fraction of reads for the max gene

        let ret = data.join( "\t" );
        format!( "{}",ret)
    }
}