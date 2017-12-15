Here you'll store your raw data files. If they are encoded in a supported file format, they'll automatically be loaded when you call `load.project()`.

Store larger datafiles in appropriate formats to facilitate use, e.g. RDS for MRexperiments. 

__Raw Data Files__
* Provide brief description of each data file in this document.

__Data Files__
_Added_

_To Add_
* E. coli MPN - csv
* MRexperiment - RDS: convert from biom
* phyloseq - RDS?- with tree and seq data
* Alpha diversity - data frame RDS
    * richness
    * evenness
* Beta diversity - matrix RDS
    * Unifrac unweighted
    * Unifrac weighted  
    * Bray-Curtis  
    * Jacard  
* Relative abundance - data frame RDS
* Pathogen genomes  
* Functional analysis

__Data Preprocessing scripts__ - script to generate each of the data files listed above
* biom to MRexperiment  
* biom and tree to phyloseq  
* Alpha diversity richness  
* Alpha diversity evenness  
* Beta diversity  
* Relative abundance  

Use the project template cache function to add data to project template