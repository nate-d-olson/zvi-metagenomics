## Create phyloseq object with qiime output
biom_file <- "raw_data/otu_table_no_pynast_failures.biom"
seq_file <- "raw_data/rep_set.fna"
tree_file <- "raw_data/rep_set_aligned_pfiltered.tre"


qiime_ps <- phyloseq::import_biom(biom_file)

## Need meta data 
#meta_df <- read_tsv("sample_metadata.tsv")
#sample_data(qiime_ps) <- meta_df

## Filtering rep seq set
seq_dat <- Biostrings::readDNAStringSet(seq_file)
names(seq_dat) <- gsub(pattern = " .*",replacement =  "",x = names(seq_dat))
biom_otus <- rownames(phyloseq::otu_table(qiime_ps))

filtered_seq_dat <- seq_dat[names(seq_dat) %in% biom_otus]

## saving filtered seq data  
filtered_seq_file <- "raw_data/rep_set_pfiltred.fna"
Biostrings::writeXStringSet(filtered_seq_dat, filtered_seq_file)

### Creating phyloseq object with sequences and phylogenetic tree
## Error with otu names - I think some might be interpreted as numeric and others as characters
ps <- phyloseq::import_biom(biom_file, 
                            # refseqfilename = filtered_seq_file,
                            treefilename = tree_file)
saveRDS(ps, "phyloseq.RDS")
