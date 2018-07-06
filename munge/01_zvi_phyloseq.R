################################################################################
##### Create phyloseq object with qiime output

library(ape)
library(Biostrings)
library(lubridate)

### Source data file names
biom_file <- "data/raw_data/otu_table_no_pynast_failures.biom"
seq_file <- "data/raw_data/rep_set.fna"
tree_file <- "data/raw_data/rep_set_aligned_pfiltered.tre"
metadata_file <- "data/raw_data/zvi_meta.csv"

################ Adding count and taxa data to phyloseq object #################
ps <- import_biom(biom_file)

################ Adding metadata to phyloseq object ############################
meta_df <- read.csv(metadata_file, stringsAsFactors = FALSE) %>%
    mutate(Collection_Date = mdy(Collection_Date)) %>% 
    column_to_rownames(var = "Sample_ID") 

## Reformat biom sample names to match metadata sample names
sample_names(ps) <- sample_names(ps) %>%
    str_replace_all("(?<=[:alpha:])\\.|\\.(?=[:alpha:])", "_") %>%
    str_replace_all("(?<=[:digit:])\\.(?=[:digit:])", "-")

## defining slot
sample_data(ps) <- meta_df 

################ Renaming OTUs #################################################
##
## Addresses issue with OTU names interpreted as numeric and character strings 
##

format_otus_names <- function(otus) {
    ## All numbers the same length
    otus <- str_pad(
        string = as.numeric(otus),
        width = max(nchar(otus)),
        side = "left",
        pad = "0"
    )
    
    ## Add characters to otu ids so they are interpreted as characters
    paste0("OTU_", otus)
}

biom_otus <- taxa_names(ps) %>% format_otus_names()
taxa_names(ps) <- biom_otus


################ Adding seq data to phyloseq object ############################

seq_dat <- readDNAStringSet(seq_file)

###### Reformating OTU names
names(seq_dat) <- names(seq_dat) %>% 
    str_replace(" .*","") %>% 
    format_otus_names()

###### Filtering rep seq set
filtered_seq_dat <- seq_dat[names(seq_dat) %in% biom_otus]


## Defining seq slot
ps@refseq <- filtered_seq_dat

################ Adding tree data to phyloseq object ###########################

tree_dat <- read.tree(tree_file)

###### Reformating OTU names
tree_dat$tip.label <- format_otus_names(tree_dat$tip.label)

## Defining tree slot
phy_tree(ps) <- tree_dat

################### Saving phyloseq object ################
ProjectTemplate::cache("ps")
