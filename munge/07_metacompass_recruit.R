library(tidyverse)
make_metacompass_df <- function(){
    ## Get list of annotation files
    id_files <- list.files("data/raw_data/metacompass",
                                   recursive = TRUE,
                                   pattern = "ids$",
                                   full.names = TRUE) %>% 
        set_names(str_extract(., "(?<=use1X/).*(?=_/asm)"))
    
    ## Define column names and types (for faster processing)
    col_names <- c("acc","len","reads","org_name")
    col_type <- "cddc"
    
    ## Read files and concatenate into a single data frame
    id_files %>% 
        map_dfr(read_tsv,  
                col_names = col_names, 
                col_types = col_type,
                .id = "Sample_ID")
}


## Generate and cache data frame
ProjectTemplate::cache("asm_id_df", {make_metacompass_df()})