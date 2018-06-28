library(tidyverse)
make_funct_anno_df <- function(){
    ## Get list of annotation files
    annotation_files <- list.files("data/raw_data/eggnog", 
                                   pattern = "annotations$",full.names = TRUE) %>% 
        set_names(str_extract(., "(?<=eggnog/).*(?=.pep.*)"))
    
    ## Define column names and types (for faster processing)
    anno_cols <- read_lines(annotation_files[[1]],n_max = 5)[[3]] %>% 
        str_remove_all("#") %>% 
        str_split(pattern = "\t") %>%  
        unlist() 
    anno_col_type <- "ccddcccccccc"
    
    ## Read files and concatenate into a single data frame
    annotation_files %>% 
        map_dfr(read_tsv, comment = "#", 
                col_names = anno_cols, 
                col_types = anno_col_type,
                .id = "Sample_ID")
}

funct_anno_df <- make_funct_anno_df()

## Generate and cache data frame
ProjectTemplate::cache("funct_anno_df", {make_funct_anno_df()})