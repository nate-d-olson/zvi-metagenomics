make_funct_mrexp <- function(funct_anno_df){
    ## Functional Counts
    wide_fun_counts <- funct_anno_df %>% 
        ## Sample/ gene counts
        group_by(Sample_ID, predicted_gene_name) %>% summarize(count = n()) %>% 
        ## Converting to matrix
        spread(Sample_ID, count, fill = 0) %>% 
        mutate(predicted_gene_name = if_else(is.na(predicted_gene_name), 
                                             "uncategorized", 
                                             predicted_gene_name))
    fun_mat <- ungroup(wide_fun_counts) %>% 
        as.data.frame() %>% 
        column_to_rownames(var = "predicted_gene_name") %>% 
        as.matrix()

    ## Sample metadata
    pheno_df <- data_frame(Sample_ID = colnames(fun_mat)) %>% 
        separate(Sample_ID, 
                 c("Sample_Code", "Collection_Date", "Sample_Type"), 
                 sep = "_", remove = FALSE) %>%   
        #Fix date format
        mutate(Collection_Date = mdy(paste0(Collection_Date,"-2017")))
    pd <- as.data.frame(pheno_df)
    rownames(pd) <- pd$Sample_ID
    
    ## Generate MRexp
    newMRexperiment(fun_mat, phenoData = AnnotatedDataFrame(pd))
}

## Generate and cache mrexp with functional data
ProjectTemplate::cache("funct_mrexp", {make_funct_mrexp(funct_anno_df)})