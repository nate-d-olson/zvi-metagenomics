# Alpha diversity richness estimates - all samples
require(breakaway)

## Calculate richness for all samples using breakaway
ps_to_breakaway_df <- function(ps) {
    freq_df <- otu_table(ps) %>%
        as.data.frame() %>%
        rownames_to_column(var = "otu_id")  %>%
        gather("sample_id", "abundance", -otu_id) %>%
        group_by(sample_id, abundance) %>%
        summarise(freq = n()) %>%
        filter(abundance > 1) %>%
        group_by(sample_id) %>%
        nest()
    
    break_fit <- freq_df$data %>%
        map(as.data.frame) %>%
        map(
            breakaway_nof1,
            print = FALSE,
            answer = TRUE,
            plot = FALSE,
            force = TRUE
        )
    
    ## Generate a data frame with breakaway richness estimates (est), standard error
    ## (seest), and confidence interval (lower, upper).
    break_est <- break_fit %>%
        set_names(freq_df$sample_id)  %>%
        map( ~ .$est) %>%
        keep(is.numeric) %>% as_tibble() %>%
        gather("sample_id", "est")
    
    break_seest <- break_fit %>%
        set_names(freq_df$sample_id)  %>%
        map( ~ .$seest) %>%
        keep(is.numeric) %>% as_tibble() %>%
        gather("sample_id", "seest")
    
    break_ci <- break_fit %>%
        set_names(freq_df$sample_id)  %>%
        map( ~ .$ci) %>%
        keep(is.numeric) %>% as_tibble() %>%
        mutate(ci = c("lower", "upper")) %>%
        gather("sample_id", "value", -ci) %>%
        spread(ci, value)
    
    break_model <- break_fit %>%
        set_names(freq_df$sample_id) %>%
        map( ~ .$name) %>%
        keep(is.character) %>% as_tibble() %>%
        gather("sample_id", "model")
    
    ## Combine into single data frame
    list(break_est, break_seest, break_ci, break_model) %>%
        Reduce(function(dtf1, dtf2)
            left_join(dtf1, dtf2), .)
    
}

ProjectTemplate::cache("richness_df", {ps_to_breakaway_df(ps)}, depends = "ps")