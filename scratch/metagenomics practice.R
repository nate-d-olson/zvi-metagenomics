# Set working directory to project directory
setwd("C:/Users/prachi.kulkarni/My Documents/R/metagenomics_practice")

# Load required packages
library(ape)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(phangorn)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)
library(phyloseq)
library(data.table)
library(readr)
library(tidyr)
library(stringr)
library(magrittr)
library(tibble)
library(knitr)
library(devtools)
library(metagenomeSeq)
library(broom)
library(iNEXT)
library(ggExtra)
library(DT)
library(purrr)
library(ggfortify)

# Load phyloseq object (ZVI phyloseq Object "zpo")
zpo = readRDS("phyloseq.RDS")

# Subset phyloseq object so that only water samples CE, RB, RW, ZW and TW are included
#zpow1 = subset_samples(zpo, Sample_Code != "NI" & Sample_Code != "BG" & Sample_Code != "NR"& Sample_Code != "NC"
                    # & Sample_Code != "RL"& Sample_Code != "ZL"& Sample_Code != "TL"
                    # & Sample_Code != "RS"& Sample_Code != "ZS"& Sample_Code != "TS")

# Subset phyloseq object so that only water samples RW, ZW and TW are included
#zpow = subset_samples(zpo, Sample_Code != "NI" & Sample_Code != "BG" & Sample_Code != "NR"& Sample_Code != "NC"
                    #  & Sample_Code != "RL"& Sample_Code != "ZL"& Sample_Code != "TL"
                     # & Sample_Code != "RS"& Sample_Code != "ZS"& Sample_Code != "TS"
                     # & Sample_Code !="CE" & Sample_Code !="RB")

# Subset phyloseq object so that only lettuce samples RL,ZL and TL are included
zpol = subset_samples(zpo, Sample_Code != "NI" & Sample_Code != "BG" & Sample_Code != "NR"& Sample_Code != "NC"
                     & Sample_Code != "RW"& Sample_Code != "ZW"& Sample_Code != "TW"& Sample_Code != "CE"& Sample_Code != "RB"
                      & Sample_Code != "RS"& Sample_Code != "ZS"& Sample_Code != "TS")

# Subset phyloseq object so that only soil samples RS,ZS and TS are included
#zpol = subset_samples(zpo, Sample_Code != "NI" & Sample_Code != "BG" & Sample_Code != "NR"& Sample_Code != "NC"
#  & Sample_Code != "RW"& Sample_Code != "ZW"& Sample_Code != "TW"& Sample_Code != "CE"& Sample_Code != "RB"
# & Sample_Code != "RL"& Sample_Code != "ZL"& Sample_Code != "TL")

# Subset phyloseq object so that only background samples BG,NI and NR are included
#zpol = subset_samples(zpo, Sample_Code != "NC"
#  & Sample_Code != "RW"& Sample_Code != "ZW"& Sample_Code != "TW"& Sample_Code != "CE"& Sample_Code != "RB"
# & Sample_Code != "RL"& Sample_Code != "ZL"& Sample_Code != "TL")

# Set up variable for summarizing sequencing depths
#zo = zpow
#zo = zpow1
zo = zpol
#zo = zpos
#zo = zpob

# Summarize data
zo
nsamples(zo)
ntaxa(zo)
sample_variables(zo)
rank_names(zo)

# Summarize sequencing depths 
zdt = data.table(as(sample_data(zo), "data.frame"),
                 TotalReads = sample_sums(zo), keep.rownames = TRUE)
pSeqDepth = ggplot(zdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
pSeqDepth + facet_wrap(~Sample_Code)


#Sequencing depth across time
ggplot(zdt, aes(Collection_Date, TotalReads)) + 
  geom_point(size = 5) + 
  geom_smooth(method = lm) +
  ggtitle("Sequencing Depth vs. Time") +
  scale_y_log10() +
  facet_grid(Sample_Code ~ .)

# Filtering taxa
# Taxa total counts histogram
tdt = data.table(tax_table(zo),
                 TotalCounts = taxa_sums(zo),
                 OTU = taxa_names(zo))
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + 
  ggtitle("Histogram of Total Counts")

# How many singleton taxa "singletons"?
tdt[(TotalCounts <= 0), .N]
tdt[(TotalCounts <= 1), .N]

# How many doubleton taxa "doubletons"?
tdt[(TotalCounts <= 2), .N]

# Taxa cumulative sum
taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]
# Define the plot
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Minimum Total Counts") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")
pCumSum

# Zoom-in
pCumSum + xlim(0, 500)

# Taxa prevalence histogram
#Prevalence here as the number of times an OTU is observed at least once. 
# That is, it is the number of samples in which each OTU was non-zero.
source("taxa_summary.R", local = TRUE)
zdt = fast_melt(zo)
prevdt = zdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]

ggplot(prevdt, aes(Prevalence)) + 
  geom_histogram() + 
  ggtitle("Histogram of Taxa Prevalence")

# How many singletons?
prevdt[(Prevalence <= 0), .N]

prevdt[(Prevalence <= 1), .N]

# How many doubletons?
prevdt[(Prevalence <= 2), .N]

# taxa cumulative sum
prevcumsum = prevdt[, .N, by = Prevalence]
setkey(prevcumsum, Prevalence)
prevcumsum[, CumSum := cumsum(N)]
pPrevCumSum = ggplot(prevcumsum, aes(Prevalence, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Prevalence") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")
pPrevCumSum

#Prevalence vs. Total Count Scatter plot
ggplot(prevdt, aes(Prevalence, TotalCounts)) + 
  geom_point(size = 4, alpha = 0.75) + 
  scale_y_log10()

# Add the top 9 phylum to the plot above
addPhylum = unique(copy(zdt[, list(TaxaID, Rank2)]))
# Join by TaxaID
setkey(prevdt, TaxaID)
setkey(addPhylum, TaxaID)
prevdt <- addPhylum[prevdt]
showPhyla = prevdt[, sum(TotalCounts), by = Rank2][order(-V1)][1:9]$Rank2
setkey(prevdt, Rank2)
ggplot(prevdt[showPhyla], 
       mapping = aes(Prevalence, TotalCounts, color = Rank2)) + 
  geom_point(size = 4, alpha = 0.75) + 
  scale_y_log10()

# Tree Plot
keepTaxa = prevdt[(Prevalence >= 10 & TotalCounts > 3), TaxaID]
mpf1 = prune_taxa(keepTaxa, zo)
tipg = tip_glom(mpf1, h = 0.05)
# Transform to relative abundance
tipg <- transform_sample_counts(tipg, function(x) x / sum(x))
ntaxa(tipg)

plot_tree(tipg, size = "Abundance",
          color = "Sample_Code",
          justify = "yes please", 
          ladderize = "left") +
  scale_size_continuous(range = c(1, 3))

#Alpha Diversity
pAlpha = plot_richness(zo,
                       color = "Sample_Code",
                       measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
                       title = "Alpha Diveristy")
pAlpha + geom_point(size = 5)

# as comparing boxplots
plot_richness(zo, x = "Sample_Code", color = "Sample_Code") + geom_boxplot()



# Store as a new data variable
alphadt = data.table(pAlpha$data)
# Subset to just the Shannon index
alphadt <- alphadt[(variable == "Shannon")]
# Order by Days
alphadt <- alphadt[order(Collection_Date)][(is.na(se))]
# Define the plot
ggplot(data = alphadt, 
       mapping = aes(Collection_Date, value,
                     color = Sample_Code)) +
  geom_point(size = 5) + 
  geom_path() +
  facet_wrap(~variable, ncol = 2, scales = "free_y") +
  ylab("Shannon Index") +
  ggtitle("Shannon Index")






# Beta Diversity
# Define taxa to keep.
keepTaxa = prevdt[(Prevalence >= 10 & TotalCounts > 3), TaxaID]
# Define new object with relative abundance
mpra = transform_sample_counts(zo, function(x) x / sum(x))
# Filter this new object
mpraf = prune_taxa(keepTaxa, mpra)
# Calculate distances
DistBC = distance(mpraf, method = "bray")
DistUF = distance(mpraf, method = "wUniFrac")
ordBC = ordinate(mpraf, method = "PCoA", distance = DistBC)
ordUF = ordinate(mpraf, method = "PCoA", distance = DistUF)
plot_scree(ordBC, "Scree Plot: Bray-Curtis MDS")
plot_scree(ordUF, "Scree Plot: Weighted UniFrac MDS")
plot_ordination(mpraf, ordBC, color = "Sample_Code") +
  ggtitle("PCoA: Bray-Curtis")
plot_ordination(mpraf, ordUF, color = "Sample_Code") + 
  ggtitle("PCoA: Weigthed Unifrac")






### Beta Diversity Analysis - Nate's Code for WWTP Community Analysis
dist_methods <- unlist(distanceMethodList)[c(8,10)]

zo_pcoa_list <- list()
zo_nmds_list <- list()
for (i in dist_methods) {
  # Calculate distance matrix
  iDist <- distance(zo, method = i)
  # Calculate ordination
  zo_pcoa_list[[i]] <- ordinate(zo, method = "PCoA", distance = iDist)
  zo_nmds_list[[i]] <- ordinate(zo, method = "NMDS", distance = iDist)
}



### Figure 3 - Ordination Plots


plot_ordination(zo, zo_pcoa_list[[1]], color= "Sample_Code", title = "Bray Curtis PCoA") + theme_bw()+stat_ellipse()



plot_ordination(zo, zo_nmds_list[[1]], color= "Sample_Code", title = "Bray Curtis NMDS") + theme_bw()+stat_ellipse()
stressplot(zo_nmds_list[[1]])



#### Jaccard

plot_ordination(zo,zo_pcoa_list[[2]], color= "Sample_Code", title = "Jaccard PCoA") + theme_bw()+stat_ellipse()
plot_ordination(zo, zo_nmds_list[[2]], color= "Sample_Code", title = "Jaccard NMDS") + theme_bw()+stat_ellipse()
stressplot(zo_nmds_list[[2]])



### Testing for Differences


iDist <- distance(zo, method = "bray")
anosim(iDist,grouping = factor(zo@sam_data$Sample_Code))






## betadisper has a post-hoc test
iDist <- distance(zo, method = "bray")
beta_fit <- betadisper(iDist,group = factor(zo@sam_data$Sample_Code))




anova(beta_fit) %>% tidy() %>% kable()




TukeyHSD(beta_fit)$group %>% as.data.frame() %>% 
  rownames_to_column(var = "comparison") %>% 
  filter(`p adj` < 0.05) %>% kable(digits = 3)



# Load mrexperiment object (ZVI mr experiment Object "zmo")
  zmo = readRDS("mrexp.RDS")
  
  # Subset mrexperiment object so that only water samples RW, ZW and TW are included
  zmow = subset_samples(zpo, Sample_Code != "NI" & Sample_Code != "BG" & Sample_Code != "NR"& Sample_Code != "NC"
   & Sample_Code != "RL"& Sample_Code != "ZL"& Sample_Code != "TL"
   & Sample_Code != "RS"& Sample_Code != "ZS"& Sample_Code != "TS"
   & Sample_Code !="CE" & Sample_Code !="RB") 

  aggregated_species <-  cumNorm(aggregateByTaxonomy(zmo, lvl="Rank7"), p = 0.75)
  
  aggregation_level <- "Rank3"
  aggregated_zmo <- aggregateByTaxonomy(zmo, lvl=aggregation_level)
  
  normed_zmo <-  cumNorm(aggregated_zmo, p = 0.75)
  zmo_sample_data <-  pData(normed_zmo)
  mod <-  model.matrix(~Sample_Code, data = zmo_sample_data)
  results_zmo <-  fitFeatureModel(normed_zmo, mod)
  
  
  
  
  count_dat <- metagenomeSeq::expSummary(zmo)
  count_dat$Seq_ID <- rownames(count_dat)
  
  count_tbl <- zmo@assayData$counts %>% as_tibble() 
  count_tbl$OTU <- zmo@featureData@data$OTU
  count_tbl <- count_tbl %>% gather("sam","count",-OTU) %>% 
    mutate(count = if_else(count != 0, 1,0)) %>% 
    group_by(OTU) %>% summarise(n_sam = sum(count)) %>% 
    filter(n_sam == 1)
  assigned_un_df <- count_tbl %>% left_join(zmo@featureData@data) %>% 
    mutate(assinged_un = if_else(grepl("OTU_",Taxonomy), "unassigned","assigned")) %>% 
    group_by(assinged_un) %>% summarise(count = n())
  
  # total seq count
  count_total <- sum(count_dat$libSize)
  # total number of samples
  num_samples <- nrow(count_total)
  
  # average seq count per sample
  avg_seq_sample <- count_total/num_samples
  
  # number of unique assigned species level OTUs
  n_unique_assigned_otus <- assigned_un_df$count[assigned_un_df$assinged_un == "assigned"]
  # number of unique unassigned species level OTUs
  n_unique_unassigned_otus <- assigned_un_df$count[assigned_un_df$assinged_un == "unassigned"]

  sample_count_summary <- count_total %>% select(Sample_Type, Sample_Code) %>% 
    group_by(Sample_Code) %>% 
    summarize(count = n())
  sample_count_summary$count %>% sum()
  sample_count_summary %>% kable()
  
  
  div_est <- zmo  %>% MRcounts()  %>% DataInfo()
  div_est$Stage <- factor(count_total$Sample_Code)
  p <- ggplot(div_est, aes(x = n, y = SC)) +
    geom_point(aes(color = Sample_Code)) + scale_x_log10() + 
    geom_vline(aes(xintercept = 100), linetype = 2, color = "grey 60") + 
    theme_bw() +
    labs(x = "Number of Sequences",y = "Coverage") +
    theme(legend.position = "bottom") 
  ggMarginal(p, type = "histogram", margins = "x")