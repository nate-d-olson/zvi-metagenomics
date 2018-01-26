setwd("C:/Users/prachi.kulkarni/My Documents/R/zvi-metagenomics")
library("phyloseq"); packageVersion("phyloseq")
library("data.table"); packageVersion("data.table")
library("ggplot2"); packageVersion("ggplot2")
mp = readRDS("phyloseq.RDS")
print(mp)
nsamples(mp)
ntaxa(mp)
sample_data(mp)
tax_table(mp)
sample_variables(mp)
rank_names(mp)
phy_tree(mp)
sample_sums(mp)
sdt = data.table(as(sample_data(mp), "data.frame"),
                 TotalReads = sample_sums(mp), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
pSeqDepth + facet_wrap(~Sample_Code)
pSeqDepth + 
    facet_grid(Sample_Type ~ Sample_Code)
tdt = data.table(tax_table(mp),
                 TotalCounts = taxa_sums(mp),
                 OTU = taxa_names(mp))
ggplot(tdt, aes(TotalCounts)) + 
    geom_histogram() + 
    ggtitle("Histogram of Total Counts")
tdt[(TotalCounts <= 0), .N]
tdt[(TotalCounts <= 1), .N]
tdt[(TotalCounts <= 2), .N]
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
pCumSum + xlim(0, 100)
source("taxa_summary.R", local = TRUE)
mdt = fast_melt(mp)
prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
source("taxa_summary.R", local = TRUE)
mdt = fast_melt(mp)
prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
ggplot(prevdt, aes(Prevalence)) + 
    geom_histogram() + 
    ggtitle("Histogram of Taxa Prevalence")
prevdt[(Prevalence <= 0), .N]
prevdt[(Prevalence <= 1), .N]
prevdt[(Prevalence <= 2), .N]
prevcumsum = prevdt[, .N, by = Prevalence]
setkey(prevcumsum, Prevalence)
prevcumsum[, CumSum := cumsum(N)]
pPrevCumSum = ggplot(prevcumsum, aes(Prevalence, CumSum)) + 
    geom_point() +
    xlab("Filtering Threshold, Prevalence") +
    ylab("OTUs Filtered") +
    ggtitle("OTUs that would be filtered vs. the minimum count threshold")
pPrevCumSum
ggplot(prevdt, aes(Prevalence, TotalCounts)) + 
    geom_point(size = 4, alpha = 0.75) + 
    scale_y_log10()
keepTaxa = prevdt[(Prevalence >= 10 & TotalCounts > 3), TaxaID]
mpf1 = prune_taxa(keepTaxa, mp)
tipg = tip_glom(mpf1, h = 0.05)
# Transform to relative abundance
tipg <- transform_sample_counts(tipg, function(x) x / sum(x))
ntaxa(tipg)
plot_tree(tipg, size = "Abundance",
          color = "Sample_Type",
          justify = "yes please", 
          ladderize = "left") +
    scale_size_continuous(range = c(1, 3))
plot_tree(tipg, size = "Abundance",
          color = "Sample_Code",
          justify = "yes please", 
          ladderize = "left") +
    scale_size_continuous(range = c(1, 3))