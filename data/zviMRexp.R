## Create MRexperiment from biom
## Use biomformat package and metagenomeSeq
ps <- readRDS("phyloseq.RDS")
mrexp <- phyloseq_to_metagenomeSeq(ps)
saveRDS(mrexp, "mrexp.RDS")
