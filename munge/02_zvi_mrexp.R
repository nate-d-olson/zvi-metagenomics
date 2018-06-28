## Create MRexperiment from biom
## Use biomformat package and metagenomeSeq
mrexp <- phyloseq_to_metagenomeSeq(ps)
ProjectTemplate::cache("mrexp", depends = c("ps"))
