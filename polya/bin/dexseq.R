# source("http://bioconductor.org/biocLite.R")
# biocLite("DEXSeq")
library(DEXSeq)
library(parallel)
# library(DESeq)
setwd("~/projects/polya/data/20130418/")

# annotation_file = "~/projects/polya/data/polya_extended_with_names.gff"

run_dexseq <- function(rownames, conditions, count_files, output){
    design = data.frame(row.names = rownames,
                        condition = conditions)
    cds = read.HTSeqCounts(countfiles = count_files,
                           design = design)
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds, minCount = 1, nCores=4, quiet=TRUE)
    cds = fitDispersionFunction(cds)
    cds = testForDEU(cds, nCores=4)
    cds = estimatelog2FoldChanges(cds)
    res = DEUresultTable(cds)
    write.table(res, file=output, sep="\t")
}

rownames = c("MP55", "MP56")
conditions = c("MP55", "MP55", "MP56", "MP56")
count_files = c("MP55.pos.counts", "MP56.pos.counts")
output = "TESTING.txt"

rownames=c("MP55", "MP55x", "MP56", "MP56x")
conditions=c("normal", "normal", "tumor", "tumor")
count_files=c("MP55.pos.counts", "MP55.pos.x.counts", "MP56.pos.counts", "MP56.pos.x.counts")


run_dexseq(c("MP55", "MP55x", "MP56", "MP56x"),
           c("normal", "normal", "tumor", "tumor"),
           c("MP55.pos.counts", "MP55.pos.x.counts", "MP56.pos.counts", "MP56.pos.x.counts"),
           "MP55_vs_MP56_pos.csv")

run_dexseq(c("nbt29","nbt29x","ts21","ts21x"),
           c("normal","normal","tumor","tumor"),
           c("nbt29.pos.counts","nbt29.posx.counts","ts21.pos.counts","ts21.posx.counts"),
           "nbt29_vs_ts21_pos.csv")
run_dexseq(c("nbt29","nbt29x","ts21","ts21x"),
           c("normal","normal","tumor","tumor"),
           c("nbt29.neg.counts","nbt29.negx.counts","ts21.neg.counts","ts21.negx.counts"),
           "nbt29_vs_ts21_neg.csv")

run_dexseq(c("nbt39","nbt39x","ts28","ts28x"),
           c("normal","normal","tumor","tumor"),
           c("nbt39.pos.counts","nbt39.posx.counts","ts28.pos.counts","ts28.posx.counts"),
           "nbt39_vs_ts28_pos.csv")
run_dexseq(c("nbt39","nbt39x","ts28","ts28x"),
           c("normal","normal","tumor","tumor"),
           c("nbt39.neg.counts","nbt39.negx.counts","ts28.neg.counts","ts28.negx.counts"),
           "nbt39_vs_ts28_neg.csv")

run_dexseq(c("nbt89","nbt89x","ts57","ts57x"),
           c("normal","normal","tumor","tumor"),
           c("nbt89.pos.counts","nbt89.posx.counts","ts57.pos.counts","ts57.posx.counts"),
           "nbt89_vs_ts57_pos.csv")
run_dexseq(c("nbt89","nbt89x","ts57","ts57x"),
           c("normal","normal","tumor","tumor"),
           c("nbt89.neg.counts","nbt89.negx.counts","ts57.neg.counts","ts57.negx.counts"),
           "nbt89_vs_ts57_neg.csv")
