# Usage:
# Rscript norm_count_table.R file_in.tab file_out.tab
 
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq")
library(DESeq)
 
args <- commandArgs(TRUE)
filein <- args[1]
fileout <- args[2]
 
ct <- read.table(filein, header=TRUE, row.names=1, sep="\t")
cds <- newCountDataSet(ct, condition=c(rep("t", ncol(ct))))
cds <- estimateSizeFactors(cds)
write.table(counts(cds, normalized=TRUE), file=fileout, sep="\t", col.names=NA)
