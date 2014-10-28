#!/usr/bin/env R
library(DESeq)
library(RColorBrewer)
library(gplots)

# Once we have all of the cell lines done we will make 
# comparisons for each drug of the most and least sensitive 
# to identify potential biomarkers. We have made similar 
# comparisons using microarray based analysis of the miRNAs
# in the NCI-60 human cancer cell line panel to generate a 
# list of potential miRNAs.
# I attached a copy of a grant that describes this process 
# a little. At this point, we just need to get relative 
# expression levels of the miRNA species present in these
# samples. Originally, we had planned to use a qRT-PCR 
# based process to get expression levels for 384 targets, 
# but decided to just go for broke and have sequencing done.
# I think the most difficult aspect at this point will be 
# identification of the miRNAs since far fewer canine miRNAs 
# have been identified and the canine genome is not as 
# complete as some.

ctable <- read.table("~/projects/canine/data/20130605/counts_matrix.txt", header=TRUE, row.names=1)

edesign <- data.frame(
    row.names = colnames(ctable),
    condition = colnames(ctable),
    libType = c(rep("1", length(colnames(ctable))))
)
cds = newCountDataSet(ctable, edesign)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds, method="blind")
vds = varianceStabilizingTransformation(cds)

select = order(rowMeans(counts(cds, normalized=TRUE)), decreasing=T)[1:1000]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vds)[select,], col=hmcol, trace="none", margin=c(10, 6))
# sample to sample distances
dists = dist(t(exprs(vds)))
mat = as.matrix(dists)
rownames(mat) = colnames(mat) = with(pData(cds), paste(edesign$condition, edesign$libType, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

# case view
plotPCA(vds, intgroup=c("condition"))
# replicate view
plotPCA(vds, intgroup=c("condition", "libType"))
