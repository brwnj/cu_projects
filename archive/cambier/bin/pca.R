#!/usr/bin/env R
library(DESeq)
library(RColorBrewer)
library(gplots)

ctable <- read.csv("~/projects/archive/cambier/data/fpkm.csv", header=TRUE, row.names=1)
# need to round to integer
ctable <- round(ctable, 0)
edesign <- data.frame(
    row.names = colnames(ctable),
    condition = c("MD4ML5","MD4","ARS_A1","BL6B220","BL6ALPHA","BL6GAMMA","BL6FO","BL6ACT","ARS_A1","BL6","MD4","MD4ML5"),
    libType = c("1","1","1","1","1","1","1","1","2","2","2","2")
)
cds = newCountDataSet(ctable, edesign)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds, method="blind")
vds = varianceStabilizingTransformation(cds)

select = order(rowMeans(counts(cds)), decreasing=T)[1:1000]
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
