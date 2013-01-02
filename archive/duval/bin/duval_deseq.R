library(DESeq)

duvalCountTable = read.csv("~/projects/duval/data/sample_counts.csv",
                           header=T,
                           row.names=1)

# A vs D
# subset(world,select=income:military)
use = duvalCountTable[,c(3,4,5,6,7,8)]
design = data.frame(row.names = colnames(use),
                    condition = c("T220A","T220A","T220A",
                                  "T220D","T220D","T220D"))
cds = newCountDataSet(use, design$condition)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
plotDispEsts(cds)
res = nbinomTest(cds, "T220A", "T220D")
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="p-val dist")
write.csv(res, file="~/projects/duval/data/T220A_to_T220D.csv")

# A vs E
use = duvalCountTable[,c(3,4,5,12,13,14)]
design = data.frame(row.names = colnames(use),
                    condition = c("T220A","T220A","T220A",
                                  "T220E","T220E","T220E"))
cds = newCountDataSet(use, design$condition)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds, fitType="local")
plotDispEsts(cds)
res = nbinomTest(cds, "T220A", "T220E")
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="p-val dist")
write.csv(res, file="~/projects/duval/data/T220A_to_T220E.csv")

# E vs WT
use = duvalCountTable[,c(12,13,14,9,10,11)]
design = data.frame(row.names = colnames(use),
                    condition = c("T220E","T220E","T220E",
                                  "WT","WT","WT"))
cds = newCountDataSet(use, design$condition)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
plotDispEsts(cds)
res = nbinomTest(cds, "T220E", "WT")
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="p-val dist")
write.csv(res, file="~/projects/duval/data/T220E_to_WT.csv")

# D vs WT
use = duvalCountTable[,c(6,7,8,9,10,11)]
design = data.frame(row.names = colnames(use),
                    condition = c("T220D","T220D","T220D",
                                  "WT","WT","WT"))
cds = newCountDataSet(use, design$condition)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds, fitType="local")
plotDispEsts(cds)
res = nbinomTest(cds, "T220D", "WT")
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="p-val dist")
write.csv(res, file="~/projects/duval/data/T220D_to_WT.csv")

# A vs WT
use = duvalCountTable[,c(3,4,5,9,10,11)]
design = data.frame(row.names = colnames(use),
                    condition = c("T220A","T220A","T220A",
                                  "WT","WT","WT"))
cds = newCountDataSet(use, design$condition)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds, fitType="local")
plotDispEsts(cds)
res = nbinomTest(cds, "T220A", "WT")
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="p-val dist")
write.csv(res, file="~/projects/duval/data/T220A_to_WT.csv")

# QC
exp_sub = duvalCountTable[,c(3,4,5,6,7,8,9,10,11,12,13,14)]
duvalDesign = data.frame(row.names = colnames(exp_sub),
                         condition = c("T220A","T220A","T220A",
                                       "T220D","T220D","T220D",
                                       "WT","WT","WT",
                                       "T220E","T220E","T220E"),
                         libType = c("0","1","2",
                                     "0","1","2",
                                     "0","1","2",
                                     "0","1","2"))
cdsFull = newCountDataSet(exp_sub, duvalDesign)
cdsFull = estimateSizeFactors(cdsFull)
# cdsFull = estimateDispersions(cdsFull)
cdsFullBlind = estimateDispersions(cdsFull, method="blind", fitType="local")
vdsFull = varianceStabilizingTransformation(cdsFullBlind)
library(RColorBrewer)
library(gplots)
select = order(rowMeans(counts(cdsFull)), decreasing=T)[1:150]
hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(exprs(vdsFull)[select,], col = hmcol, trace="none", margin=c(10, 6))
# sample to sample distances
dists = dist(t(exprs(vdsFull)))
mat = as.matrix(dists)
rownames(mat) = colnames(mat) = with(pData(cdsFullBlind), paste(duvalDesign$condition, duvalDesign$libType, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
# case view
plotPCA(vdsFull, intgroup=c("condition"))
# replicate view
plotPCA(vdsFull, intgroup=c("condition", "libType"))