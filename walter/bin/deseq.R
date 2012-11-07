walterCountTable = read.table("/Users/brownj/projects/walter/data/countsmatrix.txt", header=T, row.names=1)
walterDesign = data.frame(  
    row.names = colnames( walterCountTable ),
    condition = c("inf","uninf","inf","uninf","inf",
                  "uninf","inf","uninf","inf","uninf",
                  "inf","uninf","inf","uninf","inf",
                  "uninf","inf","uninf","inf","uninf",
                  "inf","uninf","inf","uninf"),
    libType = c("1","1","24","24","2",
                "2","8","8","1","1",
                "24","24","2","2","8",
                "8","1","1","24","24",
                "2","2","8","8"))
library(DESeq)
cds = newCountDataSet(walterCountTable, walterDesign$condition)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
plotDispEsts(cds)
res = nbinomTest(cds, "uninf", "inf")
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="p-val dist")
write.csv(res, file="/Users/brownj/projects/walter/data/result.txt")

# QC
cdsFull = newCountDataSet(walterCountTable, walterDesign)
cdsFull = estimateSizeFactors(cdsFull)
cdsFull = estimateDispersions(cdsFull)
cdsFullBlind = estimateDispersions(cdsFull, method="blind")
vsdFull = varianceStabilizingTransformation(cdsFullBlind)
library(RColorBrewer)
library(gplots)
select = order(rowMeans(counts(cdsFull)), decreasing=T)[1:150]
hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6))
# sample to sample distances
dists = dist(t(exprs(vdsFull)))
mat = as.matrix(dists)
rownames(mat) = colnames(mat) = with(pData(cdsFullBlind), paste(walterDesign$condition, walterDesign$libType, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))