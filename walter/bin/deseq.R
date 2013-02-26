library(DESeq)

walterCountTable <- read.csv("~/projects/walter/data/20121120/sample_counts.csv", header=TRUE, row.names=1)
# full experiment
walterDesign <- data.frame(
    row.names = colnames(walterCountTable),
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
cds = newCountDataSet(walterCountTable, walterDesign$condition)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
res = nbinomTest(cds, "uninf", "inf")
plotDispEsts(cds)
# png("~/projects/walter/data/ma_tb_full.png", width=2000, height=2000)
plotMA(res)
# dev.off()
# png("~/projects/walter/data/pvals_tb_full.png", width=2000, height=2000)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="p-val dist")
# dev.off()
# write.table(res, file="/Users/brownj/projects/walter/data/20121124/result_tb_full.txt", row.names=F)

use <- res[res$padj < 0.001 & !is.na(res$padj),]

t1unf <- c("E1T1_Uninf","E2T1_Uninf","E3T1_Uninf")
t1inf <- c("E1T1_Inf","E2T1_Inf","E3T1_Inf")
t2unf <- c("E1T2_Uninf","E2T2_Uninf","E3T2_Uninf")
t2inf <- c("E1T2_Inf","E2T2_Inf","E3T2_Inf")
t8unf <- c("E1T8_Uninf","E2T8_Uninf","E3T8_Uninf")
t8inf <- c("E1T8_Inf","E2T8_Inf","E3T8_Inf")
t24unf <- c("E1T24_Uninf","E2T24_Uninf","E3T24_Uninf")
t24inf <- c("E1T24_Inf","E2T24_Inf","E3T24_Inf")

colSums(walterCountTable[use$id,t1unf])
colSums(walterCountTable[use$id,t1inf])





# QC

# cdsFull = newCountDataSet(walterCountTable, walterDesign)
# cdsFull = estimateSizeFactors(cdsFull)
# cdsFull = estimateDispersions(cdsFull)
# cdsFullBlind = estimateDispersions(cdsFull, method="blind")
# vdsFull = varianceStabilizingTransformation(cdsFullBlind)
# library(RColorBrewer)
# library(gplots)
# select = order(rowMeans(counts(cdsFull)), decreasing=T)[1:1000]
# hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
# heatmap.2(exprs(vdsFull)[select,], col = hmcol, trace="none", margin=c(10, 6))
# # sample to sample distances
# dists = dist(t(exprs(vdsFull)))
# mat = as.matrix(dists)
# rownames(mat) = colnames(mat) = with(pData(cdsFullBlind), paste(walterDesign$condition, walterDesign$libType, sep=" : "))
# heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
# 
# # case view
# plotPCA(vdsFull, intgroup=c("condition"))
# # replicate view
# plotPCA(vdsFull, intgroup=c("condition", "libType"))