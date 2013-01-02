library(DESeq)

# walterCountTable = read.table("~/projects/walter/data/20121120/countsmatrix.txt", header=T, row.names=1)

# tb
# walterCountTable = read.csv("~/projects/walter/data/20121120/sample_counts.csv", header=T, row.names=1)
walterCountTable = read.csv("~/projects/walter/data/20121124/sample_counts.csv", header=T, row.names=1)
# full experiment
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
cds = newCountDataSet(walterCountTable, walterDesign$condition)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
plotDispEsts(cds)
res = nbinomTest(cds, "uninf", "inf")
png("~/projects/walter/data/ma_tb_full.png", width=2000, height=2000)
plotMA(res)
dev.off()
png("~/projects/walter/data/pvals_tb_full.png", width=2000, height=2000)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="p-val dist")
dev.off()
write.table(res, file="/Users/brownj/projects/walter/data/20121124/result_tb_full.txt", row.names=F)

# normalize, estimate, and test only per time series
# 1 hour
t1c = subset(walterCountTable, select = c("E1T1_Inf","E1T1_Uninf",
                                          "E2T1_Inf","E2T1_Uninf",
                                          "E3T1_Inf","E3T1_Uninf"))
t1d = data.frame(row.names = colnames(t1c),
                 condition = c("inf","uninf",
                               "inf","uninf",
                               "inf","uninf"))
cds = newCountDataSet(t1c, t1d$condition)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
cds = estimateDispersions(cds)
png("~/projects/walter/data/dispersions_t1.png", width=2000, height=2000)
plotDispEsts(cds)
dev.off()
res = nbinomTest(cds, "uninf", "inf")
png("~/projects/walter/data/ma_t1.png", width=2000, height=2000)
plotMA(res)
dev.off()
png("~/projects/walter/data/pvals_t1.png", width=2000, height=2000)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="p-val dist")
dev.off()
write.csv(res, file="~/projects/walter/data/result_tb_t1.csv")

# 2 hours
t2c = subset(walterCountTable, select = c("E1T2_Inf","E1T2_Uninf",
                                          "E2T2_Inf","E2T2_Uninf",
                                          "E3T2_Inf","E3T2_Uninf"))
t2d = data.frame(row.names = colnames(t2c),
                 condition = c("inf","uninf",
                               "inf","uninf",
                               "inf","uninf"))
cds = newCountDataSet(t2c, t2d$condition)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
cds = estimateDispersions(cds)
png("~/projects/walter/data/dispersions_t2.png", width=2000, height=2000)
plotDispEsts(cds)
dev.off()
res = nbinomTest(cds, "uninf", "inf")
png("~/projects/walter/data/ma_t2.png", width=2000, height=2000)
plotMA(res)
dev.off()
png("~/projects/walter/data/pvals_t2.png", width=2000, height=2000)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="p-val dist")
dev.off()
write.csv(res, file="~/projects/walter/data/result_tb_t2.csv")

# 8 hours
t8c = subset(walterCountTable, select = c("E1T8_Inf","E1T8_Uninf",
                                          "E2T8_Inf","E2T8_Uninf",
                                          "E3T8_Inf","E3T8_Uninf"))
t8d = data.frame(row.names = colnames(t8c),
                 condition = c("inf","uninf",
                               "inf","uninf",
                               "inf","uninf"))
cds = newCountDataSet(t8c, t8d$condition)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
cds = estimateDispersions(cds)
png("~/projects/walter/data/dispersions_t8.png", width=2000, height=2000)
plotDispEsts(cds)
dev.off()
res = nbinomTest(cds, "uninf", "inf")
png("~/projects/walter/data/ma_t8.png", width=2000, height=2000)
plotMA(res)
dev.off()
png("~/projects/walter/data/pvals_t8.png", width=2000, height=2000)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="p-val dist")
dev.off()
write.csv(res, file="~/projects/walter/data/result_tb_t8.csv")

# 24 hours
t24c = subset(walterCountTable, select = c("E1T24_Inf","E1T24_Uninf",
                                          "E2T24_Inf","E2T24_Uninf",
                                          "E3T24_Inf","E3T24_Uninf"))
t24d = data.frame(row.names = colnames(t24c),
                  condition = c("inf","uninf",
                                "inf","uninf",
                                "inf","uninf"))
cds = newCountDataSet(t24c, t24d$condition)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
cds = estimateDispersions(cds)
png("~/projects/walter/data/dispersions_t24.png", width=2000, height=2000)
plotDispEsts(cds)
dev.off()
res = nbinomTest(cds, "uninf", "inf")
png("~/projects/walter/data/ma_t24.png", width=2000, height=2000)
plotMA(res)
dev.off()
png("~/projects/walter/data/pvals_t24.png", width=2000, height=2000)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="p-val dist")
dev.off()
write.csv(res, file="~/projects/walter/data/result_tb_t24.csv")

save.image("~/projects/walter/data/deseq_session.RData")

# QC
cdsFull = newCountDataSet(walterCountTable, walterDesign)
cdsFull = estimateSizeFactors(cdsFull)
cdsFull = estimateDispersions(cdsFull)
cdsFullBlind = estimateDispersions(cdsFull, method="blind")
vdsFull = varianceStabilizingTransformation(cdsFullBlind)
library(RColorBrewer)
library(gplots)
select = order(rowMeans(counts(cdsFull)), decreasing=T)[1:1000]
hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(exprs(vdsFull)[select,], col = hmcol, trace="none", margin=c(10, 6))
# sample to sample distances
dists = dist(t(exprs(vdsFull)))
mat = as.matrix(dists)
rownames(mat) = colnames(mat) = with(pData(cdsFullBlind), paste(walterDesign$condition, walterDesign$libType, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

# case view
plotPCA(vdsFull, intgroup=c("condition"))
# replicate view
plotPCA(vdsFull, intgroup=c("condition", "libType"))