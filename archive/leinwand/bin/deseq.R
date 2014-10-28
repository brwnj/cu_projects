library(DESeq)

countTable = read.csv(
    "~/remote/projects/leinwand/results/common/counts_table.txt",
    header=TRUE, sep="\t",
    row.names=1)
# replace NA with zeroes
countTable[is.na(countTable)] = 0
design = data.frame(
    row.names = colnames(countTable),
    condition = c("mdx","mdx","mdx",
                  "wt","wt","wt")
    )

rs = rowSums(countTable)
rm = rowMeans(countTable)
use = countTable[rm>5,]

cds = newCountDataSet(countTable, design$condition)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
plotDispEsts(cds)
res = nbinomTest(cds, "mdx", "wt")
png("~/projects/leinwand/data/MA_mdx_vs_wt.png",
    width=2000,
    height=2000)
plotMA(res)
dev.off()
png("~/projects/leinwand/data/PVALS_mdx_vs_wt.png",
    width=2000,
    height=2000)
hist(res$pval,
     breaks=100,
     col="skyblue",
     border="slateblue",
     main="p-val dist")
dev.off()
write.table(res, file="~/projects/leinwand/data/RESULT_mdx_vs_wt.txt", sep="\t", row.names=FALSE)
save.image("~/projects/leinwand/data/deseq_session.RData")

# QC
cdsBlind = estimateDispersions(cds, method="blind")
vdsFull = varianceStabilizingTransformation(cdsBlind)
library(RColorBrewer)
library(gplots)
select = order(rowMeans(counts(cds)), decreasing=T)[1:100]
hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(exprs(vdsFull)[select,], col = hmcol, trace="none", margin=c(10, 6))

# sample to sample distances
dists = dist(t(exprs(vdsFull)))
mat = as.matrix(dists)
rownames(mat) = colnames(mat) = with(pData(cdsBlind), paste(design$condition, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

# case view
plotPCA(vdsFull, intgroup=c("condition"))

# corrplot
library(corrplot)
cdsNorm = newCountDataSet(countTable, conditions=c(rep("t",3),c(rep("c",3))))
cdsNorm = estimateSizeFactors(cdsNorm)
normcounts = counts(cdsNorm, normalized=TRUE)
write.table(normcounts, file="normalized_counts.txt", sep="\t", row.names=FALSE)

cor(normcounts)
mat = cor(normcounts)
corrplot(mat, method="circle", order="alphabet", cl.lim=c(0,1))