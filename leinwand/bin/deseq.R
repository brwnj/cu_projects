library(DESeq)

countTable = read.table(
    "~/projects/leinwand/data/sample_counts.txt",
    header=T,
    row.names=1)

design = data.frame(  
    row.names = colnames(countTable),
    condition = c("mdx","mdx","mdx",
                  "wt","wt","wt")
    )

rs = rowSums(countTable)
rm = rowMeans(countTable)
use = countTable[rm>2,]

cds = newCountDataSet(use, design$condition)
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
write.table(res, file="~/projects/leinwand/data/RESULT_mdx_vs_wt.txt", row.names=F)
save.image("~/projects/leinwand/data/deseq_session.RData")

# QC
cdsBlind = estimateDispersions(cds, method="blind")
vdsFull = varianceStabilizingTransformation(cdsBlind)
library(RColorBrewer)
library(gplots)
select = order(rowMeans(counts(cds)), decreasing=T)[1:1000]
hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(exprs(vdsFull)[select,], col = hmcol, trace="none", margin=c(10, 6))
# sample to sample distances
dists = dist(t(exprs(vdsFull)))
mat = as.matrix(dists)
rownames(mat) = colnames(mat) = with(pData(cdsBlind), paste(design$condition, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

# case view
plotPCA(vdsFull, intgroup=c("condition"))
# replicate view
plotPCA(vdsFull, intgroup=c("condition", "libType"))