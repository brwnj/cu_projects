library("DESeq2")
library("biomaRt")

# need to subset sampleInfo and countData
samples <- data.frame(row.names=c("MCF7_0_fraction1","MCF7_0_fraction2","MCF7_totalRNA", "MCF7_6_fraction1", "MCF7_6_fraction2", "MCF7_6_totalRNA", "PK12_0", "PK12_6"),
                      condition=c("fraction", "fraction", "mRNA", "fraction", "fraction", "mRNA", "fraction", "fraction"),
                      time=c("0h", "0h", "0h", "6h", "6h", "6h", "100h", "100h"))
counts <- read.table("~/projects/kabos_ribosome/data/20140403/counts.txt", header=TRUE, row.names=1)

dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = samples,
    design = ~ condition)

subdds <- dds[, dds$time == "6h" ]
subdds <- DESeq(subdds)
subres <- results(subdds)

hist(subres$pvalue, breaks=20, col="steelblue", main="pval dist for MCF7_6", xlab="pvalue")

# add gene names
subres$ensembl <- sapply(strsplit(rownames(subres), split="nn+"), "[", 1)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genemap <- getBM(attributes=c("ensembl_gene_id", "entrezgene", "hgnc_symbol"),
                  filters="ensembl_gene_id",
                  values=subres$ensembl,
                  mart=ensembl)
idx <- match(subres$ensembl, genemap$ensembl_gene_id)
subres$entrez <- genemap$entrezgene[idx]
subres$hgnc_symbol <- genemap$hgnc_symbol[idx]
write.csv(as.data.frame(subres), file="MCF7-6_DESeq.csv")

rld <- rlogTransformation(subdds)
# par(mfrow = c(1, 2))
# plot(log2(1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3)
plot(assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3, main="")

# sampleDists <- dist( t( assay(rld) ) )
# sampleDistMatrix <- as.matrix( sampleDists )
# rownames(sampleDistMatrix) <- paste( rld$treatment,
#                                      rld$patient, sep="-" )
# colnames(sampleDistMatrix) <- NULL
# library( "gplots" )
# library( "RColorBrewer" )
# colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# heatmap.2( sampleDistMatrix, trace="none", col=colours)

# ramp <- 1:3/3
# cols <- c(rgb(ramp, 0, 0),
#           rgb(0, ramp, 0))
# print(plotPCA( rld, intgroup = c( "patient", "treatment"), col=cols ))

# library( "genefilter" )
# topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
# heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
#            trace="none", dendrogram="column",
#            col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
