library("DESeq2")
library("biomaRt")

# need to subset sampleInfo and countData
samples <- data.frame(row.names=c("MCF7_0_fraction1","MCF7_0_fraction2","MCF7_totalRNA", "MCF7_6_fraction1", "MCF7_6_fraction2", "MCF7_6_totalRNA"),
                      condition=c("fraction", "fraction", "mRNA", "fraction", "fraction", "mRNA"),
                      time=c("0h", "0h", "0h", "6h", "6h", "6h"))
counts <- read.table("~/projects/kabos_ribosome/counts_deseq.txt", header=TRUE, row.names=1)
# remove unused columns
drops <- c("PK12_0","PK12_6")
counts <- counts[,!(names(counts) %in% drops)]

dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = samples,
    design = ~ time)

subdds <- dds[, dds$condition == "mRNA" ]
subdds <- DESeq(subdds)
subres <- results(subdds)

pdf("~/projects/kabos_ribosome/mRNA_0_to_6_pvaldist.pdf",width=7,height=5)
hist(subres$pvalue, breaks=20, col="steelblue", main="pval dist for mRNA 0hr to 6hr", xlab="pvalue")
dev.off()

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
write.csv(as.data.frame(subres), file="~/projects/kabos_ribosome/mRNA_0_to_6_DESeq.csv")

rld <- rlogTransformation(subdds)
# par(mfrow = c(1, 2))
# plot(log2(1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3)
pdf("~/projects/kabos_ribosome/mRNA_0_to_6_correlation.pdf",width=7,height=5)
plot(assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3, main="")
dev.off()

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
