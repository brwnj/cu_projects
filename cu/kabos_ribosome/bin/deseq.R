library("DESeq2")
library("biomaRt")

# need to subset sampleInfo and countData
samples <- data.frame(row.names=c("MCF7_0_fraction1","MCF7_0_fraction2","MCF7_6_fraction1","MCF7_6_fraction2"),
                      condition=c("fraction","fraction","fraction","fraction"),
                      time=c("0h","0h","6h","6h"))
# samples <- data.frame(row.names=c("MCF7_0_fraction1","MCF7_0_fraction2","MCF7_0_totalRNA","MCF7_6_fraction1","MCF7_6_fraction2","MCF7_6_totalRNA"),
#                       condition=c("fraction", "fraction", "mRNA", "fraction", "fraction", "mRNA"),
#                       time=c("0h", "0h", "0h", "6h", "6h", "6h"))

counts <- read.table("~/projects/kabos_ribosome/data/counts_deseq.txt", header=TRUE, row.names=1)
# remove unused columns
drops <- c("PK12_0","PK12_6","MCF7_0_totalRNA","MCF7_6_totalRNA")
counts <- counts[,!(names(counts) %in% drops)]

# counts["ENSG00000082175",]


# change the design
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = samples,
    design = ~ time)

# subdds <- dds[, dds$condition == "mRNA" ]
# subdds <- DESeq(subdds, test="Wald")
# subres <- results(subdds)

dds <- DESeq(dds)
res <- results(dds)
res["ENSG00000082175",]


# mRNA
samples <- data.frame(row.names=c("MCF7_0_totalRNA","MCF7_6_totalRNA"),
                      condition=c("mRNA","mRNA"),
                      time=c("0h","6h"))
counts <- read.table("~/projects/kabos_ribosome/data/counts_deseq.txt", header=TRUE, row.names=1)
# remove unused columns
drops <- c("PK12_0","PK12_6","MCF7_0_fraction1","MCF7_0_fraction2","MCF7_6_fraction1","MCF7_6_fraction2")
counts <- counts[,!(names(counts) %in% drops)]
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = samples,
    design = ~ time)
dds <- DESeq(dds)
res <- results(dds)
res["ENSG00000082175",]


pdf("~/projects/kabos_ribosome/mRNA_0_to_6_pvaldist.pdf",width=7,height=5)
hist(subres$pvalue, breaks=20, col="steelblue", main="pval dist for mRNA 0hr to 6hr", xlab="pvalue")
dev.off()

# add gene names to results
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
par(mfrow = c(1, 2))

# alternate transform
t = log2(1+counts(subdds, normalized=TRUE))
write.csv(t, file="~/projects/kabos_ribosome/log2_mRNA_0_to_6.csv")
t = assay(rld)
write.csv(t, file="~/projects/kabos_ribosome/rld_mRNA_0_to_6.csv")

plot(log2(1+counts(subdds, normalized=TRUE)), col="#00000020", pch=20, cex=0.3, ylab="0h", xlab="6h")
text(log2(1+counts(subdds, normalized=TRUE)), labels=row.names(subdds), pos=4, cex=0.3)

pdf("~/projects/kabos_ribosome/mRNA_0_to_6_correlation.pdf",width=7,height=5)
# appears much more normalized
plot(assay(rld), col="#00000020", pch=20, cex=0.3, main="", ylab="0h", xlab="6h")
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

library(edgeR)
counts <- read.table("~/projects/kabos_ribosome/data/counts_deseq.txt", header=TRUE, row.names=1)
# drops <- c("PK12_0","PK12_6","MCF7_0_totalRNA","MCF7_6_totalRNA")
drops <- c("PK12_0","PK12_6","MCF7_0_fraction1","MCF7_0_fraction2","MCF7_6_fraction1","MCF7_6_fraction2")
counts <- counts[,!(names(counts) %in% drops)]
# group <- factor(c("0hr","0hr","6hr","6hr"))
group <- factor(c("0hr","6hr"))
ds <- DGEList(counts=counts, group=group, genes=rownames(counts))
# ds <- calcNormFactors(ds)
# ds <- estimateCommonDisp(ds)
# ds <- estimateTagwiseDisp(ds)
# et <- exactTest(ds)

bcv <- 0.2
et <- exactTest(ds, dispersion=bcv^2)

topTags(et)

# write.csv(topTags(et, n=60000), file="~/projects/kabos_ribosome/fraction_0_to_6_edgeR.csv")
write.csv(topTags(et, n=60000), file="~/projects/kabos_ribosome/mRNA_0_to_6_edgeR.csv")

et$ensembl <- sapply(strsplit(rownames(et), split="nn+"), "[", 1)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genemap <- getBM(attributes=c("ensembl_gene_id", "entrezgene", "hgnc_symbol"),
                 filters="ensembl_gene_id",
                 values=et$ensembl,
                 mart=ensembl)
idx <- match(et$ensembl, genemap$ensembl_gene_id)
et$entrez <- genemap$entrezgene[idx]
et$hgnc_symbol <- genemap$hgnc_symbol[idx]

