#Rscript figures.R ~/project_dir/cuffdiff_dir
args = commandArgs(trailingOnly=TRUE)
working_dir=args[1]
rm(args)

library(cummeRbund)

setwd(working_dir)
cuff = readCufflinks()

png("gene_dispersion.png")
dispersionPlot(genes(cuff))
dev.off()

png("fpkm_distribution_density.png")
csDensity(genes(cuff), replicates=TRUE)
dev.off()

png("fpkm_distribution_box.png")
csBoxplot(genes(cuff), replicates=TRUE)
dev.off()

png("sample_dendrogram.png")
csDendro(genes(cuff), replicates=TRUE)
dev.off()