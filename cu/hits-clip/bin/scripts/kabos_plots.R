library(ggplot2)
library(reshape2)
setwd("~/projects/hits-clip/data/20130404/")
df <- read.table("GSE22219_values.txt", header=TRUE, sep="\t")
df <- subset(df, select=c(Gene_notest, hsa.miR.9.5p, ESR1, AR, RARA))
m <- melt(df, id.vars=c("Gene_notest", "hsa.miR.9.5p")) 
p <- ggplot(m, aes(variable, value)) + 
    geom_boxplot(aes(fill=factor(hsa.miR.9.5p)), m) +
    theme_bw() +
    labs(x="Gene", y="Intensity") +
    scale_fill_manual(values=c("blue", "red"), name="hsa-miR-9-5p", labels=c("Low", "High"))
p