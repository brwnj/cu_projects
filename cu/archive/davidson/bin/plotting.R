# in out title
# args <- commandArgs(TRUE)
# 
# filein <- args[1]
# fileout <- args[2]
# title <- args[3]
# 
# dt <- read.table(args[1],
#                 header=TRUE,
#                 sep="\t",
#                 row.names=1)
# pdf(fileout, width=8, height=11)
# # fix bottom margin
# par(mar=c(4,10,4,2) + 0.1)
# barplot(dt$coverage,
#         names.arg=rownames(dt),
#         las=1,
#         xlab="Count",
#         main=title,
#         horiz=TRUE)
# dev.off()

library(ggplot2)
setwd("~/projects/davidson/data/20130314/")
dte = read.table("246_metadata.txt", header=TRUE)
dto = read.table("135_metadata.txt", header=TRUE)

# plots for v region
p <- ggplot(dto, aes(v_region, fill=sample))
p <- p + geom_histogram(binwidth=2, position=position_dodge(), color="black")
# p <- p + facet_wrap(~ sample, ncol=1)
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 65, hjust = 1))
p <- p + labs(x="V Region", y="Counts", title="Identified Alpha Variable TCR Sequence Distribution")
# p <- p + theme(legend.position="none")
p <- p + theme(panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())
p <- p + scale_fill_discrete(name="Samples", labels=c("Alpha 1", "Alpha 3", "Alpha 5"))
p

# average coverage
ac <- ggplot(dt, aes(v_region, avg_coverage, color=sample))
ac <- ac + geom_boxplot(alpha=0.3, position=position_dodge(), outlier.colour=NA)
ac <- ac + theme_bw()
ac <- ac + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ac

# percentage of reads per contig
p <- ggplot(dt, aes(x=v_region, fill=sample))
p <- p + geom_histogram(binwidth=1, aes(weight=percent_total_reads), position=position_dodge())
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

# plots for j region
# bar plot
p <- ggplot(dte, aes(j_region, fill=sample))
p <- p + geom_histogram(binwidth=2, position=position_dodge())
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p