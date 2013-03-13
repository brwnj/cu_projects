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
dt = read.table("~/projects/davidson/data/246_metadata.txt", header=TRUE)

# bar plot
p <- ggplot(dt, aes(v_region, fill=sample))
# add the histogram
p <- p + geom_histogram(binwidth=2, position=position_dodge())
# group plots by sample
# p <- p + facet_grid(. ~ sample)
# white background
p <- p + theme_bw()
# rotate the labels
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p