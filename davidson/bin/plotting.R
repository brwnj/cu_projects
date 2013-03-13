# in out title
args <- commandArgs(TRUE)

filein <- args[1]
fileout <- args[2]
title <- args[3]

dt <- read.table(args[1],
                header=TRUE,
                sep="\t",
                row.names=1)
pdf(fileout, width=8, height=11)
# fix bottom margin
par(mar=c(4,10,4,2) + 0.1)
barplot(dt$coverage,
        names.arg=rownames(dt),
        las=1,
        xlab="Count",
        main=title,
        horiz=TRUE)
dev.off()

library(ggplot2)
dt = read.table("~/projects/davidson/data/1/1.metadata", header=TRUE, row.names=1)
# bar plot
q <- qplot(dt$v_region)
# label rotation
q + theme(axis.text.x = element_text(angle = 90, hjust = 1))

q <- qplot(dt$v_region, dt$coverage)
q + theme(axis.text.x = element_text(angle = 90, hjust = 1))

p <- ggplot(dt, aes(v_region))
p <- p + geom_histogram(binwidth=2, fill="steelblue")
p <- p + geom_point(aes(v_region, coverage), alpha=I(1/50))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

d <- ggplot(dt, aes(reads))
d + geom_histogram()
