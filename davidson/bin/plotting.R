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