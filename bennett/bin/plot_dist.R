usage <- "Usage: Rscript aa_dist_table sample_name out_pdf"

args <- commandArgs(TRUE)
if (length(args) == 0)
    stop(usage)

library(ggplot2)
library(gridExtra)

in_table = args[1]
sample = args[2]
out_pdf = args[3]

t = read.table(in_table, row.names=1)

p <- ggplot(t, aes(vleader, fill=vleader))
p <- p + geom_histogram(binwidth=2, position=position_dodge(), color="black")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 65, hjust = 1))
p <- p + labs(x="V-Leader", y="Counts", title=paste(sample, "Unique CDR3 Distribution", sep=" "))
p <- p + theme(legend.position="none")
p <- p + theme(panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())

q <- ggplot(t, aes(v_gene, fill=vleader))
q <- q + geom_histogram(binwidth=2, position=position_dodge(), color="black")
q <- q + theme_bw()
q <- q + theme(axis.text.x = element_text(angle = 65, hjust = 1))
q <- q + labs(x="IMGT V-Gene", y="Counts")
q <- q + theme(legend.position="none")
q <- q + theme(panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())

pdf(out_pdf)
grid.arrange(p, q, nrow=2)
dev.off()
