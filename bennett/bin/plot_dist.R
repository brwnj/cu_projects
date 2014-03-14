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
p <- p + geom_bar(aes(y = (..count..)/sum(..count..)))
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 65, hjust = 1))
p <- p + labs(x="V-Leader", y="Percent", title=paste(sample, "Leader", sep=" "))
p <- p + theme(legend.position="none")
p <- p + theme(panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())

q <- ggplot(t, aes(v_gene, fill=vleader))
q <- q + geom_bar(aes(y = (..count..)/sum(..count..)))
q <- q + theme_bw()
q <- q + theme(axis.text.x = element_text(angle = 65, hjust = 1))
q <- q + labs(x="IMGT V-Gene", y="Counts", title=paste(sample, "V-Gene", sep=" "))
q <- q + theme(legend.position="none")
q <- q + theme(panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())

r <- ggplot(t, aes(immunoglobulin, fill=immunoglobulin))
r <- r + geom_bar(aes(y = (..count..)/sum(..count..)))
r <- r + theme_bw()
r <- r + theme(axis.text.x = element_text(angle = 65, hjust = 1))
r <- r + labs(x="Immunoglobulin", y="Percent", title=paste(sample, "Immunoglobulin", sep=" "))
r <- r + theme(legend.position="none")
r <- r + theme(panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank())

pdf(out_pdf)
grid.arrange(p, q, r, nrow=3)
dev.off()
