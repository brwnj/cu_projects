# install.packages('corrplot')
# http://rstudio-pubs-static.s3.amazonaws.com/2107_4eb1adc1e4d44b93b6fde7eb801519fe.html

#normalized counts table
library(DESeq)
# python ../bin/data2matrix.py *MP*mirna_abundance.bed.gz > abundance_by_case.txt
ct = read.table("~/projects/hits-clip/data/abundance_by_case.txt", header=TRUE, row.names=1)
cds = newCountDataSet(ct, conditions=c(rep("t",8),rep("c",8)))
cds = estimateSizeFactors(cds)
normcounts = counts(cds, normalized=TRUE)

# correlation table
mat = cor(normcounts)
write.csv(mat, file="~/projects/hits-clip/data/correlation.csv")
write.csv((mat)^2, file="~/projects/hits-clip/data/correlation_coefficient.csv")

# plot
library(corrplot)
corrplot(mat, method="circle", order="alphabet", cl.lim=c(0,1))
corrRect(c(3,3,3,2,5), col="red", lwd=2)

ct = read.table("~/projects/hits-clip/data/20130829/pk11_pk24.txt", header=TRUE, row.names=1)
cds = newCountDataSet(ct, conditions=c(rep("t",1),rep("c",1)))
cds = estimateSizeFactors(cds)
normcounts = counts(cds, normalized=TRUE)
plot(log2(normcounts), cex=0.5)

library(ggplot2)
library(reshape2)
library(gridExtra)
pk11_24 = ggplot(log2(ct), aes(PK11, PK24)) + geom_point(size=5, family='arial') + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) + coord_cartesian(ylim = c(0,24), xlim = c(0,24))
pk11_24

ct = read.table("~/projects/hits-clip/data/20130829/peter_request.normed.txt", header=TRUE, row.names=1)
ct = log2(ct)
mt = melt(ct, id.vars=c("MCF7", "MCF7estr", "MDA231", "BT474"))

plot1 = ggplot(ct) + geom_point(aes(MCF7estr, MCF7), size=5, family = 'arial') + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) + coord_cartesian(ylim = c(0,24), xlim = c(0,24))
plot2 = ggplot(ct) + geom_point(aes(BT474, MCF7), size=5) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) + coord_cartesian(ylim = c(0,24), xlim = c(0,24))
plot3 = ggplot(ct) + geom_point(aes(MDA231, MCF7), size=5) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) + coord_cartesian(ylim = c(0,24), xlim = c(0,24))
grid.arrange(plot1, plot2, plot3, ncol=3)

ct = read.table("~/projects/hits-clip/data/20130913/mcf7_triplicate.normed.txt", header=TRUE, row.names=1)
cor(ct)
mx = cor(ct)
corrplot(mx, method="circle", order="alphabet", cl.lim=c(0,1))
ct <- log2(ct)
plot1 = ggplot(ct) + geom_point(aes(PK31, PK11), size=5, family = 'arial') + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) + coord_cartesian(ylim = c(0,24), xlim = c(0,24))
plot2 = ggplot(ct) + geom_point(aes(PK51, PK11), size=5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) + coord_cartesian(ylim = c(0,24), xlim = c(0,24))
plot3 = ggplot(ct) + geom_point(aes(PK51, PK31), size=5) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) + coord_cartesian(ylim = c(0,24), xlim = c(0,24))
grid.arrange(plot1, plot2, plot3, ncol=3)
