setwd("~/projects/hits-clip/data/")

library(gplots)
library(RColorBrewer)
MP = as.matrix(read.delim("~/projects/hits-clip/data/manoj_top50.csv", header=T, row.names=1, sep=","))
PK1 = as.matrix(read.delim("~/projects/hits-clip/data/peter1_top50.csv", header=T, row.names=1, sep=","))
PK2 = as.matrix(read.delim("~/projects/hits-clip/data/peter2_top50.csv", header=T, row.names=1, sep=","))

MP = MP - rowMeans(MP)
PK1 = PK1 - rowMeans(PK1)
PK2 = PK2 - rowMeans(PK2)

png("~/projects/hits-clip/data/manoj_R_top50.png", width=2000, height=2000)
cmap = colorRampPalette(c("lightyellow4", "yellow", "blue", "royalblue"))(128)
print(ncol(MP))
heatmap.2(MP,
          cexRow=2,
          Rowv=TRUE,
          Colv=TRUE,
          labRow=NA,
          labCol=colnames(MP),
          margins=c(5, 5),
          col=cmap,
          trace="none",
          #ColSideColors=c("gray", "slategray")[c(2, 1, 2, 1, 2, 1, 2)],
          density.info='none',
          key=FALSE
)
dev.off()
png("~/projects/hits-clip/data/peter_set1_R_top50.png", width=2000, height=2000)
heatmap.2(PK1,
          cexRow=2,
          Rowv=TRUE,
          Colv=TRUE,
          labRow=NA,
          labCol=colnames(PK1),
          margins=c(5, 5),
          col=cmap,
          trace="none",
          #ColSideColors=c("gray", "slategray")[c(2, 1, 2, 1, 2, 1, 2)],
          density.info='none',
          key=FALSE
)
dev.off()
png("~/projects/hits-clip/data/peter_set2_R_top50.png", width=2000, height=2000)
heatmap.2(PK2,
          cexRow=2,
          Rowv=TRUE,
          Colv=TRUE,
          labRow=NA,
          labCol=colnames(PK1),
          margins=c(3, 3),
          col=cmap,
          trace="none",
          #ColSideColors=c("gray", "slategray")[c(2, 1, 2, 1, 2, 1, 2)],
          density.info='none',
          key=FALSE
 )
dev.off()