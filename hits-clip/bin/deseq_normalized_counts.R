library(DESeq)
ct = read.table("~/projects/hits-clip/data/manoj_top10.csv", header=T, row.names=1, sep=",")
cds = newCountDataSet(ct, conditions=c(rep("t",8),rep("c",8)))
cds = estimateSizeFactors(cds)
write.csv(log(1 + counts( cds, normalized=TRUE )), file="~/projects/hits-clip/data/manoj_top10_norm.csv")

ct = read.table("~/projects/hits-clip/data/manoj_top25.csv", header=T, row.names=1, sep=",")
cds = newCountDataSet(ct, conditions=c(rep("t",8),rep("c",8)))
cds = estimateSizeFactors(cds)
write.csv(log(1 + counts( cds, normalized=TRUE )), file="~/projects/hits-clip/data/manoj_top25_norm.csv")

ct = read.table("~/projects/hits-clip/data/manoj_top50.csv", header=T, row.names=1, sep=",")
cds = newCountDataSet(ct, conditions=c(rep("t",8),rep("c",8)))
cds = estimateSizeFactors(cds)
write.csv(log(1 + counts( cds, normalized=TRUE )), file="~/projects/hits-clip/data/manoj_top50_norm.csv")

ct = read.table("~/projects/hits-clip/data/manoj_top100.csv", header=T, row.names=1, sep=",")
cds = newCountDataSet(ct, conditions=c(rep("t",8),rep("c",8)))
cds = estimateSizeFactors(cds)
write.csv(log(1 + counts( cds, normalized=TRUE )), file="~/projects/hits-clip/data/manoj_top100_norm.csv")

ct = read.table("~/projects/hits-clip/data/manoj_top200.csv", header=T, row.names=1, sep=",")
cds = newCountDataSet(ct, conditions=c(rep("t",8),rep("c",8)))
cds = estimateSizeFactors(cds)
write.csv(log(1 + counts( cds, normalized=TRUE )), file="~/projects/hits-clip/data/manoj_top200_norm.csv")

ct = read.table("~/projects/hits-clip/data/peter1_top50.csv", header=T, row.names=1, sep=",")
cds = newCountDataSet(ct, conditions=c(rep("t",4),rep("c",5)))
cds = estimateSizeFactors(cds)
write.csv(log(1 + counts(cds, normalized=TRUE)), file="~/projects/hits-clip/data/peter1_top50_norm.csv")

ct = read.table("~/projects/hits-clip/data/peter2_top50.csv", header=T, row.names=1, sep=",")
cds = newCountDataSet(ct, conditions=c(rep("t",5),rep("c",5)))
cds = estimateSizeFactors(cds)
write.csv(log(1 + counts(cds, normalized=TRUE)), file="~/projects/hits-clip/data/peter2_top50_norm.csv")

ct = read.table("~/projects/hits-clip/data/manoj_allsno.csv", header=T, row.names=1, sep=",")
cds = newCountDataSet(ct, conditions=c(rep("t",9),rep("c",8)))
cds = estimateSizeFactors(cds)
sizeFactors(cds)

write.csv(log(1 + counts(cds, normalized=TRUE)), file="~/projects/hits-clip/data/manoj_allsno_norm.csv")