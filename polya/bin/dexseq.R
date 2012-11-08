library(DEXSeq)
library(DESeq)
setwd("~/projects/polya/data")

MPdesign1pos = data.frame(row.names = c("MP51pos",
                                        "MP51posx",
                                        "MP52pos",
                                        "MP52posx"),
                          condition = c("control",
                                        "control",
                                        "treatment",
                                        "treatment"))
ecsMP1pos = read.HTSeqCounts(countfiles = c("MP51.umi_filtered.pos.counts",
                                            "MP51.umi_filtered.posx.counts",
                                            "MP52.umi_filtered.pos.counts",
                                            "MP52.umi_filtered.posx.counts"), 
                             design = MPdesign1pos,
                             flattenedfile = "polya.flattened.gff")
ecsMP1pos = estimateSizeFactors(ecsMP1pos)
ecsMP1pos = estimateDispersions(ecsMP1pos,
                                minCount = 10)
ecsMP1pos = fitDispersionFunction(ecsMP1pos)
ecsMP1pos = testForDEU(ecsMP1pos)
ecsMP1pos = estimatelog2FoldChanges(ecsMP1pos)
resMP1pos = DEUresultTable(ecsMP1pos)
write.csv(resMP1pos, file="MP51pos_MP52pos_dexseq.csv")

MPdesign1neg = data.frame(row.names = c("MP51neg", "MP51negx", "MP52neg", "MP52negx"),
                          condition = c("control", "control", "treatment", "treatment"))
ecsMP1neg = read.HTSeqCounts(countfiles = c("MP51.umi_filtered.neg.counts", "MP51.umi_filtered.negx.counts", "MP52.umi_filtered.neg.counts", "MP52.umi_filtered.negx.counts"), 
                             design = MPdesign1neg)
ecsMP1neg = estimateSizeFactors(ecsMP1neg)
sizeFactors(ecsMP1neg)
ecsMP1neg = estimateDispersions(ecsMP1neg, minCount = 10)
ecsMP1neg = fitDispersionFunction(ecsMP1neg)
ecsMP1neg = testForDEU(ecsMP1neg)
ecsMP1neg = estimatelog2FoldChanges(ecsMP1neg)
resMP1neg = DEUresultTable(ecsMP1neg)
write.csv(resMP1neg, file="MP51neg_MP52neg_dexseq.csv")


MPdesign2pos = data.frame(row.names = c("MP51pos", "MP51posx", "MP53pos", "MP53posx"),
                          condition = c("control", "control", "treatment", "treatment"))

MPdesign2neg = data.frame(row.names = c("MP51neg", "MP51negx", "MP53neg", "MP53negx"),
                          condition = c("control", "control", "treatment", "treatment"))

PKdesignpos = data.frame(row.names = c("PK61pos","PK61posx","PK62pos","PK62posx"),
                         condition = c("control","control","treatment","treatment"))

PKdesignneg = data.frame(row.names = c("PK61neg","PK61negx","PK62neg","PK62negx"),
                         condition = c("control","control","treatment","treatment"))



ecsManoj2 = read.HTSeqCounts(countfiles = c("test51.counts", "test51x.counts", "test53.counts", "test53x.counts"), 
                             design = manojDesign2,
                             flattenedfile = "dexseq.anno.unsorted.gff")

ecsPeter = read.HTSeqCounts(countfiles = c("PK61.counts", "PK61x.counts", "PK62.counts", "PK62x.counts"), 
                            design = peterDesign, 
                            flattenedfile = "dexseq.anno.unsorted.gff")



ecsManoj2 = estimateSizeFactors(ecsManoj2)
ecsPeter = estimateSizeFactors(ecsPeter)


ecsManoj2 = estimateDispersions(ecsManoj2, minCount = 3)
ecsPeter = estimateDispersions(ecsPeter, minCount = 3)


ecsManoj2 = fitDispersionFunction(ecsManoj2)
ecsPeter = fitDispersionFunction(ecsPeter)


ecsManoj2 = testForDEU(ecsManoj2)
ecsPeter = testForDEU(ecsPeter)


ecsManoj2 = estimeatelog2FoldChanges(ecsManoj2)
ecsPeter = estimatelog2FoldChanges(ecsPeter)


manojres2 = DEUresultTable(ecsManoj2)
peterres = DEUresultTable(ecsPeter)

table(manojres1$padjust<0.1)
table(peterres$padjust<0.1)


write.csv(peterres, file="PK61_PK62.dexseq.csv")

DEXSeqHTML(ecsManoj1, path="MP51_MP52")

plotDEXSeq(ecsManoj1, "KMO", displayTranscripts=T, cex.axis=1.2, cex=1.3, lwd=2, legend=T)
plotDEXSeq(ecsManoj1, "TNFAIP1", displayTranscripts=T, cex.axis=1.2, cex=1.3, lwd=2, legend=T)