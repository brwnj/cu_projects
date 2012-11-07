library(DEXSeq)
library(DESeq)
setwd("~/projects/polya/data")

manojDesign1 = data.frame(
    row.names = c("MP51", "MP51x", "MP52", "MP52x"),
    condition = c("control", "control", "treatment", "treatment"))

manojDesign2 = data.frame(
    row.names = c("MP51", "MP51x", "MP53", "MP53x"),
    condition = c("control", "control", "treatment", "treatment"))

peterDesign = data.frame(
    row.names = c("PK61","PK61x","PK62","PK62x"),
    condition = c("control","control","treatment","treatment"))

ecsManoj1 = read.HTSeqCounts(countfiles = c("MP51.counts", "MP51x.counts", "MP52.counts", "MP52x.counts"), 
                        design = manojDesign1,
                        flattenedfile = "dexseq.anno.unsorted.gff")

ecsManoj2 = read.HTSeqCounts(countfiles = c("test51.counts", "test51x.counts", "test53.counts", "test53x.counts"), 
                             design = manojDesign2,
                             flattenedfile = "dexseq.anno.unsorted.gff")

ecsPeter = read.HTSeqCounts(countfiles = c("PK61.counts", "PK61x.counts", "PK62.counts", "PK62x.counts"), 
                            design = peterDesign, 
                            flattenedfile = "dexseq.anno.unsorted.gff")

ecsManoj1 = estimateSizeFactors(ecsManoj1)
ecsManoj2 = estimateSizeFactors(ecsManoj2)
ecsPeter = estimateSizeFactors(ecsPeter)

ecsManoj1 = estimateDispersions(ecsManoj1, minCount = 3)
ecsManoj2 = estimateDispersions(ecsManoj2, minCount = 3)
ecsPeter = estimateDispersions(ecsPeter, minCount = 3)

ecsManoj1 = fitDispersionFunction(ecsManoj1)
ecsManoj2 = fitDispersionFunction(ecsManoj2)
ecsPeter = fitDispersionFunction(ecsPeter)

ecsManoj1 = testForDEU(ecsManoj1)
ecsManoj2 = testForDEU(ecsManoj2)
ecsPeter = testForDEU(ecsPeter)

ecsManoj1 = estimatelog2FoldChanges(ecsManoj1)
ecsManoj2 = estimeatelog2FoldChanges(ecsManoj2)
ecsPeter = estimatelog2FoldChanges(ecsPeter)

manojres1 = DEUresultTable(ecsManoj1)
manojres2 = DEUresultTable(ecsManoj2)
peterres = DEUresultTable(ecsPeter)

table(manojres1$padjust<0.1)
table(peterres$padjust<0.1)

write.csv(manojres1, file="MP51_MP52.dexseq.csv")
write.csv(peterres, file="PK61_PK62.dexseq.csv")

DEXSeqHTML(ecsManoj1, path="MP51_MP52")

plotDEXSeq(ecsManoj1, "KMO", displayTranscripts=T, cex.axis=1.2, cex=1.3, lwd=2, legend=T)
plotDEXSeq(ecsManoj1, "TNFAIP1", displayTranscripts=T, cex.axis=1.2, cex=1.3, lwd=2, legend=T)