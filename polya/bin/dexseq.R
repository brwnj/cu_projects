library(DEXSeq)
# library(DESeq)
setwd("~/projects/polya/data")

annotation_file = "~/projects/polya/data/polya_extended.gff"

# 1 to 2 positive strand
mp1 = data.frame(
    row.names = c("MP51.umi_filtered.counts","MP51.umi_filteredx.counts",
                  "MP52.umi_filtered.counts","MP52.umi_filteredx.counts"),
    condition = c(rep("control",2),
                  rep("treatment",2)),
    stringsAsFactors = TRUE,
    check.names = FALSE
    )
ecs = read.HTSeqCounts(
    design = mp1,
    countfiles = c("MP51.umi_filtered.counts",
                   "MP51.umi_filteredx.counts",
                   "MP52.umi_filtered.counts",
                   "MP52.umi_filteredx.counts")
    )
ecsMP1pos = estimateSizeFactors(ecsMP1pos)
ecsMP1pos = estimateDispersions(ecsMP1pos, minCount = 10)
ecsMP1pos = fitDispersionFunction(ecsMP1pos)
ecsMP1pos = testForDEU(ecsMP1pos)
ecsMP1pos = estimatelog2FoldChanges(ecsMP1pos)
resMP1pos = DEUresultTable(ecsMP1pos)
write.csv(resMP1pos, file="MP51pos_MP52pos_dexseq.csv")

# 1 to 2 negative strand
MPdesign1neg = data.frame(row.names = c("MP51neg","MP51negx",
                                        "MP52neg","MP52negx"),
                          condition = c(rep("control", 2),
                                        rep("treatment", 2)))
ecsMP1neg = read.HTSeqCounts(countfiles = c("MP51.umi_filtered.neg.counts",
                                            "MP51.umi_filtered.negx.counts",
                                            "MP52.umi_filtered.neg.counts",
                                            "MP52.umi_filtered.negx.counts"), 
                             design = MPdesign1neg)
ecsMP1neg = estimateSizeFactors(ecsMP1neg)
# sizeFactors(ecsMP1neg)
ecsMP1neg = estimateDispersions(ecsMP1neg, minCount = 10)
ecsMP1neg = fitDispersionFunction(ecsMP1neg)
ecsMP1neg = testForDEU(ecsMP1neg)
ecsMP1neg = estimatelog2FoldChanges(ecsMP1neg)
resMP1neg = DEUresultTable(ecsMP1neg)
write.csv(resMP1neg, file="MP51neg_MP52neg_dexseq.csv")

# 1 to 3 positive
MPdesign2pos = data.frame(row.names = c("MP51pos","MP51posx",
                                        "MP53pos","MP53posx"),
                          condition = c("control","control",
                                        "treatment","treatment"))
ecsMP2pos = read.HTSeqCounts(countfiles = c("MP51.umi_filtered.pos.counts",
                                            "MP51.umi_filtered.posx.counts",
                                            "MP53.umi_filtered.pos.counts",
                                            "MP53.umi_filtered.posx.counts"), 
                             design = MPdesign2pos)
ecsMP2pos = estimateSizeFactors(ecsMP2pos)
ecsMP2pos = estimateDispersions(ecsMP2pos, minCount = 10)
ecsMP2pos = fitDispersionFunction(ecsMP2pos)
ecsMP2pos = testForDEU(ecsMP2pos)
ecsMP2pos = estimatelog2FoldChanges(ecsMP2pos)
resMP2pos = DEUresultTable(ecsMP2pos)
write.csv(resMP2pos, file="MP51pos_MP53pos_dexseq.csv")

# 1 to 3 negative
MPdesign2neg = data.frame(row.names = c("MP51neg","MP51negx",
                                        "MP53neg","MP53negx"),
                          condition = c("control","control",
                                        "treatment","treatment"))
ecsMP2neg = read.HTSeqCounts(countfiles = c("MP51.umi_filtered.neg.counts",
                                            "MP51.umi_filtered.negx.counts",
                                            "MP53.umi_filtered.neg.counts",
                                            "MP53.umi_filtered.negx.counts"), 
                             design = MPdesign2neg)
ecsMP2neg = estimateSizeFactors(ecsMP2neg)
# sizeFactors(ecsMP2neg)
ecsMP2neg = estimateDispersions(ecsMP2neg, minCount = 10)
ecsMP2neg = fitDispersionFunction(ecsMP2neg)
ecsMP2neg = testForDEU(ecsMP2neg)
ecsMP2neg = estimatelog2FoldChanges(ecsMP2neg)
resMP2neg = DEUresultTable(ecsMP2neg)
write.csv(resMP2neg, file="MP51neg_MP53neg_dexseq.csv")

# peter positive
PKdesignpos = data.frame(row.names = c("PK61pos","PK61posx",
                                       "PK62pos","PK62posx"),
                         condition = c("control","control",
                                       "treatment","treatment"))
ecsPKpos = read.HTSeqCounts(countfiles = c("PK61.umi_filtered.pos.counts",
                                           "PK61.umi_filtered.posx.counts",
                                           "PK62.umi_filtered.pos.counts",
                                           "PK62.umi_filtered.posx.counts"), 
                            design = PKdesignpos)
ecsPKpos = estimateSizeFactors(ecsPKpos)
ecsPKpos = estimateDispersions(ecsPKpos, minCount = 10)
ecsPKpos = fitDispersionFunction(ecsPKpos)
ecsPKpos = testForDEU(ecsPKpos)
ecsPKpos = estimatelog2FoldChanges(ecsPKpos)
resPKpos = DEUresultTable(ecsPKpos)
write.csv(resPKpos, file="PK61pos_PK62pos_dexseq.csv")

# peter negative
PKdesignneg = data.frame(row.names = c("PK61neg","PK61negx",
                                       "PK62neg","PK62negx"),
                         condition = c("control","control",
                                       "treatment","treatment"))
ecsPKneg = read.HTSeqCounts(countfiles = c("PK61.umi_filtered.neg.counts",
                                           "PK61.umi_filtered.negx.counts",
                                           "PK62.umi_filtered.neg.counts",
                                           "PK62.umi_filtered.negx.counts"), 
                            design = PKdesignneg)
ecsPKneg = estimateSizeFactors(ecsPKneg)
ecsPKneg = estimateDispersions(ecsPKneg, minCount = 10)
ecsPKneg = fitDispersionFunction(ecsPKneg)
ecsPKneg = testForDEU(ecsPKneg)
ecsPKneg = estimatelog2FoldChanges(ecsPKneg)
resPKneg = DEUresultTable(ecsPKneg)
write.csv(resPKneg, file="PK61neg_PK62neg_dexseq.csv")

# QC
table(res$padjust<0.1)

plotDEXSeq(ecsManoj1, "KMO", displayTranscripts=T, cex.axis=1.2, cex=1.3, lwd=2, legend=T)
plotDEXSeq(ecsManoj1, "TNFAIP1", displayTranscripts=T, cex.axis=1.2, cex=1.3, lwd=2, legend=T)