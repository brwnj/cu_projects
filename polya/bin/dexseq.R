source("http://bioconductor.org/biocLite.R")
biocLite("DEXSeq")
library(DEXSeq)
# library(DESeq)
setwd("~/projects/polya/data")

annotation_file = "~/projects/polya/data/polya_extended_with_names.gff"

# 1 to 2 positive strand
MPdesign1pos = as.matrix(data.frame(
    row.names = c("MP51","MP51x",
                  "MP52","MP52x"),
    condition = c(rep("control",2),
                  rep("treatment",2)),
    stringsAsFactors = TRUE,
    check.names = FALSE
    ))

ecsMP1pos = read.HTSeqCounts(
    design = MPdesign1pos,
    countfiles = c("MP51.pos.counts",
                   "MP51.posx.counts",
                   "MP52.pos.counts",
                   "MP52.posx.counts")
    )

ecsMP1pos = estimateSizeFactors(ecsMP1pos)
ecsMP1pos = estimateDispersions(ecsMP1pos, minCount = 1)
ecsMP1pos = fitDispersionFunction(ecsMP1pos)
ecsMP1pos = testForDEU(ecsMP1pos)
ecsMP1pos = estimatelog2FoldChanges(ecsMP1pos)
resMP1pos = DEUresultTable(ecsMP1pos)
write.csv(resMP1pos, file="MP51pos_MP52pos_dexseq_cutoff1.csv")

# 1 to 2 negative strand
MPdesign1neg = data.frame(row.names = c("MP51neg","MP51negx",
                                        "MP52neg","MP52negx"),
                          condition = c(rep("control", 2),
                                        rep("treatment", 2)))
ecsMP1neg = read.HTSeqCounts(countfiles = c("MP51.neg.counts",
                                            "MP51.negx.counts",
                                            "MP52.neg.counts",
                                            "MP52.negx.counts"), 
                             design = MPdesign1neg)
ecsMP1neg = estimateSizeFactors(ecsMP1neg)
# sizeFactors(ecsMP1neg)
ecsMP1neg = estimateDispersions(ecsMP1neg, minCount = 1)
ecsMP1neg = fitDispersionFunction(ecsMP1neg)
ecsMP1neg = testForDEU(ecsMP1neg)
ecsMP1neg = estimatelog2FoldChanges(ecsMP1neg)
resMP1neg = DEUresultTable(ecsMP1neg)
write.csv(resMP1neg, file="MP51neg_MP52neg_dexseq_cutoff1.csv")

# 1 to 3 positive
MPdesign2pos = data.frame(row.names = c("MP51pos","MP51posx",
                                        "MP53pos","MP53posx"),
                          condition = c("control","control",
                                        "treatment","treatment"))
ecsMP2pos = read.HTSeqCounts(countfiles = c("MP51.pos.counts",
                                            "MP51.posx.counts",
                                            "MP53.pos.counts",
                                            "MP53.posx.counts"), 
                             design = MPdesign2pos)
ecsMP2pos = estimateSizeFactors(ecsMP2pos)
ecsMP2pos = estimateDispersions(ecsMP2pos, minCount = 1)
ecsMP2pos = fitDispersionFunction(ecsMP2pos)
ecsMP2pos = testForDEU(ecsMP2pos)
ecsMP2pos = estimatelog2FoldChanges(ecsMP2pos)
resMP2pos = DEUresultTable(ecsMP2pos)
write.csv(resMP2pos, file="MP51pos_MP53pos_dexseq_cutoff1.csv")

# 1 to 3 negative
MPdesign2neg = data.frame(row.names = c("MP51neg","MP51negx",
                                        "MP53neg","MP53negx"),
                          condition = c("control","control",
                                        "treatment","treatment"))
ecsMP2neg = read.HTSeqCounts(countfiles = c("MP51.neg.counts",
                                            "MP51.negx.counts",
                                            "MP53.neg.counts",
                                            "MP53.negx.counts"), 
                             design = MPdesign2neg)
ecsMP2neg = estimateSizeFactors(ecsMP2neg)
# sizeFactors(ecsMP2neg)
ecsMP2neg = estimateDispersions(ecsMP2neg, minCount = 1)
ecsMP2neg = fitDispersionFunction(ecsMP2neg)
ecsMP2neg = testForDEU(ecsMP2neg)
ecsMP2neg = estimatelog2FoldChanges(ecsMP2neg)
resMP2neg = DEUresultTable(ecsMP2neg)
write.csv(resMP2neg, file="MP51neg_MP53neg_dexseq_cutoff1.csv")

# 2 to 3 positive strand
MPdesign3pos = as.matrix(data.frame(
    row.names = c("MP52","MP52x",
                  "MP53","MP53x"),
    condition = c(rep("control",2),
                  rep("treatment",2)),
    stringsAsFactors = TRUE,
    check.names = FALSE
))
ecsMP3pos = read.HTSeqCounts(
    design = MPdesign3pos,
    countfiles = c("MP52.pos.counts",
                   "MP52.posx.counts",
                   "MP53.pos.counts",
                   "MP53.posx.counts")
)
ecsMP3pos = estimateSizeFactors(ecsMP3pos)
ecsMP3pos = estimateDispersions(ecsMP3pos, minCount = 1)
ecsMP3pos = fitDispersionFunction(ecsMP3pos)
ecsMP3pos = testForDEU(ecsMP3pos)
ecsMP3pos = estimatelog2FoldChanges(ecsMP3pos)
resMP3pos = DEUresultTable(ecsMP3pos)
write.csv(resMP3pos, file="MP52pos_MP53pos_dexseq_cutoff1.csv")

# 2 to 3 negative
MPdesign3neg = data.frame(row.names = c("MP52neg","MP52negx",
                                        "MP53neg","MP53negx"),
                          condition = c("control","control",
                                        "treatment","treatment"))
ecsMP3neg = read.HTSeqCounts(countfiles = c("MP52.neg.counts",
                                            "MP52.negx.counts",
                                            "MP53.neg.counts",
                                            "MP53.negx.counts"), 
                             design = MPdesign3neg)
ecsMP3neg = estimateSizeFactors(ecsMP3neg)
# sizeFactors(ecsMP2neg)
ecsMP3neg = estimateDispersions(ecsMP3neg, minCount = 1)
ecsMP3neg = fitDispersionFunction(ecsMP3neg)
ecsMP3neg = testForDEU(ecsMP3neg)
ecsMP3neg = estimatelog2FoldChanges(ecsMP3neg)
resMP3neg = DEUresultTable(ecsMP3neg)
write.csv(resMP3neg, file="MP52neg_MP53neg_dexseq_cutoff1.csv")

ecsMP1pos = estimateSizeFactors(ecsMP1pos)
ecsMP1pos = estimateDispersions(ecsMP1pos, minCount = 1)
ecsMP1pos = fitDispersionFunction(ecsMP1pos)
ecsMP1pos = testForDEU(ecsMP1pos)
ecsMP1pos = estimatelog2FoldChanges(ecsMP1pos)
resMP1pos = DEUresultTable(ecsMP1pos)
write.csv(resMP1pos, file="MP51pos_MP52pos_dexseq_cutoff1.csv")


# peter positive
PKdesignpos = data.frame(row.names = c("PK61pos","PK61posx",
                                       "PK62pos","PK62posx"),
                         condition = c("control","control",
                                       "treatment","treatment"))
ecsPKpos = read.HTSeqCounts(countfiles = c("PK61.pos.counts",
                                           "PK61.posx.counts",
                                           "PK62.pos.counts",
                                           "PK62.posx.counts"), 
                            design = PKdesignpos)
ecsPKpos = estimateSizeFactors(ecsPKpos)
ecsPKpos = estimateDispersions(ecsPKpos, minCount = 1)
ecsPKpos = fitDispersionFunction(ecsPKpos)
ecsPKpos = testForDEU(ecsPKpos)
ecsPKpos = estimatelog2FoldChanges(ecsPKpos)
resPKpos = DEUresultTable(ecsPKpos)
write.csv(resPKpos, file="PK61pos_PK62pos_dexseq_cutoff1.csv")

# peter negative
PKdesignneg = data.frame(row.names = c("PK61neg","PK61negx",
                                       "PK62neg","PK62negx"),
                         condition = c("control","control",
                                       "treatment","treatment"))
ecsPKneg = read.HTSeqCounts(countfiles = c("PK61.neg.counts",
                                           "PK61.negx.counts",
                                           "PK62.neg.counts",
                                           "PK62.negx.counts"), 
                            design = PKdesignneg)
ecsPKneg = estimateSizeFactors(ecsPKneg)
ecsPKneg = estimateDispersions(ecsPKneg, minCount = 1)
ecsPKneg = fitDispersionFunction(ecsPKneg)
ecsPKneg = testForDEU(ecsPKneg)
ecsPKneg = estimatelog2FoldChanges(ecsPKneg)
resPKneg = DEUresultTable(ecsPKneg)
write.csv(resPKneg, file="PK61neg_PK62neg_dexseq_cutoff1.csv")

# QC
table(res$padjust<0.1)

plotDEXSeq(ecsManoj1, "KMO", displayTranscripts=T, cex.axis=1.2, cex=1.3, lwd=2, legend=T)
plotDEXSeq(ecsManoj1, "TNFAIP1", displayTranscripts=T, cex.axis=1.2, cex=1.3, lwd=2, legend=T)