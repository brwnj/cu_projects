library(ChIPpeakAnno)
table = read.table(
    "~/projects/artinger/data/2Som_chip1_GCCAAT_L006_R1_001_peaks.bed",
    header=FALSE)
bed = BED2RangedData(table)

# annotation
mart = useMart(biomart="ensembl", dataset="drerio_gene_ensembl")
tssanno = getAnnotation(mart, "TSS")

# annotate peaks
annotatedPeaks = annotatePeakInBatch(bed, AnnotationData=tssanno)
# write.csv(as.data.frame(annotatedPeaks), file="~/remote/projects/artinger/results/common/2Som_chip1_GCCAAT_L006_R1_001/2Som_chip1_GCCAAT_L006_R1_001_peaks.annotated.bed")

# enriched go terms
# source("http://bioconductor.org/biocLite.R")
# biocLite("org.Dr.eg.db")
library(org.Dr.eg.db)
go = getEnrichedGO(
    annotatedPeaks,
    orgAnn="org.Dr.eg.db",
    multiAdj=TRUE,
    minGOterm=10,
    multiAdjMethod="BH"
    )

ap_withsymbol = addGeneIDs(
    annotatedPeaks,
    "org.Dr.eg.db",c("symbol")
    )

write.table(
    as.data.frame(ap_withsymbol),
    file="~/projects/artinger/data/2Som_chip1_GCCAAT_L006_R1_001_peaks.txt",
    sep="\t",
    row.names=FALSE
    )

table = read.table(
    "~/projects/artinger/data/2Som_chip2_GTCCGC_L006_R1_001_peaks.bed",
    header=FALSE)
bed = BED2RangedData(table)

# annotate peaks
annotatedPeaks = annotatePeakInBatch(bed, AnnotationData=tssanno)
# write.csv(as.data.frame(annotatedPeaks), file="~/remote/projects/artinger/results/common/2Som_chip1_GCCAAT_L006_R1_001/2Som_chip1_GCCAAT_L006_R1_001_peaks.annotated.bed")

# enriched go terms
go = getEnrichedGO(
    annotatedPeaks,
    orgAnn="org.Dr.eg.db",
    multiAdj=TRUE,
    minGOterm=10,
    multiAdjMethod="BH"
)

ap_withsymbol = addGeneIDs(
    annotatedPeaks,
    "org.Dr.eg.db",c("symbol")
)

write.table(
    as.data.frame(ap_withsymbol),
    file="~/projects/artinger/data/2Som_chip2_GTCCGC_L006_R1_001_peaks.txt",
    sep="\t",
    row.names=FALSE
)

table = read.table(
    "~/projects/artinger/data/31hpt_Chip1_CAGATC_L006_R1_001_peaks.bed",
    header=FALSE)
bed = BED2RangedData(table)

# annotate peaks
annotatedPeaks = annotatePeakInBatch(bed, AnnotationData=tssanno)
# write.csv(as.data.frame(annotatedPeaks), file="~/remote/projects/artinger/results/common/2Som_chip1_GCCAAT_L006_R1_001/2Som_chip1_GCCAAT_L006_R1_001_peaks.annotated.bed")

# enriched go terms
go = getEnrichedGO(
    annotatedPeaks,
    orgAnn="org.Dr.eg.db",
    multiAdj=TRUE,
    minGOterm=10,
    multiAdjMethod="BH"
)

ap_withsymbol = addGeneIDs(
    annotatedPeaks,
    "org.Dr.eg.db",c("symbol")
)

write.table(
    as.data.frame(ap_withsymbol),
    file="~/projects/artinger/data/31hpt_Chip1_CAGATC_L006_R1_001_peaks.txt",
    sep="\t",
    row.names=FALSE
)

table = read.table(
    "~/projects/artinger/data/31hpt_Chip2_ACAGTG_L006_R1_001_peaks.bed",
    header=FALSE)
bed = BED2RangedData(table)

# annotate peaks
annotatedPeaks = annotatePeakInBatch(bed, AnnotationData=tssanno)
# write.csv(as.data.frame(annotatedPeaks), file="~/remote/projects/artinger/results/common/2Som_chip1_GCCAAT_L006_R1_001/2Som_chip1_GCCAAT_L006_R1_001_peaks.annotated.bed")

# enriched go terms
go = getEnrichedGO(
    annotatedPeaks,
    orgAnn="org.Dr.eg.db",
    multiAdj=TRUE,
    minGOterm=10,
    multiAdjMethod="BH"
)

ap_withsymbol = addGeneIDs(
    annotatedPeaks,
    "org.Dr.eg.db",c("symbol")
)

write.table(
    as.data.frame(ap_withsymbol),
    file="~/projects/artinger/data/31hpt_Chip2_ACAGTG_L006_R1_001_peaks.txt",
    sep="\t",
    row.names=FALSE
)