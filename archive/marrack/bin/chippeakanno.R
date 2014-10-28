library(ChIPpeakAnno)
table = read.table("~/projects/marrack/data/RS_tbet_CTTGTA_L005_R1_001_unique_peaks.clipped.bed", header=F)
df = data.frame(peaksbed)
bed = BED2RangedData(peaks)

# annotation
mart = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
tssanno = getAnnotation(mart, "TSS")

# annotate peaks
annotatedPeaks = annotatePeakInBatch(bed, AnnotationData=tssanno)
write.csv(as.data.frame(annotatedPeaks), file="~/projects/marrack/data/tbet_annotated_peaks.csv")

# enriched go terms
library(org.Mm.eg.db)
go = getEnrichedGO(annotatedPeaks,
                   orgAnn="org.Mm.eg.db",
                   multiAdj=T,
                   minGOterm=10,
                   multiAdjMethod="BH")
ap_withsymbol = addGeneIDs(annotatedPeaks, "org.Mm.eg.db",c("symbol"))
write.csv(as.data.frame(ap_withsymbol), file="~/projects/marrack/data/tbet_annotated_peaks.csv")