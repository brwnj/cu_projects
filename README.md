#Bennett

+ cutoff for umi to only take the most abundant sequences accounting for the majority of what's present
+ to be implemented after initial assessment scripts are written
+ viewing in a tree may be helpful for manual picking
+ barcode to cell type conversion table

#Davidson

# Duval

+ polyphen
+ indels mostly, annotate with descriptions if available
+ EXOME
+ canine exome; bladder cancer
+ experimental design

###Expecting a list for this portion
+ chipseq; human; run through homer pipeline; will receive sample definitions at some
+ 900 miRNAs NCI60 something
+ map our numbers and match names to those existing
+ accession number of geo - normalized binned miRs to compare to other geo sets
+ list of 261 from the drug resistance testing
+ take a closer look at the park sample. it has many more positive counts than the others
+ reduce the number to nci60 and mirbase. annotate; blast; etc
+ gse26375
+ email: should receive list of oncomiRs that should also
+ annotate to mirbase; should be seeing mirbase names being counted

#Hits-Clip

+ tracks for PK61, PK62, PK63
+ heatmap of miRNA, something that you can zoom into
+ matrix2png

#Poly(A)

###20131024
+ chrM UCSC gene model, fisher testing on the combined set (only chrM genes) across samples
+ maybe put that track up on UCSC along with Peter's others
+ also interested to see peaks called in chrM among the HITS-CLIP samples
+ possibly build comprehensive trackhub to later be expanded with the HITS-CLIP data

###20131017
+ chromosome M: is it in the reference, is it in the exon model, do we have classified peaks there, are they just not significant changes
+ pooled counts from bedgraphs
+ automate the pooled processing portion
+ update the hub with the new comparisons
+ exome reference for peak classification
    + convert current coord to exon coord and pull sequence from exome
    + building exome reference and handling different transcripts

###Mouse
+ map to mouse, then to human for peter's previous samples (samples with an 'a')
+ 86-...
+ map to human, map to mouse, what's the overlap
    + map some human to mouse and see how many map (percentage wise)
    + possible to mask human genome from mouse genome based on a certain length
    + testing mouse -- for each polya grab 25, 50, 75, 100 from the site
        + align to mouse and see what maps for each length
        + are there gains in using longer length reads
