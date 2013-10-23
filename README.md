#Bennett

+ cutoff for umi to only take the most abundant sequences accounting for the majority of what's present
    + to be implemented after initial assessment scripts are written

+ group by v-gene, then j-gene, then d-gene
    + then list unique cdr3 sequences
    + for groups of cdr3 sequences, give ability to group similar based on number of mismatches
+ in the output, print simple grouped table; stats can later be done on this table

+ would likely have changes in the cdr3 while being close to match with very similar variable region or j region
+ report dimer (short mer) per barcode across runs
+ counts of translated productive reads per barcode
+ counts of translated unique productive reads per barcode

+ same variable and j segment being used
+ going back to blood likely be able to pick up where the change occurs

+ same v region and j segment family while having...
+ up to 3 mismatches in cdr3

+ viewing in a tree may be helpful for manual picking
+ 53 different v
+ 6 j
+ barcode to cell type conversion table

+ match protein sequences to imgt results
+ design best output format for this task

#Davidson

# Duval

+ polyphen
+ indels mostly, annotate with descriptions if available
+ EXOME

+ canine exome; bladder cancer
+ experimental design

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

+ chromosome M: is it in the reference, is it in the exon model, do we have classified peaks there, are they just not significant changes
+ pooled counts from bedgraphs
+ automate the pooled processing portion
+ update the hub with the new comparisons

+ exome reference for peak classification
    + convert current coord to exon coord and pull sequence from exome
    + building exome reference and handling different transcripts

+ map to mouse, then to human for peter's previous samples (samples with an 'a')
+ 86-...
+ map to human, map to mouse, what's the overlap

    + map some human to mouse and see how many map (percentage wise)
    + possible to mask human genome from mouse genome based on a certain length

    + testing mouse -- for each polya grab 25, 50, 75, 100 from the site
        + align to mouse and see what maps for each length
        + are there gains in using longer length reads
