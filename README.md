# Dev

+ script to read bams, output normalized counts table

#Bennett

+ same variable region or same j region can eb considered the same family
+ would likely have changes in the cdr3 while being close to match with very similar variable region or j region
+ report dimer (short mer) per barcode across runs
+ counts of reads per barcode
+ counts of reads after umi filtering per barcode
+ counts of reads after joining per barcode
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

+ human islet rnaseq data
+ converet database into ORFs
+ ORF possibilities; Howard may follow up
+ 8-10 mer might be coming from different ORFs
+ 25 might be a better cutoff for an ORF
+ all possible ORFs at 25mer or another length
+ also need info on where the read was aligned; like gene or isoform

+ nonredundant list of ORFs
+ matching protein sequence to protein sequence

+ swissprot fasta header for parsing
+ translate known and novel protein seqs in known ORFs

# Duval

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

## tracks for PK61, PK62, PK63

+ still need to run PK61-63
+ heatmap of miRNA, something that you can zoom into
+ matrix2png

#Poly(A)

+ move to hg19

+ currently picking up genes where expression is similar and shift is observed
+ cpox is an example across manoj's data
+ normalize per gene
+ get the gene count; calculate scale factor; divide each polya site by the scale factor
+ fisher test at that polya site, output something similar to dexseq; normalize by total count
+ qvalue calculation
+ no need for fold change, but insert the values of the two samples being compared

+ tables for dexseq results; sites with pval, padj, basemean

+ map to mouse, then to human for peter's previous samples (samples with an 'a')
    + map some human to mouse and see how many map (percentage wise)
    + possible to mask human genome from mouse genome based on a certain length

    + testing mouse -- for each polya grab 25, 50, 75, 100 from the site
        + align to mouse and see what maps for each length
        + are there gains in using longer length reads

#Walter

+ levin classification; score
+ given list with direction, compute a score for each patient along with their group (TB, noinf, etc.)
+ calculate threshold for plot

+ adding total intensity at up-regulated transcripts and subtracting total intensity at all down-regulated transcripts

#Zhao

+ deliver mutation list
