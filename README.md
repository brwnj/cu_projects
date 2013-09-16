#Bennett

+ same variable region or same j region can be considered the same family
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

+ coverage tracks for pooled samples
+ add peter to private projects folder

+ compare ERP to ERN in pooled datasets

+ map to mouse, then to human for peter's previous samples (samples with an 'a')
+ 86-...
+ map to human, map to mouse, what's the overlap

    + map some human to mouse and see how many map (percentage wise)
    + possible to mask human genome from mouse genome based on a certain length

    + testing mouse -- for each polya grab 25, 50, 75, 100 from the site
        + align to mouse and see what maps for each length
        + are there gains in using longer length reads

#Walter

+ levin classification; score
