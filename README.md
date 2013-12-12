# Bennett
+ umi cutoff now takes top matches down to cutoff of 5
+ viewing in a tree may be helpful for manual picking
+ barcode to cell type conversion table!
+ markus really liked the local alignment results more so than the current way
+ figure out why CG1 is missing the match shown in 9/16 email
    + rerun alignments using local aligner and give those results

+ sequencing in gao's core rather than pollock

#Cohrs
+ design saved in starred email
+ check the QC, should need to trim back to 220bp
+ also quality trim with seqtk

#Duval
+ mapping stats and fastqc results if specified
+ deliver freebayes variants
+ overlap tumor and normal, take unique
+ filter by coverage, mutation type

Tumor   Normal
500     172160
353     353M2
400     166408

Sample 353 is tagged GGCTAC
Sample 353M2 is tagged GAAACC

#HITS-CLIP
+ liftover to hg19
+ ucsc hub

#Leung
+ map the data i have transferred
+ gatk up to realigned bam
http://seqanswers.com/wiki/How-to/exome_analysis

# Poly(A)
+ were peaks called in the chrM in the hits-clip data? peter would like to know
+ chrM UCSC gene model, fisher testing on the combined set (only chrM genes) across samples

### Mouse
+ map to mouse, then to human for peter's previous samples (samples with an 'a')
+ 86-...
+ map to human, map to mouse, what's the overlap
    + map some human to mouse and see how many map (percentage wise)
    + possible to mask human genome from mouse genome based on a certain length
    + testing mouse -- for each polya grab 25, 50, 75, 100 from the site
        + align to mouse and see what maps for each length
        + are there gains in using longer length reads


#Hanson
+ delivered

#Nicoli
+ liftover data to hg19

#sequence_clustering
+ filtering rRNA
+ we're seing 5' noise prior to the mIR sequence
+ the contaminant type will obviously depend on the library's composition
+ how do we know when the bin has crap on the edge of it
+ ends could be actually 3' or 5' real mutations
+ if something matches ribosomal RNA; continue
+ if illumina adapter is present, trim it
+ blat all sequences; prefiltering to rRNA database
+ align all sequences to genome, annotate to mapped position
+ annotate upstream and downstream to nearest gene as well
+ soft clipping will help trim those illumina adapters
+ vervet monkey just wants to know what the sequence is; annotation
+ crystal's project wants to know everything, including the garbage things
+ clean up counts by sno or mir or whatever
+
