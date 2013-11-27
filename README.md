# Bennett
+ cutoff for umi to only take the most abundant sequences accounting for the majority of what's present
+ viewing in a tree may be helpful for manual picking
+ barcode to cell type conversion table
+ markus really liked the local alignment results more so than the current way
+ figure out why CG1 is missing the match shown in 9/16 email
    + rerun alignments using local aligner and give those results

#Cohrs
+ design saved in starred email
+ check the QC, should need to trim back to 220bp
+ also quality trim with seqtk

# Duval
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

# Hits-Clip
+ ucsc hub

#Leung
+ map the data i have transferred
+ gatk up to realigned bam
http://seqanswers.com/wiki/How-to/exome_analysis

# Poly(A)
+ were peaks called in the chrM in the hits-clip data? peter would like to know
+ chrM UCSC gene model, fisher testing on the combined set (only chrM genes) across samples
+ add UCSC chrM to refseq model and re-run pipeline from step 4.
+ refresh track hub

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
+ snpeff
+ deliver variants within:
    + Murderball  chr7:21000000-35000000
    + Nomo        chr15:95000000-104000000