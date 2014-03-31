# Bennett
+ redo using new primer
+ all pairwise comparisons between cell types for a given subject
+ run through pRESTO

+ meeting with jeff:
    + umi same; collapse into single cdr3
    + r1_primer; build consensus sequence from within UMI bin where UMI+primer ---
    + average sequence
    + plot the UMI; frequency distribution
    + filter 2 stddev; UMI total reads has to be above 2 stddev of observed UMI bin totals

+ overlaps between the csf and blood
+ overlaps between blood
+ align all unique sequences for a patient
+ be able to color by vh family

#Cohrs
+ HOMER pipeline with design

#Duval
+ genes with most mutations across samples
+ non-normal, non-synonymous mutation frequency
+ something like:
    s1  s2  s3  total
g1  1   1   0   2
g2  0   0   1   1
g3  1   1   1   3

#HITS-CLIP
+ new samples through pipeline
+ ucsc hub
+ 1-58
+ 1-50 + 1, 51, 52-57 +1, 58

#Leung
TODO: review data for download; repeating Dan's Cytoscape work

+ annotate non-syn variants everything; place in lookup table for all samples
+ be able to lookup changes for each gene to see the nt and AA effect

# Poly(A)
# Ribosome stuffs
+ read trimming, sequence content
+ linker sequence
+ 27-33bp
+ trim, align

+ analyzing what is different between time points
+ determine what is transcribed vs what is regulated prior to reaching the ribosome
+ set of genes that are changing
+ validate at protein level

#Nicoli
+ annotations for genes sent in email
+ cytoscape networks

#Williams;chipseq
#Tesla
+ set up experiment to test execution time to cpus used
+ test common aligners gsnap, novoalign, etc using 1,2,4,6,8,12 cpus
+ plot results
+ gsnap also has -B and could test 0,1,2,3,4,5
