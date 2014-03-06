# Bennett
+ stats for v, d, and j families so they know the distribution of the composition
+ this has to do with vh families and what not
+ tables with dists, they'll run any stats
+ individual sequences that are duplicated many times (would be across)

+ same v and j region, unlikely to have different CDR3 sequence
+ need to have another column in the metadata counting those unique sequences in this manner
+ don't care about the actual CDR3 sequence at this point
+ need population dist for these sequences as well

+ MS and NMO see higher perc. within vh2 and vh4

+ overlapping > binning > alignments

+ overlaps between the csf and blood
+ tree needs to be colored in order to identify the transition; to maybe id the origin of the sequence


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

#Leung
+ annotate variants with snpsift

# Poly(A)

# Ribosome stuffs
+ read trimming, sequence content
+ linker sequence
+ 27-33bp
+ trim, align

#Nicoli
+ annotations for genes sent in email
+ cytoscape networks

#Williams;chipseq
+ peaks that don't change between normal and mitotic
+ differential peaks of mitotics
+ take ones that are up in mitotic
+ subtract differentials from everything else
+ ratio of peak size between conditions; conditions diff bound would be significant;

+ interested in peaks that remain the same between mitotic and non-mitotic conditions
+ interested in peaks that are differentially bound between conditions

#Tesla
+ set up experiment to test execution time to cpus used
+ test common aligners gsnap, novoalign, etc using 1,2,4,6,8,12 cpus
+ plot results
+ gsnap also has -B and could test 0,1,2,3,4,5
