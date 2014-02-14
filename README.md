# Bennett
+ viewing in a tree may be helpful for manual picking
+ 93 plasmablast, align to 9-3 A and B populations

+ proteimoics sheet, align to whole IMGT resultant seq with fewer mismatches, write table
+ compare within a patient or across all depending on time to process
+ only compare CDR3 within length +- 2

+ overlaps between clones
+ overlaps of CDR3 sequences of all vs all

+ stats for v, d, and j families so they know the distribution of the composition
+ this has to do with vh families and what not
+ tables with dists, they'll run any stats
+ individual sequences that are duplicated many times (would be across )

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
+ UCSC session with alignments loaded
+ ingenuity for ensembl gene name to gene symbol along with any additional details
+ coverage measures for each sample across capture intervals
+ along with coverage, count the number of (not synonymous) variants
+ set up the script to filter variant count based on coverage of individual mutations
+ homozygousity rate over mutation rate
+ LD block to identify regions of mass deletions

#HITS-CLIP
+ new samples through pipeline
+ ucsc hub
+ 6 hour, remap using iterative trimming; they'll need tracks to see whether ER reads come back

#Leung
+ call variants on everything

# Poly(A)

### Mouse
+ map to mouse, then to human for peter's previous samples (samples with an 'a')
+ 86-...
+ map to human, map to mouse, what's the overlap
    + map some human to mouse and see how many map (percentage wise)
    + possible to mask human genome from mouse genome based on a certain length
    + testing mouse -- for each polya grab 25, 50, 75, 100 from the site
        + align to mouse and see what maps for each length
        + are there gains in using longer length reads


# Ribosome stuffs

+ read trimming, sequence content
+ linker sequence
+ 27-33

+ tables, genomedata, rpy2
+ trim, align

#Nicoli
+ liftover data to hg19
+ continue pipeline
+ take seed sequences from miRBase
+ map back to the mRNA enriched regions to come up with the linkages and miRNA quantification

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

#Projects
+ williams chipseq
+ cohrs chipseq
+ duval canine exome
+ kabos polya
+ kabos hitsclip
+ kabos ribosome profiling
+ kabos lincRNA expression profiling
+ bennett ig repertoire in neuromyelitis optica (NMO)
+ bennett ig repertoire in multiple sclerosis (MS)
+ nicoli miRNA quantification and networks of interactions

