# Bennett
+ viewing in a tree may be helpful for manual picking
+ 93 plasmablast, align to 9-3 A and B populations

+ proteimoics sheet, align to whole IMGT resultant seq with fewer mismatches, write table
+ compare within a patient or across all depending on time to process
+ only compare CDR3 within length +- 2

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
+ liftover to hg19
+ ucsc hub

$ python ../../bin/scripts/trim_adapter.py --reverse-complement -a CCGCTGGAAGTGACTGACAC MCF10A.untrimmed.fastq.gz | gzip -c > MCF10A.fastq.gz
Total reads: 32976107
Passing reads: 8056294
Average passing read length: 26.184786578
Failed reads due to length: 9701918
Average failed read length: 1.78702066952
Failed reads due to missing seed sequence: 15217895

$ python ../../bin/scripts/trim_adapter.py --reverse-complement -a CCGCTGGAAGTGACTGACAC PK12-0.untrimmed.fastq.gz | gzip -c > PK12-0.fastq.gz
Total reads: 44194116
Passing reads: 13935940
Average passing read length: 23.0765316871
Failed reads due to length: 8077624
Average failed read length: 0.911654714307
Failed reads due to missing seed sequence: 22180552

$ python ../../bin/scripts/trim_adapter.py --reverse-complement -a CCGCTGGAAGTGACTGACAC PK12-6.untrimmed.fastq.gz | gzip -c > PK12-6.fastq.gz
Total reads: 28544477
Passing reads: 8839771
Average passing read length: 23.1146355488
Failed reads due to length: 4303307
Average failed read length: 1.0208658132
Failed reads due to missing seed sequence: 15401399

$ python ../../bin/scripts/trim_adapter.py --reverse-complement -a CCGCTGGAAGTGACTGACAC PK12-24.untrimmed.fastq.gz | gzip -c > PK12-24.fastq.gz
Total reads: 39612559
Passing reads: 10884219
Average passing read length: 23.3281243238
Failed reads due to length: 5658198
Average failed read length: 1.0089576222
Failed reads due to missing seed sequence: 23070142

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


#Nicoli
+ liftover data to hg19
+ continue pipeline
+ take seed sequences from miRBase
+ map back to the mRNA enriched regions to come up with the linkages and miRNA quantification

#Williams;chipseq
+ do things

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

