PROJECTID=pillai_kabos_polya
SAMPLES=(MP51 MP52 MP53 MP54 MP55 MP56 MP57 MP58 MP59 MP60              #10
        MP61 MP62 MP63 PK61 PK62 PK63 PK64 PK65 PK66 PK67               #20
        PK68 PK69 PK70 PK71 PK72 PK73 PK74 PK75 PK76 PK77               #30
        PK78 PK79 PK80 PK81 PK82 PK83 PK84 PK85 PK86 PK87               #40
        PK88 PK89 PK90 PK91 PK92 MP70 MP71 MP72 MP73 MP74               #50
        MP75 PK93 PK94 PK95 PK96 PK97 PK98 PK99 PK100 PK101             #60
        PK102 PK103 PK104 PK105 PK106 PK107 PK108 PK109 PK110)          #69
NOVOIDX=$HOME/ref/hg19/hg19.9606.novoidx
RESULT=$HOME/projects/polya/results/common
DATA=$HOME/projects/polya/data/common
CHROM_SIZES=$HOME/ref/hg19/hg19.sizes
BIN=$HOME/devel/polya
UMI=NNNNNV
FASTA=$HOME/ref/hg19/hg19.fa
EXONS=$BIN/ref/refseq.exons.hg19.trimmedname.bed.gz
XREF=$BIN/ref/refseq.xref.txt.gz
WHOLEGENE=$BIN/ref/refseq.wholegene.hg19.trimmedname.bed.gz
SITESHIFTS=$RESULT/site_shifts
POLYASITES=$RESULT/polya_sites
HUB=$RESULT/hub
FISHERRESULTS=$RESULT/fisher_results
METADATA=$HUB/metadata.tsv
POOLEDRESULTS=$RESULT/pooled_results
