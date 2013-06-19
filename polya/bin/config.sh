PROJECTID=pillai_kabos_polya
SAMPLES=(MP51 MP52 MP53 PK61 PK62 PK100-3 PK86-7 PK93-1 PK93-3 PKm+c    #10
    PKm+e2 PKnbt29 PKnbt39 PKnbt89 PKts21 PKts28 PKts57 MP54 MP55 MP56  #20
    MP57 MP58 MP59 MP60 MP61 MP62 MP63 PK93-1a PK86-7a PK93-3a          #30
    PK100-3a PK89-5a PK89-9a PKNbt94 PKTs61 PKNbt89 PKTs57 PKNbt102 PKTs68) #39
NOVOIDX=$HOME/projects/hits-clip/data/common/novoalign/hg18
RESULT=$HOME/projects/polya/results/common
DATA=$HOME/projects/polya/data/common
CHROM_SIZES=$HOME/ref/hg18/hg18.sizes
BIN=$HOME/devel/polya
# REFBED=$DATA/polya_3utr_slop.bed
UMI=NNNNNV
FASTA=$HOME/ref/hg18/hg18.fa
EXONS=$HOME/ref/hg18/refseq.exon.bed.gz
RUNDEXSEQ=$BIN/run_dexseq.R
