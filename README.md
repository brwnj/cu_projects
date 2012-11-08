==============================================================================
Duval
==============================================================================

* Counts and DESeq

==============================================================================
Hits-Clip
==============================================================================

* stats for snoRNA
* ayb with shorter and shorter (less priority)

==============================================================================
Marrack
==============================================================================

    ==> RS_input_CCGTCC_L005_R1_001/mapping_stats.txt <==
    genome size: 2,725,765,481
    number of bases covered by unique mappers: 32,581,696 (1.19%)
    number of bases covered by non-unique mappers: 5,794,538 (0.21%)
    Number of reads: 7,856,284
    ------
    UNIQUE MAPPERS: 3,757,338 (47.82%)
    NON-UNIQUE MAPPERS: 222,422 (2.8%)
    -----
    TOTAL: 3,979,760 (50.6%)

    ==> RS_iso_ATGTCA_L005_R1_001/mapping_stats.txt <==
    genome size: 2,725,765,481
    number of bases covered by unique mappers: 3,404,077 (0.12%)
    number of bases covered by non-unique mappers: 1,102,603 (0.04%)
    Number of reads: 9,994,627
    ------
    UNIQUE MAPPERS: 4,225,504 (42.27%)
    NON-UNIQUE MAPPERS: 178,947 (1.7%)
    -----
    TOTAL: 4,404,451 (44%)

    ==> RS_tbet_CTTGTA_L005_R1_001/mapping_stats.txt <==
    genome size: 2,725,765,481
    number of bases covered by unique mappers: 3,544,680 (0.13%)
    number of bases covered by non-unique mappers: 923,512 (0.03%)
    Number of reads: 4,848,946
    ------
    UNIQUE MAPPERS: 1,252,496 (25.83%)
    NON-UNIQUE MAPPERS: 57,767 (1.1%)
    -----
    TOTAL: 1,310,263 (27%)

The ISO is an off-target IP control and Tbet is the target IP that you want to examine
input sample as a random background
then use the ISO sample to subtract any non-specific pull-downs that are seen in both ISO and Tbet

* Peaks
    * control background = RS_input_CCGTCC_L005_R1_001
* filter peaks by subtracting iso from tbet
* deliver tbet only peaks

==============================================================================
Poly(A)
==============================================================================

    samtools view bam | python dexseq_count.py flattened.gff - out.counts
    awk 'BEGIN{OFS=FS="\t"}{if($2!=0){n=int(rand()*5);print $1,$2+n}else{print $0}}' out.counts > outx.counts

* neg and pos split everything; annotation, bam, etc
* tracks have to be shown only using 5' read ends
* re-run dexseq