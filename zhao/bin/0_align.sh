#!/usr/bin/env bash
#BSUB -J align[1-3]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 8
#BSUB -P zhao

<<DOC
	                    ACATCG	GCCTAA	TGGTCA
Paired Reads	        2469370	3748472	3050090
Pairs Aligned	        2273604	3440160	2866330
Read Sequences	        4938740	7496944	6100180
Aligned	                4600208	7014260	5788755
Unique Alignment	    3947426	6031922	5010330
Gapped Alignment	    221906	263355	186405
Quality Filter	        0	    0	    0
Homopolymer Filter	    38	    70	    55
Perc Aligned	        93.14	93.56	94.89
Perc Uniquely Aligned	79.92	80.45	82.13
DOC


set -o nounset -o pipefail -o errexit -x
source $HOME/projects/zhao/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
r1=$DATA/${sample}_R1.fastq.gz
r2=$DATA/${sample}_R2.fastq.gz
results=$RESULTS/$sample
bam=$results/$sample.bam
if [[ ! -d $results ]]; then
    mkdir -p $results
fi
if [[ ! -f $bam ]]; then
    novoalign -d $NOVOIDX -f $r1 $r2 -o SAM "@RG\tID:$sample\tSM:$sample\tPU:Illumina\tLB:PE" -r None -i 250 100 -c 8 -k \
        2> $sample.alignment.txt \
        | samtools view -ShuF4 - \
        | samtools sort -o -m 8G - $sample.temp \
        > $bam
fi
