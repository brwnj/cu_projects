#!/usr/bin/env bash
#BSUB -J igg_repertoire[1-18]
#BSUB -e igg_repertoire.%J.%I.err
#BSUB -o igg_repertoire.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P bennett

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/bennett/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
bin=$HOME/projects/bennett/bin
jobids=""

# prepend the umi back onto the reads
i1=$READS/${sample}_I1.fastq.gz
r1_prependumi_in=$READS/${sample}_R1.fastq.gz
r2_prependumi_in=$READS/${sample}_R2.fastq.gz
r1_prependumi_out=$READS/${sample}_R1.umi.fastq.gz
r2_prependumi_out=$READS/${sample}_R2.umi.fastq.gz

if [[ ! -f $r1_prependumi_out ]]; then
    job=$(echo "python $bin/prepend_umi.py -b 6 -e 14 $i1 $r1_prependumi_in | gzip -c > $r1_prependumi_out" | bsez prepend_umi -P $PI)
    jobid=$(echo $job | sed -rn 's/.*<([0-9]+)>.*/\1/p')
    jobids="$jobids $jobid"
fi
if [[ ! -f $r2_prependumi_out ]]; then
    job=$(echo "python $bin/prepend_umi.py -b 6 -e 14 $i1 $r2_prependumi_in | gzip -c > $r2_prependumi_out" | bsez prepend_umi -P $PI)
    jobid=$(echo $job | sed -rn 's/.*<([0-9]+)>.*/\1/p')
    jobids="$jobids $jobid"
fi

python -m bsub $jobids

# sort the reads according to umi sequence
r1_sortumi_out=$READS/${sample}_R1.umi.sorted.fastq.gz
r2_sortumi_out=$READS/${sample}_R2.umi.sorted.fastq.gz

if [[ ! -f $r1_sortumi_out ]]; then
    job=$(echo "python $bin/process_umi.py sort $r1_prependumi_out $UMILENGTH | gzip -c > $r1_sortumi_out" | bsez sort_umi -P $PI)
    jobid=$(echo $job | sed -rn 's/.*<([0-9]+)>.*/\1/p')
    jobids="$jobids $jobid"
fi
if [[ ! -f $r2_sortumi_out ]]; then
    job=$(echo "python $bin/process_umi.py sort $r2_prependumi_out $UMILENGTH | gzip -c > $r2_sortumi_out" | bsez sort_umi -P $PI)
    jobid=$(echo $job | sed -rn 's/.*<([0-9]+)>.*/\1/p')
    jobids="$jobids $jobid"
fi

python -m bsub $jobids

# collapse and filter out redundant sequences per umi
r1_umifiltered_out=$READS/${sample}_R1.umi.sorted.umifiltered.fastq.gz
r2_umifiltered_out=$READS/${sample}_R2.umi.sorted.umifiltered.fastq.gz

if [[ ! -f $r1_umifiltered_out ]] || [[ ! -f $r2_umifiltered_out ]]; then
    python $bin/process_umi.py collapse $r1_sortumi_out $r2_sortumi_out $UMISEQ
fi

# trim the adapters and append known primer to read name
r1_rmadptr_out=$READS/${sample}_R1.umi.sorted.umifiltered.rmadptr.fastq.gz
r2_rmadptr_out=$READS/${sample}_R2.umi.sorted.umifiltered.rmadptr.fastq.gz

if [[ ! -f $r1_rmadptr_out ]] || [[ ! -f $r2_rmadptr_out ]]; then
    python $BIN/trim_adapters.py $r1_umifiltered_out $r2_umifiltered_out $R1PRIMERS $R2PRIMERS
fi

# join r1 and r2 based on local alignment of the overlap
join_out=$READS/${sample}.joined.fastq.gz

if [[ ! -f $join_out ]]; then
    echo "python $bin/join_reads.py -t 8 $r1_rmadptr_out $r2_rmadptr_out | gzip -c > $join_out" | bsez join_reads -P $PI -n 8 -R span[hosts=1]
fi
