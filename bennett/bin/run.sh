#!/usr/bin/env bash
#BSUB -J igg_repertoire[1-16]
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
jobids="1"

# prepend the umi back onto the reads
i1=$READS/${sample}_I1.fastq.gz
r1_in=$READS/${sample}_R1.fastq.gz
r2_in=$READS/${sample}_R2.fastq.gz
r1_prependumi=${r1_in/.fastq.gz/.umi.fastq.gz}
r2_prependumi=${r2_in/.fastq.gz/.umi.fastq.gz}
r1_sortumi=${r1_in/.fastq.gz/.umi.sorted.fastq.gz}
r2_sortumi=${r2_in/.fastq.gz/.umi.sorted.fastq.gz}
r1_umifiltered=${r1_in/.fastq.gz/.umi.sorted.umifiltered.fastq.gz}
r2_umifiltered=${r2_in/.fastq.gz/.umi.sorted.umifiltered.fastq.gz}
r1_rmadptr=${r1_in/.fastq.gz/.umi.sorted.umifiltered.rmadptr.fastq.gz}
r2_rmadptr=${r2_in/.fastq.gz/.umi.sorted.umifiltered.rmadptr.fastq.gz}
joined=$READS/${sample}.joined.fastq.gz

if [[ ! -f $r1_prependumi ]]; then
    job=$(echo "python $bin/prepend_umi.py -b 6 -e 14 $i1 $r1_in | gzip -c > $r1_prependumi" | bsez prepend_umi -P $PI)
    jobid=$(echo $job | sed -rn 's/.*<([0-9]+)>.*/\1/p')
    jobids="$jobids $jobid"
fi
if [[ ! -f $r2_prependumi ]]; then
    job=$(echo "python $bin/prepend_umi.py -b 6 -e 14 $i1 $r2_in | gzip -c > $r2_prependumi" | bsez prepend_umi -P $PI)
    jobid=$(echo $job | sed -rn 's/.*<([0-9]+)>.*/\1/p')
    jobids="$jobids $jobid"
fi

python -m bsub $jobids

# sort the reads according to umi sequence
if [[ ! -f $r1_sortumi ]]; then
    job=$(echo "python $bin/process_umi.py sort $r1_prependumi $UMILENGTH | gzip -c > $r1_sortumi" | bsez sort_umi -P $PI)
    jobid=$(echo $job | sed -rn 's/.*<([0-9]+)>.*/\1/p')
    jobids="$jobids $jobid"
fi
if [[ ! -f $r2_sortumi ]]; then
    job=$(echo "python $bin/process_umi.py sort $r2_prependumi $UMILENGTH | gzip -c > $r2_sortumi" | bsez sort_umi -P $PI)
    jobid=$(echo $job | sed -rn 's/.*<([0-9]+)>.*/\1/p')
    jobids="$jobids $jobid"
fi

python -m bsub $jobids

# collapse and filter out redundant sequences per umi
if [[ ! -f $r1_umifiltered ]] || [[ ! -f $r2_umifiltered ]]; then
    python $bin/process_umi.py collapse $r1_sortumi $r2_sortumi $UMISEQ
fi

# trim the adapters and append known primer to read name
if [[ ! -f $r1_rmadptr ]] || [[ ! -f $r2_rmadptr ]]; then
    python $BIN/trim_adapters.py $r1_umifiltered $r2_umifiltered $R1PRIMERS $R2PRIMERS
fi

# join r1 and r2 based on local alignment of the overlap
if [[ ! -f $joined ]]; then
    echo "python $bin/join_reads.py -t 8 $r1_rmadptr $r2_rmadptr | gzip -c > $joined" | bsez join_reads -P $PI -n 8 -R span[hosts=1]
fi
