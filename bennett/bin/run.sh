#!/usr/bin/env bash
#BSUB -J igg_repertoire[1-64]
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

# prepend the umi back onto the reads
i1=$READS/${sample}_I1.fastq.gz
r1_in=$READS/${sample}_R1.fastq.gz
r2_in=$READS/${sample}_R2.fastq.gz
r1_prependumi=${r1_in/.fastq.gz/.umi.fastq.gz}
r2_prependumi=${r2_in/.fastq.gz/.umi.fastq.gz}
r1_sortumi=${r1_in/.fastq.gz/.umisorted.fastq.gz}
r2_sortumi=${r2_in/.fastq.gz/.umisorted.fastq.gz}
r1_umifiltered=${r1_in/.fastq.gz/.umifiltered.fastq.gz}
r2_umifiltered=${r2_in/.fastq.gz/.umifiltered.fastq.gz}
r1_trimadapter=${r1_in/.fastq.gz/.trimmed.fastq.gz}
r2_trimadapter=${r2_in/.fastq.gz/.trimmed.fastq.gz}
joined=$READS/${sample}.joined.fastq.gz

out=prepend_umi.$LSB_JOBID.$LSB_JOBINDEX.out
err=prepend_umi.$LSB_JOBID.$LSB_JOBINDEX.err
if [[ ! -f $r1_prependumi ]]; then
    # check for samples from Gao...
    if [[ condition ]]; then
        bsub -J prepend_umi -o $out -e $err -P $PI -K "python $bin/prepend_umi.py $i1 $r1_in | gzip -c > $r1_prependumi" &
    else
        bsub -J prepend_umi -o $out -e $err -P $PI -K "python $bin/prepend_umi.py -b 6 -e 14 $i1 $r1_in | gzip -c > $r1_prependumi" &
    fi
fi

if [[ ! -f $r2_prependumi ]]; then
    if [[ condition ]]; then
        bsub -J prepend_umi -o $out -e $err -P $PI -K "python $bin/prepend_umi.py $i1 $r2_in | gzip -c > $r2_prependumi" &
    else
        bsub -J prepend_umi -o $out -e $err -P $PI -K "python $bin/prepend_umi.py -b 6 -e 14 $i1 $r2_in | gzip -c > $r2_prependumi" &
    fi
fi
wait

# sort the reads according to umi sequence
out=sort_umi.$LSB_JOBID.$LSB_JOBINDEX.out
err=sort_umi.$LSB_JOBID.$LSB_JOBINDEX.err
if [[ ! -f $r1_sortumi ]]; then
    bsub -J sort_umi -o $out -e $err -P $PI -K "python $bin/process_umi.py sort $r1_prependumi $UMILENGTH | gzip -c > $r1_sortumi" &
fi

if [[ ! -f $r2_sortumi ]]; then
    bsub -J sort_umi -o $out -e $err -P $PI -K "python $bin/process_umi.py sort $r2_prependumi $UMILENGTH | gzip -c > $r2_sortumi" &
fi
wait

# collapse and filter out redundant sequences per umi
if [[ ! -f $r1_umifiltered ]] || [[ ! -f $r2_umifiltered ]]; then
    python $bin/process_umi.py collapse $r1_sortumi $r2_sortumi $UMISEQ
fi

# trim the adapters and append known primer to read name
if [[ ! -f $r1_trimadapter ]] || [[ ! -f $r2_trimadapter ]]; then
    python $bin/trim_adapters.py --mismatches 3 $r1_umifiltered $r2_umifiltered $R1PRIMERS $R2PRIMERS
fi

# join r1 and r2 based on local alignment of the overlap
if [[ ! -f $joined ]]; then
    echo "python $bin/join_reads.py -t 8 $r1_trimadapter $r2_trimadapter | gzip -c > $joined" | bsez join_reads -P $PI -n 8 -R span[hosts=1]
fi
