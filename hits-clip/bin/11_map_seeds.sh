#! /usr/bin/env bash

function usage () {
    echo "No arguments detected.
         "
    echo "Usage: sh $0 <seed_fastas> <seed_lengths>"
    echo "Example: sh $0 hsa.seeds.2.8.fa hsa.seeds.2.9.fa 7 8"
    exit 1
}

if [ $# -lt 2 ]; then
    usage
fi

ARGS=($@)
SEED_FASTAS=""
SEED_LENGTHS=""
for (( i = 0; i < ${#ARGS[@]}; i++ )); do

    # the first half of the parameters are fastas
    if [ $i -gt $(expr ${#ARGS[@]} / 2 - 1) ]; then
        SEED_LENGTHS="$SEED_LENGTHS ${ARGS[$i]}"
    else
        SEED_FASTAS="$SEED_FASTAS ${ARGS[$i]}"
    fi
done

SAMPLES=(MCF7 MCF7estr MDA231 BT474 BT474estr BT474herc HS27A HS5 hMSC BMEC HUVEC HELA)

GDARCHIVE=$HOME/projects/hits-clip/results/20120917/gd_20120917.rmd
ANNOTATION=$HOME/ref/hg18/refseq.wholegene.bed.gz

# likely gone.
SRC=$HOME/devel/cliptools/cliptools/

for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    SAMPLE=${SAMPLES[$i]}
    DATA="$HOME/projects/hits-clip/results/common/$SAMPLE"

    POSPEAK="$DATA/$SAMPLE.pos.peaks.trimmed.bed.gz"
    NEGPEAK="$DATA/$SAMPLE.neg.peaks.trimmed.bed.gz"

    SIF="$SAMPLE.network.unfiltered.sif"
    FUNCTION="$SAMPLE.function.noa"

    RUN=$SAMPLE.mirnaseeds.sh

    if [ ! -f $SIF ]; then
        echo "#!/usr/bin/env bash" > $RUN
        echo "#BSUB -J $SAMPLE.mirnaseeds" >> $RUN
        echo "#BSUB -o %J.out" >> $RUN
        echo "#BSUB -e %J.err" >> $RUN
        echo "

FASTAS=($SEED_FASTAS)
LENGTHS=($SEED_LENGTHS)

for STRAND in pos neg; do
    for (( i = 0; i < \${#FASTAS[@]}; i++ )); do
        SEED_FASTA=\${FASTAS[\$i]}
        SEED_LENGTH=\${LENGTHS[\$i]}

        POSNOTSIF=\"$SAMPLE.pos.mirna.\$SEED_LENGTH.notSIF.gz\"
        NEGNOTSIF=\"$SAMPLE.neg.mirna.\$SEED_LENGTH.notSIF.gz\"

        if [ \"\$STRAND\" == \"pos\" ]; then
            python $SRC/aggregate_peaks.py $POSPEAK $ANNOTATION \$SEED_FASTA $GDARCHIVE | gzip -c > \$POSNOTSIF
            awk -v len=\"\$SEED_LENGTH\" -c header '{split(\$1,seed,\"|\"); print seed[1]\"\t\"len\"\t\"\$5}' \$POSNOTSIF > $SIF
        else
            python $SRC/aggregate_peaks.py $NEGPEAK $ANNOTATION \$SEED_FASTA $GDARCHIVE | gzip -c > \$NEGNOTSIF
            awk -v len=\"\$SEED_LENGTH\" -c header '{split(\$1,seed,\"|\"); print seed[1]\"\t\"len\"\t\"\$5}' \$NEGNOTSIF >> $SIF
        fi

    done
done

echo \"functionalCategory\" > $FUNCTION
bioawk -c header '{print \$1\" = miRNA\n\"\$3\" = Gene\"}' $SIF | sort | uniq >> $FUNCTION

        " >> $RUN

        bsub < $RUN
    fi
done
