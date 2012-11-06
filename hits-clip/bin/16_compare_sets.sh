#! /usr/bin/env bash
#BSUB -J compare_sets
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q normal

<<DOC
Running the requested comparisons between test cases to identify unique miRNA
and the intersection between the cases.
DOC

COMPARISONS=("HS5 vs HS27A"
             "HS5+HS27A vs hMSC"
             "HS5 vs hMSC"
             "HS27A vs hMSC"
             "BMEC vs hMSC"
             "BMEC vs HS27A+HS5"
             "MCF7 vs MCF7estr"
             "BT474 vs BT474estr"
             "MCF7 vs BT474"
             "MCF7 vs MDA231"
             "BT474 vs MDA231"
)

DATA=$HOME/projects/hits-clip/results/common
#EXT=abundance.noa
#EXT=intensity.eda.8.intensity.eda
EXT=7.sif

BIN=$HOME/projects/hits-clip/bin

STATS=network.overlaps.txt

# loop over the array
for ((i=0; i<${#COMPARISONS[@]}; i++)); do
    
    # cut the string by spaces and save each side of the "vs"
    ANAME=$(echo ${COMPARISONS[$i]} | cut -f1 -d" ")
    BNAME=$(echo ${COMPARISONS[$i]} | cut -f3 -d" ")
    
    CASEA=""
    CASEB=""
    
    # account for multiple cases on either side of the "vs" separated by "+"
    for a in $(echo $ANAME | tr "+" "\n"); do
        #CASEA="$CASEA $DATA/$a/$a.$EXT"    # from common
        CASEA="$CASEA $a.$EXT"              # from cwd
    done
    
    for b in $(echo $BNAME | tr "+" "\n"); do
        #CASEB="$CASEB $DATA/$b/$b.$EXT"    # from common
        CASEB="$CASEB $b.$EXT"              # from cwd
    done
    
    # trim whitespace from vars -- completely unnecessary
    CASEA=$(echo $CASEA)
    CASEB=$(echo $CASEB)
    
    python $BIN/compare_sets.py \
            -a $CASEA -aname $ANAME \
            -b $CASEB -bname $BNAME \
            -mode sif >> $STATS
done