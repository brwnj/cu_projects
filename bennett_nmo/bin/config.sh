#!/usr/bin/env bash
set -o nounset

SAMPLES=(
control_ACATCG
control_ATTGGC
control_CACTGT
control_CGTGAT
control_GATCTG
control_GCCTAA
control_TCAAGT
control_TGGTCA
MS13-02_gc1blood
MS13-02_gc1csf
MS13-02_memoryblood
MS13-02_memorycsf
MS13-02_naiveblood
MS13-02_naivecsf
MS13-02_plasmablastblood
MS13-02_plasmablastcsf
ON07-05A_gc1
ON07-05A_memorya
ON07-05A_memoryb
ON07-05A_naive
ON07-05A_plasmablasta
ON07-05A_plasmablastb
ON07-05B_memory
ON07-05B_plasmablasts
ON08-08_gc
ON08-08_memory
ON08-08_naive
ON08-08_plasmablast
ON09-03A_gc1
ON09-03A_gc2
ON09-03A_memory
ON09-03A_naive
ON09-03A_plasmablast
ON09-03B_gc1
ON09-03B_gc2
ON09-03B_memory
ON09-03B_naive
ON09-03B_plasmablast
ON10-01A_gc1
ON10-01A_gc2
ON10-01A_memory
ON10-01A_naive
ON10-01A_plasmablast
ON10-01B_gc1
ON10-01B_gc2
ON10-01B_memory
ON10-01B_naive
ON10-01B_plasmablast
ON10-03A_gc1
ON10-03A_gc2
ON10-03A_memory
ON10-03A_naive
ON10-03A_plasmablast
ON10-03B_gc1
ON10-03B_gc2
ON10-03B_memory
ON10-03B_naive
ON10-03B_plasmablast
ON11-04_gc
ON11-04_memory
ON11-04_naive
ON11-04_plasmablasts
TUM09-527_gc
TUM09-527_plasmablast
MS13-02B_NaiveBlood
MS13-02B_GCBlood
MS13-02B_MemoryBlood
MS13-02B_PlasmablastBlood
MS13-02B_NaiveCSF
MS13-02B_GCCSF
MS13-02B_MemoryCSF
MS13-02B_PlasmablastCSF
)

PI=bennett
RESULTS=$HOME/projects/bennett/results/common
READS=$HOME/projects/bennett/data/common
BIN=$HOME/devel/repertoire/repertoire
R1PRIMERS=$READS/r1_primers.fasta
R2PRIMERS=$READS/r2_primers.fasta
UMILENGTH=8
UMISEQ=NNNNNNNN